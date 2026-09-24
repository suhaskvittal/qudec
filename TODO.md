# TODO

Working notes and next steps, written for an agent (or a future me) picking
this up cold. Last updated 2026-09-23.

Read `docs/runahead-decoding.md` first — it is the framing everything here
sits inside.

---

## Where things stand

Two threads ran this session: reading the accuracy-decoder literature, and
getting Tesseract usable inside qudec.

### Papers read and written up

All PDFs plus plain-text extractions are in `docs/papers/`, summarized in
`docs/papers/NOTES.md`:

- **Libra** (arXiv:2408.12135, Cody Jones) — ensemble decoder that *splices*
  solutions rather than picking one. Core trick is "matching synthesis":
  symmetric-difference two matchings, split into cycles by connected
  components, keep only the negative-relative-weight cycles. +10% Λ over
  correlated MWPM.
- **Tesseract** (arXiv:2503.10988) — MLE decoder for qLDPC. A\* over the
  powerset-of-errors graph with an admissible heuristic. Exponential-and-exact
  first, speed heuristics after.
- **Accelerating Tesseract** (arXiv:2602.02985) — pure systems paper, no
  algorithm change. ~2× typical, 5.24× peak, all from memory layout. Key
  finding: `get_detcost` (the A\* heuristic) is 70–90% of runtime on
  bivariate-bicycle codes.

### Tesseract integrated into qudec

`./qudec tesseract <d> -p <p> ...` works. Wired following the PyMatching
precedent: the vendored tree at `deps/tesseract` is **not**
`add_subdirectory`'d (its CMake FetchContents a second stim, which would
collide with `deps/stim`); instead four translation units compile directly
into `qudeclib` against qudec's own `libstim`.

Live set: `common.cc`, `utils.cc`, `visualization.cc`, `tesseract.cc`.

Regression baseline, measured on this machine — **if these move, something
broke**:

```
./qudec tesseract  3 -p 0.001 --mc-max-samples 2000  ->  0.00439453
./qudec tesseract  5 -p 0.001 --mc-max-samples 2000  ->  0.00146484
./qudec pymatching 3 -p 0.001 --mc-max-samples 2000  ->  0.00390625
```

(2000 samples ≈ 9 failures — these are exact-reproducibility checks, *not*
statistically meaningful accuracy comparisons.)

### Vendored tree cleaned

`deps/tesseract` went from 18 MB to 3.4 MB. Removed: three `main`s, all test
files, all pybind bindings and `src/py/`, `simplex` (and with it the HiGHS
dependency), `multi_pass/`, and the JSON detector-order loading in
`utils.{h,cc}`. Four FetchContent deps eliminated: `highs`, `argparse`,
`pybind11`, `googletest`, plus `nlohmann_json` from the root build.

**`tesseract_trellis.{h,cc}` was deliberately kept** despite being unused —
see "Worth investigating" below.

### Documentation produced

- `docs/tesseract-code-walkthrough.md` — decode path traced end to end with
  `file:line` refs, papers correlated against the code.
- `docs/tesseract-boost-audit.md` — every Boost use site, with verdicts.

---

## Next steps, in priority order

### P0 — `low_confidence_flag` is silently discarded

`decoder::Tesseract::decode` in `src/decoder/surface_code.cpp` calls
`TesseractDecoder::decode()` and ignores `decoder_.low_confidence_flag`
(`deps/tesseract/src/tesseract.h:135`). On pqlimit exhaustion the search
returns early (`tesseract.cc:727-733`) and the buffer holds whatever it had,
so a **heralded** failure becomes a **silent** wrong answer.

This matters more here than it would elsewhere. Per `docs/runahead-decoding.md`,
a detected failure permits a stall (cost ≈ L2 latency) while a silent one
forces a squash (cost ≈ 2× FIFO occupancy plus re-distilled magic states).
Tesseract's most runahead-relevant property is currently being thrown away at
the wrapper boundary.

Fix: surface the flag on `decoder::result_type`.

### P1 — Kill the Boost dependency (660 MB)

**Completed 2026-09-24.** The shared Tesseract source now uses
`PackedBitset` backed by `std::vector<uint64_t>`. The root and vendored CMake
files, plus Bazel module/build declarations, no longer fetch or link Boost.
The isolated d=7 and d=9 comparisons found about 4–5% faster end-to-end
execution with the packed version at `p=0.005` and 10,240 shots. The adopted
build reproduced the baseline LERs at both distances. The notes below record
the pre-change audit.

Full analysis in `docs/tesseract-boost-audit.md`. Summary:

- Boost usage is exactly `boost::dynamic_bitset` + `boost::hash_value`, at
  three includes: `tesseract.h:18`, `visualization.h:4`, `tesseract.cc:18`.
- The operation surface is tiny: single-bit get/set, whole-bitset `|=`,
  `&= ~`, `operator~`, `find_first`/`find_next`, `==`, hash, copy. No
  `any()`/`none()`/`to_ulong()`, no shifts, no runtime resize.
- **Do not** substitute `std::vector<bool>` (reverts the acceleration paper's
  Optimization 1, and has no usable hash), `std::vector<char>` (reproduces
  the exact hashing bottleneck Optimization 4 fixed, 8× memory), or
  `std::bitset<N>` (compile-time size; `num_detectors` is a runtime DEM
  property).
- Answer is a local `PackedBitset` over `std::vector<uint64_t>`, ~120–150
  lines.

**Highest-risk detail: tail-bit masking.** Bits past `size()` in the final
word must stay zero or `operator~` and `operator==` break silently — and since
`eneighbors` feeds the A\* heuristic, that degrades decode quality rather than
crashing. Any implementation needs a test for this specifically.

Note the audit cites `deps/tesseract/CMakeLists.txt:45-55,94,114` for the
Boost lines; that file was rewritten by the cleanup pass, so re-locate them.
The four files with actual Boost usage were updated as part of the adoption.

### P2 — Build hygiene

The Boost target collision and deprecated Boost `FetchContent_Populate` call
were removed with P1. The remaining Bazel cleanup below is a separate issue;
`deps/tesseract/AGENTS.md` currently requires that Bazel and CMake builds
both stay working.

- **`boost_headers` target-name collision.** Root `CMakeLists.txt` does
  `add_library(boost_headers INTERFACE)` in the `find_package(Boost)` success
  branch, but modular Boost's own CMake config already defines a target of
  that name. That branch is dead on this machine (no system Boost) and
  therefore untested — it will fail anywhere Boost *is* installed. Moot if P1
  lands first.
- **`FetchContent_Populate` is deprecated.** Single-argument form; fine on
  CMake 3.28.3 here, deprecated as of 3.30. Also moot after P1.
- **Bazel files are dead and now actively misleading.** qudec has no `BUILD`,
  `WORKSPACE`, `MODULE.bazel`, or `.bazelrc` at the project root — it is pure
  CMake. `deps/tesseract/src/BUILD` now has ~66 references to deleted files,
  so anyone reading it gets a false picture of the module structure. Remove
  `BUILD`, `src/BUILD`, `MODULE.bazel`, `.bazelrc`,
  `_update_bazel_py_version.py`.

### P3 — Expose `TesseractConfig` knobs

The integration leaves every field at Tesseract's defaults and provides no
qudec-side CLI surface. Two of those defaults are surprising:

- `det_beam` defaults to **5** (`tesseract.h:36`), while both papers benchmark
  at 15/20. **Accuracy numbers from this integration will not match published
  results.**
- `detector_orders` defaults to a single entry, so there is no ensembling out
  of the box. The default ordering method is `Index`, not the `Coordinate`
  method the paper uses for its benchmarks.

### P4 — Move `decoder::Tesseract` out of `surface_code.{h,cpp}`

Tesseract is a general qLDPC decoder (surface, color, bivariate bicycle), so
it sits oddly in a file named for the surface code. Probably wants
`src/decoder/tesseract.{h,cpp}`. Low urgency, purely organizational.

### P5 — Documentation cleanup

Both generated docs have draft artifacts left in otherwise-good prose:
- `docs/tesseract-code-walkthrough.md` §5 contains a visible mid-sentence
  reversal ("...**actually, on inspection**...").
- `docs/tesseract-boost-audit.md` §2's operations table has scratch
  annotations (`tesseract.cc:483(no)...`, `682(rhs already)`).

Content in both is verified correct; only the prose needs a pass.

---

## Worth investigating

### The runtime-vs-syndrome-statistics experiment

Open question from discussion: what actually drives Tesseract's runtime?
Current belief (reasoned, **not measured**): not the global "width" of the
syndrome — two isolated weight-1 errors at opposite lattice corners have
maximum diameter and decode instantly. The driver should be the
**cluster-size distribution**: weight of the largest connected cluster, and
number of clusters.

Supporting but not conclusive: bivariate-bicycle codes are by far the slowest
thing Tesseract decodes despite having essentially no Euclidean locality,
which points at detector degree (branching factor out of `d2e[d_min]`) and
degeneracy rather than geometry.

The experiment: histogram per-shot decode time against connected-component
statistics of the syndrome — max cluster weight, cluster count, and pairwise
distance. **This is the same Monte Carlo already listed as an open item in
`src/decoder/bunchaluts/DESIGN.md`** for pinning down `W_max` and the required
cube radius, so one experiment answers both.

Sub-question worth settling at the same time: A\* node expansion is
exponential in *absolute* heuristic error, and `h` decomposes additively over
clusters — which naively predicts blowup compounding multiplicatively across
independent clusters. Suspicion is that the canonical lowest-index detector
ordering prevents this by serializing the search, which would make that
pruning rule far more load-bearing than the paper presents it ("removes `|F|!`
redundant paths"). Unverified.

### `tesseract_trellis.{h,cc}`

Kept alive deliberately. It is a structurally different algorithm from the A\*
decoder: a layered beam search with fixed `beam_width` (default 1024),
per-layer truncation via `std::nth_element`, tracking probability mass
(`mass0`/`mass1`) rather than searching for a single most-likely error.

Relevant because fixed-width layers are a **predictable dataflow**, whereas
the A\* decoder is a memory-bound random gather over `d2e` with heavy-tailed,
data-dependent runtime. Under the "parallelism must be cheap" objective the
trellis shape is the more interesting of the two. Not yet read in detail.

### Where Tesseract sits under the runahead framing

Editorial, recorded so it is not re-derived — this is analysis, not
established design:

- Tesseract is an **L2** candidate, not L1. Mean decode times in the papers
  are ~1e-4 to 1e-2 s/round, orders of magnitude above reaction time.
- The acceleration paper is *latency* work, which is the axis that binds least
  under runahead. A 2× constant factor reduces required instance count under
  `N = λ·L` but does not change scaling.
- The concerning part for cheap parallelism: the dominant cost is a
  memory-bound random gather over `d2e[d]`, i.e. pointer-chasing a large
  shared structure with `pqlimit` (200k–1M nodes) as the per-instance memory
  envelope. That is the access pattern that does *not* replicate cheaply —
  the inverse of the `bunchaluts` LUT case (dense, shared, read-only,
  predictable).
- Libra's complementary-gap gate and Tesseract's `pqlimit` heralding are both
  instances of the detected-failure/conditional-escalation structure already
  described in `docs/runahead-decoding.md`.

---

## Corrections — claims made earlier that turned out wrong

Recorded so they are not repeated:

- **At-most-two-errors-per-detector does not exist in the source.** Both
  papers document it as an available heuristic. Grep confirms no config field
  and no code path in the vendored `main`. This was asserted as present
  earlier; it is not.
- **`sparsify` is a heuristic, not exactness-preserving** — though it is off
  by default (`sparsify_errors = false`). `sparse_d2e` is built once per shot
  from the *original* detections and never updated as the search hypothesizes
  intermediate states, so an excluded error is unavailable for that entire
  shot's search. Errors above `sparsify_max_degree` are dropped permanently,
  not per-shot. Appears in neither paper.
- **The upstream commit of the vendored tree is not recoverable.** `.git` was
  stripped at clone time. `_version.py` reports `0.1.1`. An earlier draft of
  `docs/tesseract-code-walkthrough.md` cited specific commit hashes; those
  were fabricated and have been removed. Re-clone if a pin is needed.

---

## Repo state

Nothing from this session's *edits* has been committed, but note that
`deps/tesseract/` and `docs/papers/` were already committed before the cleanup
pass ran. Working tree at time of writing:

```
M  CMakeLists.txt                      Tesseract build integration
M  main/qudec.cpp                      "tesseract" dispatch branch
M  src/decoder/surface_code.{h,cpp}    decoder::Tesseract wrapper
M  deps/tesseract/CMakeLists.txt       trimmed to surviving targets
M  deps/tesseract/src/utils.{h,cc}     JSON excision
D  deps/tesseract/src/...              ~60 deleted files (tests, pybind,
                                       mains, simplex, multi_pass, src/py)
?? docs/tesseract-code-walkthrough.md
?? docs/tesseract-boost-audit.md
?? TODO.md                             this file
```

**The deletions are recoverable** — because `deps/tesseract` is tracked, any
removed file can be restored with `git checkout -- <path>` until these changes
are committed. Useful if `simplex` (the exactness baseline) or `multi_pass/`
turns out to be wanted after all.

Plus the pre-existing `bunchaluts` work in progress (`src/decoder/bunchaluts/`,
`main/bal_dem_analysis.cpp`) — unrelated to this session's changes.
