# Paper notes: Libra and Tesseract

PDFs (and extracted text) in this directory. Read 2026-09-23. Both are Google
Quantum AI accuracy-oriented decoders; both are naturally **L2** candidates in
the runahead framing (`docs/runahead-decoding.md`).

---

## Libra — "Improved accuracy for decoding surface codes with matching synthesis"

Cody Jones, Google Quantum AI, arXiv:2408.12135 (Aug 2024). `libra-2408.12135.pdf`

### What it does

An **ensemble** decoder for surface codes that does not pick a winner from the
ensemble — it *splices* the ensemble's solutions together. Claimed result:
Λ improves ~10% over correlated MWPM (Λ₇,₁₁ 3.64 → 4.02), saturating at
ensemble size ≈ 60.

### Problem framing

Decoding = **minimum-weight hypergraph perfect matching** (MWHPM), NP-hard.
Hypergraph `G=(V,E)`: vertices are detectors, hyperedges are error mechanisms,
`w(e) = log((1-p)/p)`. Syndrome `S`; find `min w(e)` s.t. `s(e) = S` where `s`
is mod-2 parity of incident vertices. Each hyperedge also carries an observable
list `o(e)`; `o` defines the equivalence class of a solution.

### The core algorithm: matching synthesis

Given two matchings `e`, `f` for the same syndrome:

1. `d = e ⊕ f` — symmetric difference. `s(d) = ∅`.
2. Split `d` into **cycles** `{dᵢ}` — pieces with `s(dᵢ) = ∅` *and*
   `o(dᵢ) = ∅` (the observable condition matters: it keeps the equivalence
   class fixed). Splitting is done by **connected components** of `d`, which is
   cheap and works well on 3D-local surface-code hypergraphs. (Not guaranteed
   minimal — two cycles touching at a vertex merge into one component.)
3. Compute each cycle's **relative weight**
   `w(c|e) = w(c\e) − w(c∩e)`. Note this costs O(|c|), *not* O(|matching|) —
   this is what makes the whole thing affordable.
4. Keep only the improving cycles (`w(c|e) < 0`) and form
   `g = e ⊕ (⊕_{c∈N} c)`.

When `N` is a strict nonempty subset of the cycles, `g` is a **new** matching
better than both inputs. That's the whole trick.

Pieces with `s(dᵢ)=∅` but `o(dᵢ)≠∅` are logical operators, not cycles. They're
used two ways: to hop equivalence class, and pairwise (`dᵢ ⊕ dⱼ` for
`o(dᵢ)=o(dⱼ)`) to manufacture cycles connected components would have missed.

### The Libra decoder proper

- **Init:** complementary matching (correlated MWPM, unperturbed graph) gives
  one representative per equivalence class — two, for a memory experiment.
- **Gate:** the complementary gap is computed here. If gap ≥ 20 dB, stop and
  return correlated MWPM. The ensemble only runs on "hard" (small-gap) shots,
  whose frequency scales with the logical error rate and so falls off
  exponentially in `d` below threshold. **Average runtime is therefore basically
  correlated MWPM's**, even with 100 ensemble members.
- **Ensemble:** correlated MWPM with per-mechanism probabilities perturbed
  `p'ᵢ = pᵢ·rᵢ`, `rᵢ ~ Lognorm(0,σ²)`. Half the members at `σ=ln 2`, half at
  `σ=ln 4` (inhomogeneous helped slightly). Perturbed weights are used *only*
  to generate diverse matchings; all weight evaluation uses true weights.
- **Synthesize:** fold ensemble members in one at a time (Fig. 3), or as a
  **log-depth binary tree** (Fig. 4) — the paper flags the tree as the
  parallel-friendly option.
- **Degeneracy ("Libra-degen"):** small *positive*-relative-weight cycles are
  kept in a fixed-size max-heap (30). Each generates a near-best alternative
  solution, so `n` non-overlapping cycles encode `2ⁿ` solutions. Libra sums
  relative probabilities per connected component of cycles to approximate the
  equivalence class probability, truncating to low-order combinations when a
  component has ≥10 cycles. The heap is cleared whenever an improving cycle is
  found, so synthesis is run **twice** — pass 2 mostly exists to collect
  positive cycles. Worth ~1% extra Λ.

### Why it works (author's conjecture, §VI)

Ensemble members differ from optimal by a random distribution of cycles with
**positive mean** relative weight. Picking the globally-lowest-weight member
("global ensembling") requires one member to get lucky across the *whole*
syndrome volume — and by CLT that tail shrinks as the volume `~d²r` grows.
Hence global ensembling **degrades with `d` and `r`** in their data, while Libra
is flat: it picks the best cycle in each local neighborhood independently.
Windowing/partitioning would mitigate the volume effect for global ensembling,
but the minimum window is `~d³`.

### Runahead-relevant observations

- The complementary-gap gate is exactly the **detected-failure → conditional
  escalation** structure in `docs/runahead-decoding.md` §"Implications for L1
  design", implemented at the decoder level rather than the architecture level.
- **Parallelism cost:** the ensemble is embarrassingly parallel and synthesis is
  log-depth — good. But each instance is a *full correlated MWPM*, so the
  marginal cost of the N-th instance is high. The paper's own discussion (§VIII)
  proposes exactly the fix our framing wants: swap correlated MWPM for
  **uncorrelated matching / Union-Find / clustering / RG decoders**, accept a
  larger ensemble, and win on depth via the log-tree on FPGA/GPU. That is a
  cheap-parallelism argument in all but name.
- Cycle evaluation being O(cycle size) and cycles being local is the property
  that makes the synthesis step itself cheap to parallelize.

---

## Tesseract — "A Search-Based Decoder for Quantum Error Correction"

Beni, Higgott, Shutty; Google Quantum AI, arXiv:2503.10988v2 (Aug 2025).
`tesseract-2503.10988.pdf`. Open source: github.com/quantumlib/tesseract-decoder

### What it does

A **Most-Likely Error** decoder for qLDPC codes generally (surface, color,
bivariate bicycle). Inverts the usual strategy: rather than starting polynomial
and adding accuracy heuristics, it **starts exponential-and-exact and adds
speed heuristics**. ~5× faster than an integer-program decoder at matching
accuracy; 1–2 orders of magnitude better LER than BP+OSD on bicycle codes.

### Problem framing

Same DEM vocabulary as Stim. `w(e) = -log(p/(1-p))`. MLE problem:
`argmin_{F⊆E, D(F)=S} w(F)`.

Reformulated as **shortest path on the powerset graph** `G = (2^E, T, w)`:
nodes are error subsets, an edge `F → F∪{e}` has weight `w(e)`, START = `∅`,
EXIT = `{F : D(F) = S}`. Trivially: path cost = total error weight, so
MLE ≡ shortest START→EXIT path (Thm 3.1). The graph is exponentially large;
everything interesting is in how the search is pruned.

### The algorithm

**A\*** over that graph with an **admissible heuristic** — so exactness is
preserved. The heuristic (Alg. 2/3):

```
h(F) = Σ_{d ∈ x}  min_{e ∈ E(d), e ∉ J}  w(e) / |x ∩ D(e)|      x = S ⊕ D(F)
```

i.e. for each residual detection event, the cheapest per-detector-amortized
mechanism that could explain it. For LDPC codes this updates in O(1) per added
error.

**Pruning predicate `P(F,F')`** (Alg. 1) — restricts which errors may be added:

- only errors incident to the **lowest-index** residual detector `d_min`;
- minus `GetForbiddenErrorsByPrecedence(F)` — canonicalizes the insertion
  order so there is exactly one path to each `F` instead of `|F|!`. **Exactness
  preserved**, and it turns `G` into a *tree*, so no visited-set is needed.
- minus `GetForbiddenErrorsAtMostTwo(F)` — forbids >2 errors on any one
  detector. **Breaks exactness**, but prunes hard.

**Speed heuristics on top:**

| heuristic | what it does | exact? |
|---|---|---|
| **beam** | don't visit `F` with `r(F) > r_min + beam`, where `r(F)=\|S⊕D(F)\|`. ~20 works. | no |
| **pqlimit** | cap priority-queue insertions; on exceeding, declare **low-confidence** | n/a (herald) |
| **ensemble reordering** | Alg. 1 needs a total detector order; sample `z~N(0,1)^t`, sort detectors by `⟨coord, z⟩`; take min-cost solution over orderings | yes |
| **beam climbing** | try every beam in `{0..B}`, different random ordering each | no |
| **no-revisit detections** | skip nodes with an already-seen residual `D(F)⊕S` | no |
| **detection penalty** | add `c·\|D(F)⊕S\|` to the priority key | no |

Benchmark settings: "short beam" = beam 15 + climbing + 16 orderings, pqlimit
200k; "long beam" = beam 20 + climbing + 21 orderings, pqlimit 1M.

### Results worth remembering

- Accuracy ≈ integer program at `p ≤ 0.001`, ~5× faster. Error floor appears at
  `p=0.002` for large `d`.
- 1–2 orders of magnitude better than BP+OSD on bivariate bicycle codes.
- Appendix A: **uncorrelated** BP+OSD beats correlated BP+OSD, because Y errors
  create 4-cycles `(S_X, Y_i, S_Z, Y_j)` in the full Tanner graph plus extra
  degeneracy → trapping sets. Tesseract's advantage over BP+OSD is largely that
  it exploits Y errors properly.
- Code comparison: gross code `[[288,12,10]]` beats `d=13`–`15` surface codes
  under Tesseract (14–19× qubit saving) vs only 10× under BP+OSD/correlated
  matching. With 10× noisier long-range couplers (NLR10) that collapses to 4×
  (Tesseract) / 2× (BP+OSD).
- Concurrent independent work: Ott/Hetényi/Beverland "Decision-Tree Decoder"
  (arXiv:2502.16408) — same A* idea; their "syndrome height" *is* an admissible
  heuristic. Tesseract adds canonical path ordering and beam cutoffs; DTD
  explores fancier heuristics and BP-guided search.
- Follow-up exists: arXiv:2602.02985, "Accelerating the Tesseract Decoder"
  (not read).

### Runahead-relevant observations

- **`pqlimit` is a heralded failure.** The paper counts low-confidence outcomes
  as logical errors, but in a runahead system a heralded failure is exactly the
  cheap case — stall instead of squash. Same lever as Libra's complementary gap.
- **Parallelism cost is the problem.** A* with a priority queue is a large,
  mutable, per-instance state with a data-dependent, heavy-tailed runtime.
  `pqlimit` of 200k–1M nodes is the per-instance memory envelope. Under
  `N = λ·L`, replicating that `N` times is exactly the expensive kind of
  parallelism — the opposite of the shared-read-only-table case. Its ensemble
  axes (reordering, beam climbing) parallelize trivially, but they multiply
  instance count rather than reduce instance cost.
- Mean decode time per round is ~1e-4 to 1e-2 s in Fig. 2 — many orders of
  magnitude above reaction-time requirements, which is the regime runahead
  exists to make usable. The open question for both papers is whether their
  *throughput* can be bought at acceptable hardware cost, not their latency.

---

## Side-by-side

| | Libra | Tesseract |
|---|---|---|
| target codes | surface (topological, matchable) | general qLDPC |
| base primitive | correlated MWPM ensemble | A\* on the powerset graph |
| optimality | heuristic, improves on MWPM | exact if only precedence-pruning + A\*; heuristic with beam |
| escalation gate | complementary gap < 20 dB | — (beam/pqlimit instead) |
| heralded failure | no (falls back to MWPM) | yes, via `pqlimit` |
| per-instance state | a full MWPM solve | priority queue up to `pqlimit` nodes |
| parallel structure | ensemble + log-depth synthesis tree | ensemble over detector orderings/beams |
| accuracy claim | +10% Λ vs correlated MWPM | ≈ integer program, 5× faster; 100× vs BP+OSD on BB codes |

---

## Accelerating the Tesseract Decoder — arXiv:2602.02985

Grbic (Rice, intern), Aghababaie Beni, Shutty (Google Quantum AI), v2 Feb 2026.
`accel-tesseract-2602.02985.pdf`

### What it is

**Not an algorithms paper.** Zero changes to the search, the heuristic, or the
accuracy — every optimization is verified against the integer-program decoder to
be output-identical. It is a pure systems/profiling case study on the existing
open-source C++ implementation. Result: ~2× across most code families, >2.5×
often, **5.24× peak** on the heaviest bivariate-bicycle config (36,853 s →
7,048 s for 1000 shots).

Method: profile with HPCToolkit, then `perf`; four optimizations applied
incrementally and measured one at a time; validated on three microarchitectures
(Xeon W-2135, Cascade Lake, Sapphire Rapids) and two beam configs.

### The four optimizations

| # | change | why |
|---|---|---|
| 1 | `std::vector<bool>` → `std::vector<char>` | bit-packing forces proxy objects doing shift+mask on every access; blocked-error flags and syndrome patterns are touched constantly. Trading memory for byte-addressability wins. +13–42% |
| 2 | **SoA → AoS** in `get_detcost`: fold the blocked-errors vector and the fired-detector-count vector into one `DetectorCostTuple {uint32_t, uint32_t}` array | the two were always co-accessed at the same arbitrary index → two scattered cache lines per probe. One 8-byte struct → one line, plus prefetcher-friendly. **Biggest single win.** |
| 3 | **early exit** in `get_detcost` | precompute each error's cost *lower bound* = `likelihood_cost / max possible fired detectors`, sort each detector's error list by it; break as soon as `min_cost ≤` next bound |
| 4 | `boost::dynamic_bitset` for syndrome patterns + `boost::hash_value` | the no-revisit heuristic hashes the residual pattern constantly; a naive `seed*31+el` loop over `vector<char>` was hot. dynamic_bitset gives hardware bitwise ops *and* claws back Opt-1's memory growth. Biggest effect on surface codes / transversal CNOT, where short beams make "bushy" search graphs with many redundant paths. |

### The key finding

**`get_detcost` — the A\* heuristic — is the whole story.** It is 70–90% of
decode time for bivariate-bicycle codes, >60% for color codes, ~40% for surface
codes and transversal CNOT. Opt 2 alone gave 2.75× on the worst BB config.
Measured effect: LLC miss rate inside the kernel drops by up to **~90%**
(e.g. color code SI1000 p=1e-3 d=7), L1 by 16–30%.

Speedups are smaller under the long-beam config (peak 4.31× Cascade Lake,
2.8× Sapphire Rapids) than short-beam, but hold across all three machines —
the authors read this as evidence the wins are memory-hierarchy-fundamental,
not hardware quirks.

### Are the optimizations online? **Yes — all merged into `main`.**

Verified 2026-09-23 against `github.com/quantumlib/tesseract-decoder`:

| paper opt | PR | merged |
|---|---|---|
| 1, `vector<bool>` → `vector<char>` | [#25](https://github.com/quantumlib/tesseract-decoder/pull/25) | 2025-05-28 |
| 2 + 3, AoS + early exit | [#34](https://github.com/quantumlib/tesseract-decoder/pull/34) | 2025-06-19 |
| 4, `boost::dynamic_bitset` hashing | [#57](https://github.com/quantumlib/tesseract-decoder/pull/57) | 2025-07-30 |

Present in `main` today: `DetectorCostTuple` and `ErrorCost` structs and the
new `get_detcost` signature in `src/tesseract.h`; `boost::dynamic_bitset`
include, a `std::hash` specialization delegating to `boost::hash_value`, and
`edets_bitsets` in `src/tesseract.cc`.

**`main` has in fact moved past the paper.** The published Listing 2 still
divides (`ec.likelihood_cost / dct.detectors_count`) inside the hot loop.
Current `main` has eliminated division entirely via cross-multiplication,
carrying `min_det_cost_det_count` alongside `min_cost` and comparing
`a.cost * b.count  vs  b.cost * a.count`; the single divide is deferred to the
return. So the real upstream code is a further step ahead of what's written up.

### Runahead-relevant observations

- This paper is **latency work, and under runahead that is the axis that
  matters least for L2.** 2× off a 1e-2 s/round decode is still ~5 orders of
  magnitude from reaction time. It does not move Tesseract toward being an L1.
- It *does* matter for **throughput per unit hardware**, which is the axis that
  does bind. Under `N = λ·L`, a 2× constant-factor win is a 2× reduction in the
  required instance count — real, but it does not change the scaling.
- The deeper reading is a warning about Tesseract as an L2 under a
  cheap-parallelism objective: the dominant cost is a **memory-bound random
  gather** over `d2e[d]`, i.e. pointer-chasing a large shared structure. That
  is precisely the access pattern that does *not* replicate cheaply, and
  precisely the opposite of the `bunchaluts` LUT case (dense, shared,
  read-only, predictable). 70–90% of runtime in a cache-missing kernel says the
  bottleneck is the memory hierarchy, not arithmetic — so N parallel instances
  contend for exactly the resource that is already saturated.
- Useful methodological precedent regardless: profile-first, verify
  bit-identical output against an exact baseline, report per-optimization
  deltas. Reasonable template if we ever benchmark `bunchaluts` seriously.
