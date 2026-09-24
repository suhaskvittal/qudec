# Tesseract Boost audit: can the 660 MB dependency go away?

**Update (2026-09-24):** The local `PackedBitset` recommendation was adopted
in `deps/tesseract/src/packed_bitset.h`; the build declarations no longer
fetch Boost. The analysis below documents the original Boost-based snapshot
and retains its historical source references.

Written for an agent picking this up cold. Read-only analysis, no code
changed. Scope: every Boost usage in `deps/tesseract/src/` as of 2026-09-23
(same vendored snapshot documented in `docs/tesseract-code-walkthrough.md`).
Cross-referenced against `docs/papers/accel-tesseract-2602.02985.txt`
(§4.1, §4.4), which documents the exact history of the data structure now
under review. All `file:line` references were read directly; none inferred
unless marked "inference."

**Note on concurrent work**: another agent is deleting dead code from
`deps/tesseract` (executables, tests, pybind, JSON, possibly `simplex` /
`multi_pass/`) concurrently with this audit. None of the files this audit
depends on (`tesseract.{h,cc}`, `visualization.{h,cc}`) are candidates for
that deletion pass, and grepping the rest of `src/` turned up zero other
Boost references, so this audit's file set is stable regardless of what the
other agent removes.

## Summary table

| Facility | Sites | Verdict |
|---|---|---|
| `boost::dynamic_bitset<>` — syndrome/detector bit-vector | `tesseract.h:164`; `tesseract.cc:350-374, 474, 523, 551, 563` (+ every subsequent read/write of those variables to `tesseract.cc:724`); `visualization.h:16`, `visualization.cc:38` | **[LOCAL-IMPL]** — write a ~100-150 line local packed-bitset class exposing exactly the operations used below. Neither raw STL option is safe (see rationale). |
| `boost::hash_value` / `boost::functional/hash.hpp` (`std::hash<boost::dynamic_bitset<>>` specialization) | `tesseract.cc:18, 86-93` | **[LOCAL-IMPL]** — replaced automatically once the bitset itself is local; give the local type its own `hash()` (word-combine, no `std::hash` specialization games needed if it's used as a class member or a free function object). |
| `Boost::boost` / vendored Boost headers (CMake) | `CMakeLists.txt:60-82, 126`; `deps/tesseract/CMakeLists.txt:45-55, 94, 114` | **[REPLACE]** (conditional on the above two landing) — the entire 660 MB `FetchContent` fetch of Boost 1.83.0 exists *solely* to provide `<boost/dynamic_bitset.hpp>` and `<boost/functional/hash.hpp>`. No other Boost header is referenced anywhere in `deps/tesseract/src/` or the rest of the qudec tree (verified by grep, see §0). |

Bottom line up front, expanded at the end: **yes, removable**, but not by
swapping in `std::vector<bool>`, `std::vector<char>`, or `std::bitset<N>` —
each fails a requirement the paper already established by measurement. The
correct move is a small hand-written packed-bitset class matching the
handful of operations actually exercised (enumerated in §2). Estimated
effort: one new header (~120-150 lines), edits to 4 files, no algorithmic
changes.

## 0. Confirming Boost's footprint is exactly this one facility

```
$ grep -rn "boost" deps/tesseract/src/
BUILD:185,209                @boost//:dynamic_bitset          (Bazel target, unused by the CMake build)
tesseract.h:18                #include <boost/dynamic_bitset.hpp>
tesseract.h:164               boost::dynamic_bitset<>& detectors  (parameter)
visualization.h:4             #include <boost/dynamic_bitset.hpp>
visualization.h:16            const boost::dynamic_bitset<>&  (parameter)
visualization.cc:38           const boost::dynamic_bitset<>&  (parameter)
tesseract.cc:18               #include <boost/functional/hash.hpp>
tesseract.cc:87-91             std::hash<boost::dynamic_bitset<>> specialization, boost::hash_value(bs)
tesseract.cc:350,351,359,369,474,521,523,551,563   construction/use sites (detailed in §2)
```

No `boost::multi_index`, no `boost::graph`, no Boost.Asio, no
`boost::optional`/`variant` (all superseded by `<optional>`/`<variant>` in
this codebase already), nothing else. The CMake comment at
`CMakeLists.txt:60` and `deps/tesseract/CMakeLists.txt:45` both say as much
explicitly ("Boost (`boost::dynamic_bitset`, `boost::hash_range`)"). This is
a single-facility dependency: **removing `boost::dynamic_bitset` removes all
of Boost.**

## 1. What `boost::dynamic_bitset<>` is doing here, in context

Two disjoint roles, both inside `TesseractDecoder`'s A* search
(`tesseract.cc`, decode path documented in `docs/tesseract-code-walkthrough.md`
§1, §6):

### 1a. One-time setup: neighbor-detector computation (`initialize_structures`, `tesseract.cc:321-412`)

- `edets_bitsets` (`tesseract.cc:350-351, 354`): one bitset per error
  (`num_errors` of them), each `num_detectors` bits wide, with a bit set for
  every detector that error's symptom activates. Built once per decoder
  construction (amortized across all shots decoded by that instance).
- `neighbor_set` (`tesseract.cc:359-372`): per error, OR together the
  `edets_bitsets` of every error sharing a detector with it
  (`neighbor_set |= edets_bitsets[oei]`), then subtract the error's own
  detector set (`neighbor_set &= ~edets_bitsets[ei]`), then enumerate the
  remaining set bits via `find_first()`/`find_next()` to populate
  `eneighbors[ei]` (a `vector<int>`, used later purely as a plain index
  list — the bitset itself is discarded after this loop).

This is **not** hot-path (runs once per decoder, not per shot, not per
search node), so its performance is not what the paper's Optimization 4
targeted. But it does exercise the full bitwise-algebra surface (`|=`, `&=`,
unary `~`, `find_first`, `find_next`, `npos` sentinel) that any replacement
must support.

### 1b. Hot path: per-search-node syndrome tracking (`decode_to_errors_with_graph`, `tesseract.cc:504-728`)

This is exactly what §4.4 of the accel paper describes and optimized.

- `initial_detectors` (`tesseract.cc:523`): one bitset, `num_detectors` bits,
  set from the shot's `detections` list (`initial_detectors[d] = true`,
  `tesseract.cc:532`).
- `detectors` (`tesseract.cc:563`): copy-constructed from `initial_detectors`
  **once per priority-queue pop** (i.e., once per A* node expanded), then
  mutated bit-by-bit in `flip_detectors_and_block_errors`
  (`tesseract.cc:489`, `detectors[d] = !detectors[d]`) while walking the
  parent chain of that node.
- `next_detectors` (`tesseract.cc:551`, default-constructed empty, then
  assigned `next_detectors = detectors` at `tesseract.cc:674` — full-bitset
  copy-assign, once per candidate error `ei` considered for each popped
  node), then flipped bit-by-bit per detector in that error's symptom
  (`tesseract.cc:681`, `next_detectors[d] = !next_detectors[d]`).
- `visited_detectors` (`tesseract.cc:521`):
  `unordered_map<size_t, unordered_set<boost::dynamic_bitset<>>>`, keyed by
  number of activated detectors. `boost::dynamic_bitset<>` is the **key type
  of an `unordered_set`**, which requires both `operator==` (provided by
  Boost) and a specialized `std::hash` (hand-written at `tesseract.cc:86-93`,
  delegating to `boost::hash_value`). Exercised via:
  - `.insert(detectors)` (`tesseract.cc:602`) — the "No-Revisit Detection
    Events" pruning check: if this exact detector pattern was already
    visited at this depth, skip the node.
  - `.find(next_detectors) != .end()` (`tesseract.cc:691-692`) — same check
    for a not-yet-committed candidate state.
  - `.clear()` on individual map buckets (`tesseract.cc:634`) as the beam
    window slides — no bitset-level op, just map/set housekeeping.
- Per-bit read in a plain loop for debug/visualization output
  (`tesseract.cc:582-586, 621-627`, `visualization.cc:42-46`) — only under
  `config.verbose` / `config.create_visualization`, not hot path.

**Frequency**: `detectors` is copied once per node popped from the priority
queue (bounded by `pqlimit`, default 200,000, `tesseract.h:37`); `next_detectors`
is copy-assigned and hashed once per candidate error considered per popped
node (fan-out = number of errors touching the node's minimum detector,
`active_d2e[min_detector]`, `tesseract.cc:661`). For real workloads this is
the dominant per-shot cost the accel paper measured (§4.4: "the function
responsible for hashing these patterns consumed a significant portion of
decoding time").

## 2. Exhaustive list of operations exercised

From walking every use site above, the *complete* operation surface required
of any replacement type is:

| Operation | Site(s) | Notes |
|---|---|---|
| Construct with `(n_bits)` or `(n_bits, false)` | `tesseract.cc:351, 359, 523` | all-zero init |
| Default-construct (empty) then later copy-assign | `tesseract.cc:551, 674` | `next_detectors` starts empty, sized on first assignment |
| Copy-construct | `tesseract.cc:563` (`detectors = initial_detectors`, at decl) | full-word-range copy |
| Copy-assign | `tesseract.cc:674` (`next_detectors = detectors`) | same size both sides in this codebase — no true resize ever needed at runtime beyond the two setup calls |
| `operator[]` read | `tesseract.cc:483(no)... 582, 596(no), 621-627, 640-641, 652, 666, 682(rhs already), 696, 707`, `visualization.cc:43` | boolean test of one bit |
| `operator[]` write (`= true`, `= !x`) | `tesseract.cc:354, 489(via ref), 532, 681` | single-bit set/flip |
| `operator|=` (bitset OR) | `tesseract.cc:363` | whole-bitset, word-parallel in Boost |
| `operator&=` with `~other` (AND-NOT) | `tesseract.cc:367` | whole-bitset |
| `operator~` (unary complement) | `tesseract.cc:367` | produces temporary, whole-bitset |
| `find_first()` / `find_next(pos)` / `npos` sentinel | `tesseract.cc:369-370` | word-skipping bit iteration |
| `operator==` (via `unordered_set`/`unordered_map` key use) | implicit in `.insert()`/`.find()` at `tesseract.cc:602, 691` | must be exact, whole-word compare |
| Hash (`std::hash` specialization → `boost::hash_value`) | `tesseract.cc:86-93` | used as `unordered_set<bitset>` key hash |
| No use of: `count()`, `any()`, `none()`, `all()`, `to_ulong()`, `to_string()`, `resize()` (beyond initial sizing), `set(pos,len)` range form, `<<`/`>>` shift, `test_set()` | — | confirmed absent by grep across `tesseract.cc`/`.h` and `visualization.*` |

This is a narrow surface: single-bit get/set, whole-bitset OR/AND-NOT/NOT,
first/next-set-bit iteration, equality, hash, copy. Nothing resize-heavy,
nothing exotic.

## 3. Why the two "obvious" STL swaps are traps

This is the part the task brief specifically warns about, so it's argued
directly against the paper's evidence rather than asserted.

### `std::vector<bool>` — **do not use**

This is literally what Tesseract started with, and Optimization 1 in the
accel paper (`docs/papers/accel-tesseract-2602.02985.txt:219-266`) replaced
it precisely because its bit-packing proxy references made per-element
read/write slow (`operator[]` returns a proxy object requiring shift+mask on
every access, and the paper shows this dominated runtime — not memory
traffic). §1b above shows `detectors[d]`/`next_detectors[d]` single-bit
read/write happening on the hot path (once per detector touched, per
candidate error, per popped node) — exactly the access pattern Optimization
1 was written to eliminate. Reintroducing `vector<bool>` here would revert a
measured, documented win. There is also no portable, efficient way to hash a
`vector<bool>` (no `std::hash` specialization exists; a correct one would
have to either iterate bit-by-bit through the proxy, i.e. exactly the
overhead being avoided, or reach into libstdc++/libc++-specific internals to
get the underlying word storage, which is non-portable and not something to
build on for a project already citing Boost portability as a concern).
**Verdict: rejected, not merely deprioritized.**

### `std::vector<char>` — **do not use for the hot-path bitsets**

This was Tesseract's *intermediate* state after Optimization 1 and before
Optimization 4 (`accel-tesseract-2602.02985.txt:371-442`). It fixed the
proxy-access problem (byte-addressable, direct access) but:
- **No word-at-a-time bitwise ops.** `|=`/`&=`/`~` would each need an
  explicit per-byte loop (or the compiler would have to auto-vectorize a
  byte-wise loop, which is a much weaker guarantee than an intentionally
  packed 64-bit-word implementation).
- **Hashing was the specific, measured bottleneck the paper eliminated.**
  §4.4 shows the "unoptimized" hash (Listing 3 in the paper, a per-byte
  accumulation loop over `vector<char>`) was slow enough that the authors
  went looking for a "highly optimized, built-in hashing function" and
  switched to `boost::dynamic_bitset` + `boost::hash_value` specifically to
  fix it, citing this as the highest-impact optimization for Surface Codes
  and Transversal CNOT protocols with long beams (`accel-tesseract-...:440-442`).
  A `vector<char>` reintroduction reproduces that exact bottleneck.
- **8x worse memory footprint** than a packed representation: 1 byte/detector
  vs. 1 bit/detector. `visited_detectors` can hold on the order of
  `pqlimit` (200,000 default) live bitset instances simultaneously during a
  single decode; at, say, 1000 detectors that's ~200 MB in `vector<char>`
  vs. ~25 MB packed (128 bytes/instance + `vector` overhead). The paper
  makes the same observation in reverse when justifying Optimization 4
  (`accel-tesseract-...:437-439`: "mitigated the increased memory footprint
  from our initial optimization with `std::vector<char>`").

**Verdict: rejected for the hot-path bitsets** (`detectors`, `next_detectors`,
`initial_detectors`, and anything used as an `unordered_set` key). It remains
a fine choice elsewhere in the codebase for things that are *not* bitsets
(and indeed is used generically elsewhere per the paper) — this verdict is
scoped to the `boost::dynamic_bitset` replacement only.

### `std::bitset<N>` — **structurally disqualified**

`std::bitset` requires a compile-time constant size template parameter. The
task's own framing note is correct and worth stating precisely: `num_detectors`
is read from the DEM at decoder-construction time
(`initialize_structures(size_t num_detectors)`, `tesseract.h:159`,
called from the constructor per `docs/tesseract-code-walkthrough.md` §1 step
8) — it is a **runtime** value that varies by code family and distance
(quadratic-ish growth with distance for surface codes; QLDPC/BB codes have
their own detector counts). There is no single `N` to pick. Workarounds
(instantiate `std::bitset<N>` for a handful of `N` values and dispatch at
runtime via a variant/switch, or pick one large `N` and waedge every DEM into
it with unused high bits) are both worse than writing a dynamically-sized
type directly: the first is combinatorial-blowup ugly and still needs a
fallback, the second wastes memory proportional to
`max_supported_detectors - num_detectors` for every one of the potentially
hundreds of thousands of live `visited_detectors` entries. The accel paper
reaches the identical conclusion for the identical reason
(`accel-tesseract-...:422-427`: "`std::bitset` ... requires a static size
determined at compile-time, making it unsuitable for Tesseract's dynamically
sized syndrome patterns"). **Not a candidate at all**, regardless of
performance — it is disqualified on interface grounds alone.

## 4. The actual recommendation: `[LOCAL-IMPL]`

Given the narrow operation surface in §2, and that a correct, fast answer
needs (a) runtime-determined size, (b) word-parallel bitwise ops, (c) a
cheap, well-distributed hash usable as an `unordered_set` key, and (d) no
proxy-object element access on the hot path — write a small local type. This
is exactly what `boost::dynamic_bitset` *is*, minus the parts of its API
Tesseract never touches (stream I/O, `to_string`, `to_ulong`, range-`set`,
shifts, `count`, comparison operators other than `==`, iterators, etc.).

Sketch (not code to be written by this audit, since this task is read-only,
but concrete enough to size the work):

```cpp
class PackedBitset {
  std::vector<uint64_t> words_;
  size_t num_bits_;
 public:
  explicit PackedBitset(size_t num_bits, bool value = false);
  PackedBitset() : num_bits_(0) {}          // matches next_detectors's default-ctor use
  bool operator[](size_t i) const;           // word/bit test, no proxy needed for read
  void set(size_t i, bool v);                // explicit setter (replaces write-through operator[])
  void flip(size_t i);                       // used at tesseract.cc:489, 681
  PackedBitset& operator|=(const PackedBitset&);   // word-parallel OR, tesseract.cc:363
  PackedBitset& and_not(const PackedBitset&);      // fuses &= and ~ used together at tesseract.cc:367
  size_t find_first() const;                 // ctz-based, tesseract.cc:369
  size_t find_next(size_t pos) const;        // ctz-based, tesseract.cc:370
  static constexpr size_t npos = SIZE_MAX;
  bool operator==(const PackedBitset&) const = default;  // vector<uint64_t>'s == is word-parallel already
  size_t hash() const;                       // word-combine (e.g. FNV-1a or boost::hash_combine's mix, folded over words_)
};
```

- `operator[]` read-only avoids ever needing a proxy reference type — every
  write site in the actual code (`tesseract.cc:354, 489, 532, 681`) is either
  a known `= true` or a `flip`, both expressible as named methods instead of
  `operator[]=`. This sidesteps the entire class of proxy-object overhead
  that motivated Optimization 1 in the first place, for both read *and*
  write.
- `find_first`/`find_next` are `std::countr_zero` (C++20, `<bit>`) over
  non-zero words — same algorithmic shape as Boost's implementation, no
  library dependency.
- `operator==` on `std::vector<uint64_t>` is already word-parallel
  (effectively `memcmp` for trivial element types under most
  implementations) — this is a **strict improvement** over anything
  byte-wise.
- Hash: fold `std::hash<uint64_t>` (or raw multiply-xor) over `words_`,
  something in the shape of `boost::hash_range` (which is literally what
  `tesseract.cc:18`'s comment says is being borrowed for `boost::hash_value`
  on a `dynamic_bitset` — it hashes over the internal blocks). This is a
  ~5-line loop, not a research problem.
- **No `std::hash` specialization gymnastics needed.** The current code
  specializes `std::hash<boost::dynamic_bitset<>>` in `namespace std`
  (`tesseract.cc:85-94`) because `boost::dynamic_bitset` doesn't provide a
  `std::hash` itself, only `boost::hash_value` (Boost's own hashing
  customization point, ADL-found `hash_value` free function). A local class
  can instead provide its own `hash()` method and either (a) specialize
  `std::hash<PackedBitset>` to call it (same pattern as today, one function,
  trivial), or (b) pass an explicit hasher functor as the `unordered_set`'s
  second template argument. Either is a direct, drop-in replacement for
  `tesseract.cc:86-93` and `tesseract.cc:521`'s `unordered_set` declaration.
- **Memory footprint**: `sizeof(PackedBitset)` ≈ `sizeof(vector<uint64_t>)`
  (24 bytes: ptr/size/capacity) + `sizeof(size_t)` (8 bytes for `num_bits_`)
  + `ceil(num_bits/64)*8` bytes of heap storage — this is essentially
  identical to `boost::dynamic_bitset<>`'s own layout (a `vector<block_type>`
  plus a bit-count field), so there is no footprint regression relative to
  the status quo, and it is the ~8x improvement over `vector<char>` noted
  in §3.
- **Estimated size**: ~120-150 lines including the constructor, the six
  methods above, `operator=`/copy semantics (defaulted, since `words_` is a
  plain `vector`), and light documentation. This is a genuinely small,
  self-contained header — call it `packed_bitset.h` alongside `common.h`/
  `utils.h`.

**Caveats / regression risk** (be honest, per the task's own instruction):
- The AND-NOT and OR at `tesseract.cc:363, 367` run over `num_errors` pairs
  during one-time setup; not perf-critical, but must be bit-exact — any
  off-by-one in word/tail-bit masking (e.g. the last, partially-used word
  when `num_bits_` isn't a multiple of 64) will silently corrupt
  `eneighbors`, which feeds the A* heuristic. **Tail-bit masking after every
  mutating op is the single highest-risk detail** — if not zeroed correctly,
  `operator==`/`hash()` will disagree with logical bitset equality
  (garbage high bits differing between two logically-equal patterns), which
  would silently break the "No-Revisit Detection Events" pruning
  (`tesseract.cc:602`) — false cache misses, not wrong answers, but a
  performance regression that would be very hard to notice without a targeted
  test.
- This is a decode-correctness-adjacent data structure sitting under an A*
  search; any change here should be validated against the existing
  `tesseract.test.cc` (if it survives the concurrent dead-code cleanup) or an
  equivalent regression test comparing decoded outputs before/after on a
  representative DEM set, not just unit-tested in isolation.
- The hash function's exact bit-mixing does not need to match
  `boost::hash_value`'s bit-for-bit (nothing depends on absolute hash values,
  only on consistent hashing == equality), so this is a genuinely free choice,
  not a compatibility constraint.

## 5. Bottom line

**Yes, the 660 MB Boost dependency can be removed entirely.** It exists
solely to provide `boost::dynamic_bitset` and `boost::hash_value`
(`tesseract.h:18`, `tesseract.cc:18`, `visualization.h:4`) for one
narrowly-scoped purpose: dynamically-sized, word-packed detector/syndrome
bit-vectors with fast bitwise algebra and fast hashing. No other Boost
facility is used anywhere in `deps/tesseract/src/` or the rest of qudec.

What it takes:
1. Write a local `PackedBitset` (or similarly named) class covering exactly
   the operation surface in §2 (~120-150 lines, sketch in §4).
2. Replace `boost::dynamic_bitset<>` with it at the 5 declaration/parameter
   sites (`tesseract.h:164`; `tesseract.cc:350-351, 359, 523, 551, 563`;
   `visualization.h:16`; `visualization.cc:38`) and update every call site
   that uses `operator[]=`/`find_first`/`find_next`/`|=`/`&= ~` to the new
   method names.
3. Replace the `std::hash<boost::dynamic_bitset<>>` specialization
   (`tesseract.cc:85-94`) with an equivalent for the new type, or an explicit
   hasher functor on the `unordered_set` (`tesseract.cc:521`).
4. Remove `#include <boost/dynamic_bitset.hpp>` / `<boost/functional/hash.hpp>`
   from `tesseract.h`, `tesseract.cc`, `visualization.h`.
5. Delete the `boost_headers` `FetchContent` blocks and `Boost::boost`/
   `boost_headers` link lines from `CMakeLists.txt:60-82,126` and
   `deps/tesseract/CMakeLists.txt:45-55,94,114`, and the two `@boost//:dynamic_bitset`
   Bazel targets in `deps/tesseract/src/BUILD:185,209` if the Bazel build is
   still maintained (it is unclear from this audit whether the Bazel path is
   actually built by qudec's CMake-based flow — **not verified, flagged as
   an open question** for whoever does the mechanical edit).
6. Re-run/validate against `tesseract.test.cc` (if still present after the
   concurrent cleanup) plus a broader decode-output regression check, given
   the tail-bit-masking risk called out in §4.

Do **not** take the shortcut of `std::vector<bool>`, `std::vector<char>`, or
`std::bitset<N>` — each is disqualified for a documented, specific reason
(§3), not merely "less good."

## Open questions / not determined by this audit

- Whether the Bazel `BUILD` file (`deps/tesseract/src/BUILD:185,209`) is
  actually part of qudec's build (the CMake path looks canonical based on
  `CMakeLists.txt`/`deps/tesseract/CMakeLists.txt`, but this was not
  independently confirmed by checking a Bazel invocation anywhere in the
  qudec repo tooling).
- Exact current value(s) of `num_detectors`/`num_errors` for the code
  families qudec exercises (needed to turn the §3/§4 memory-footprint
  estimates from order-of-magnitude into exact numbers) — not looked up here;
  the 1000-detector example in §3/§4 is illustrative, not measured.
- Whether `tesseract.test.cc` / `tesseract_trellis.test.cc` survive the
  concurrent dead-code-deletion pass mentioned in this task's instructions —
  not re-checked after the initial file listing in §0.
