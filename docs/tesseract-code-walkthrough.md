# Tesseract decoder: code walkthrough

Written for an agent picking this up cold. Describes the vendored source at
`deps/tesseract/src/` (from `https://github.com/quantumlib/tesseract-decoder`,
`main` as of a shallow clone on 2026-09-23, git metadata stripped;
`_version.py` reports `0.1.1`. The upstream commit is *not* recoverable from
the vendored tree — re-clone if you need to pin it). Read 2026-09-23. This is
a **read-only trace of the actual decode path**, correlated
against the two source papers in `docs/papers/`:

- `tesseract-2503.10988.txt` — original paper (Beni, Higgott, Shutty), Algorithms 1-3.
- `accel-tesseract-2602.02985.txt` — follow-up systems paper (Grbic, Beni, Shutty),
  four `get_detcost` optimizations.
- `docs/papers/NOTES.md` — prior summary. Its Tesseract section is **accurate**
  on everything checked here, including the division-elimination point in §4
  below, which it already called out. No contradictions found between the code
  and either paper's *description of the algorithm* (the pruning predicate,
  A* heuristic, beam, pqlimit, etc. all match). The two things that don't
  appear in either paper are the sparsify pre-filter (§5) and the fact that
  current `get_detcost` has gone beyond Listing 2 (§4).

All line numbers below were read directly from the files; they are not
inferred.

---

## 1. The decode path, end to end

### Entry points

`TesseractDecoder` (`tesseract.h:93`) implements the `Decoder` interface
(`decoder.h:52`), whose only virtual method is
`decode_result(detections) -> DecodeResult` (`decoder.h:55`). Two overloads
exist in `tesseract.h:117-122`; the version taking `ErrorProbabilityUpdate`s
temporarily reweights errors (used for multipass correlated decoding, see §5)
and restores original weights before returning, via `ScopeExit`
(`tesseract.cc:29-41`, RAII rollback, used at `tesseract.cc:812`).

`decode_result` (`tesseract.cc:793-803`) is a thin wrapper: it calls
`decode_to_errors(detections)`, then reads out `predicted_errors_buffer`,
`low_confidence_flag`, and `cost_from_errors(...)` into a `DecodeResult`.

### Construction (`TesseractDecoder::TesseractDecoder`, `tesseract.cc:168-207`)

Per-decoder (amortized once, shared across shots) setup:

1. `config.dem = common::flatten(config.dem)` (`tesseract.cc:172`) — expands
   `REPEAT`/`SHIFT_DETECTORS` if present (`common.cc:171-173`).
2. Builds `dem_error_map`, identity at first (`tesseract.cc:174-175`).
3. If `merge_errors` (default true, `tesseract.h:48`):
   `common::merge_indistinguishable_errors` (`common.cc:175-230`) collapses
   DEM_ERROR instructions with identical `Symptom` (same detector set +
   observable set) into one, combining probabilities via
   `merge_weights` (`common.cc:154-159`, the log-likelihood-ratio combination
   formula for independent mechanisms with the same effect). Chains the index
   map (`chain_error_maps`, `common.cc:264-270`).
4. `common::remove_zero_probability_errors` (`common.cc:232-262`) drops
   probability-0 DEM_ERROR instructions, chains the map again.
5. `dem_error_to_error` / `error_to_dem_error` (`tesseract.h:137-138`) are the
   two directions of the original-DEM-index ↔ retained-decoder-error-index
   mapping; `invert_error_map` (`common.cc:272-283`) builds the inverse,
   taking the *first* original index that maps to each retained index.
6. `resolve_detector_orders` (`tesseract.cc:190`, defined `utils.cc:363-380`)
   materializes every requested `DetectorOrder` against the concrete DEM (see
   §3, ensemble reordering).
7. `errors = get_errors_from_dem(config.dem)` (`utils.cc:441-449`) — one
   `common::Error` per retained DEM_ERROR instruction. `common::Error`'s
   constructor (`common.cc:65-101`) computes
   `likelihood_cost = -log(p/(1-p))` and canonicalizes the detector/observable
   sets by XOR-into-a-`std::set` (so repeated targets in one DEM_ERROR cancel).
8. `initialize_structures(num_detectors)` (`tesseract.cc:321-412`, §6 below)
   builds `d2e`, `edets`, `eneighbors`, `error_costs`, and (if
   `sparsify_errors`) the mandatory/optional error partition.

Everything built here is **read-only for the lifetime of the decoder** except
under `decode_result(detections, probability_updates)`, which temporarily
mutates `errors[i].likelihood_cost` and the sort order of the affected rows
of `d2e`, then restores both (`apply_error_probability_updates` /
`restore_error_probabilities`, `tesseract.cc:281-319`).

### Per-shot setup and outer loop (`decode_to_errors`, `tesseract.cc:414-471`)

`decode_to_errors(detections)` is the per-shot entry point. Per-shot mutable
state is `predicted_errors_buffer`, `low_confidence_flag`, and (if sparsify is
on) `sparse_d2e` / `sparse_error_active` — all decoder-instance fields, not
thread-local statics, so one `TesseractDecoder` instance is not safely usable
from two threads concurrently on different shots (there is one shared,
mutated-per-call state block per instance).

If `config.sparsify_errors`, `build_sparse_d2e(detections)` runs first
(`tesseract.cc:420-422`, see §5) and `active_d2e` becomes an alias for
`sparse_d2e` instead of `d2e` (`tesseract.cc:423`) for the rest of the shot.

The outer loop is the **ensemble** over `config.detector_orders` (§3.
"Ensemble Reordering" in the paper). Two modes:

- **Beam climbing off** (default, `config.beam_climbing = false`,
  `tesseract.h:44`): `tesseract.cc:453-468` iterates `order_index` over every
  entry in `config.detector_orders`, calling
  `decode_to_errors_with_graph(detections, order_index, config.det_beam,
  active_d2e)` for each, and keeping the lowest-cost non-low-confidence
  result.
- **Beam climbing on**: `tesseract.cc:431-452` iterates a fixed number of
  trials (`max(det_beam+1, num_orders)`), cycling `beam` through
  `{0..det_beam}` and `order_index` through the order list independently
  (mod their respective sizes) each trial — this is exactly the paper's
  description ("try once for each beam value in `{0,...,B}}`, combined with
  ensemble reordering by using a different ordering for each beam setting",
  `tesseract-2503.10988.txt:213-216`).

In both modes the final answer is `best_errors` / `best_cost`
(`tesseract.cc:425-426`), i.e. the minimum-cost successful trial across the
ensemble — "if multiple valid solutions are obtained, output the minimum cost
decoding" (paper, line 211-212). `low_confidence_flag` is only cleared if at
least one trial converged (`tesseract.cc:470`).

### The A* search itself (`decode_to_errors_with_graph`, `tesseract.cc:504-744`)

One call = one A* search with one fixed detector order and one fixed beam
value.

**Setup** (`tesseract.cc:507-554`):
- `detector_at_position = config.detector_orders.at(order_index).get_order()`
  — the absolute detector ordering used by the pruning predicate (Algorithm 1
  needs this).
- `error_chain_arena` is cleared and (bounded) reserved to `pqlimit`
  (`tesseract.cc:512-518`) — this is the **backing store for search-tree
  nodes' parent pointers**, described below.
- `initial_detectors` (a `boost::dynamic_bitset<>`, one bit per detector) is
  set from the shot's `detections` (`tesseract.cc:523,526-536`); this is the
  residual syndrome `x = S ⊕ D(F)` at the root, i.e. `F = ∅` so `x = S`.
  `initial_detector_cost_tuples[ei].detectors_count` is simultaneously
  populated: for every fired detector `d`, every error `ei` incident to `d`
  gets its count bumped — this is the per-error "how many currently-fired
  detectors does this error touch" count that `get_detcost` needs (§6).
- `initial_cost = Σ_d get_detcost(d, ...)` over the fired detectors
  (`tesseract.cc:538-541`) — this is exactly Algorithm 2's `h(∅)`. If it's
  `INF` (no error can explain some fired detector at all — an inconsistent or
  malformed DEM/shot), the shot is immediately declared low-confidence
  (`tesseract.cc:543-546`).
- `min_num_dets = |detections|`, `max_num_dets = min_num_dets + detector_beam`
  — the beam-cutoff bookkeeping (`tesseract.cc:548-549`; see §3, "beam
  cutoff").
- Root node pushed: `pq.push({initial_cost, min_num_dets, depth=0,
  error_chain_idx=-1})` (`tesseract.cc:554`).

**The priority queue and node representation.** `Node` (`tesseract.h:65-75`)
has four fields: `cost` (the A* priority, `f = g + h`),
`num_dets` (residual detector-set popcount `r(F)`, used for the beam and as a
tiebreaker), `depth` (path length, used to size the output buffer), and
`error_chain_idx` (an index into `error_chain_arena`, or `-1` at the root).
`Node::operator>` (`tesseract.cc:131-133`) orders by `cost` ascending, ties
broken by *larger* `num_dets` first (`std::greater` + this comparator makes
the underlying `std::priority_queue` a min-heap on `(cost, -num_dets)`). The
queue itself is `std::priority_queue<Node, std::vector<Node>, std::greater<Node>>`
(`tesseract.cc:520`) — a **binary heap over `sizeof(Node)`-byte value types**,
not over pointers; nodes are small and copied by value.

**What state is carried per node.** A `Node` does *not* carry the residual
syndrome bitset, the per-error blocked flags, or the per-error fired-detector
counts — those are **reconstructed on pop** by walking the node's error chain
back to the root:

- `flip_detectors_and_block_errors` (`tesseract.cc:473-493`) walks
  `error_chain_arena` from `error_chain_idx` up through `parent_idx` chains,
  XOR-ing each ancestor error's detector set into a fresh copy of
  `initial_detectors` (`tesseract.cc:488-490`) and marking every error at or
  before that ancestor's `min_detector` position in `active_d2e[min_detector]`
  as blocked in a fresh `detector_cost_tuples` array (`tesseract.cc:483-486`).
  This reconstruction is **O(depth)** per pop — the arena is a compact
  append-only structure (`common::ErrorChainNode { error_index, min_detector,
  parent_idx }`, `common.h:61-65`), effectively a persistent/functional
  singly-linked list shared across the whole search tree, so each `Node` is a
  cheap `O(1)`-sized handle into that shared structure rather than an O(number
  of errors) or O(number of detectors) copy. This is the standard trick for
  making A* over an implicit exponential graph tractable in memory: **per-node
  incremental state is O(1) words; full state is reconstructed lazily and only
  when a node is actually expanded**, not when it's merely queued.

**Main loop** (`tesseract.cc:557-736`):

1. Pop the top node. If `node.num_dets > max_num_dets`, skip it
   (`tesseract.cc:561`) — a node can become beam-stale between being pushed
   and being popped, if `max_num_dets` shrank in the meantime (see beam
   climbing/tightening below).
2. Reconstruct `detectors` (residual syndrome bitset) and
   `detector_cost_tuples` (per-error blocked+count) via
   `flip_detectors_and_block_errors` (`tesseract.cc:563-566`).
3. **Termination / EXIT detection**: `if (node.num_dets == 0)`
   (`tesseract.cc:568`) — the residual syndrome is empty, i.e. `D(F) = S`, i.e.
   `F ∈ EXIT`. Because `h` is admissible and this is the first `EXIT` node
   popped off a consistent-cost priority queue, this is the A* optimality
   argument: the search returns the reconstructed error chain
   (`tesseract.cc:592-599`, walking `error_chain_arena` again, this time to
   fill `predicted_errors_buffer` with `error_to_dem_error[...]` indices in
   root-to-leaf order) and returns immediately.
4. **No-revisit-detections check** (`tesseract.cc:602-603`): if enabled and
   this exact `detectors` bitset was already visited at this `num_dets`
   level, skip expanding (Algorithm doesn't apply — this state's already been
   explored via a cheaper-or-equal-cost path, since it's popped later).
   `visited_detectors` is `unordered_map<size_t num_dets, unordered_set<dynamic_bitset<>>>`
   (`tesseract.cc:521`); dynamic_bitset hashing uses the `std::hash`
   specialization at `tesseract.cc:85-94` which delegates to
   `boost::hash_value` (§4, optimization 4).
5. **Beam tightening** (`tesseract.cc:630-638`): if this node's `num_dets` is
   a new minimum, `min_num_dets` drops, `max_num_dets` is retightened to
   `min(max_num_dets, min_num_dets + detector_beam)`, and stale
   `visited_detectors` buckets above the new `max_num_dets` are cleared
   (`tesseract.cc:632-635`) — a memory-bounding step, not just correctness.
6. **Pick the branching detector**: `min_detector` is the first detector in
   `detector_at_position` order (the resolved `DetectorOrder`) that is
   currently fired (`tesseract.cc:649-656`) — this *is* `dmin :=
   Minimum(x)` in Algorithm 1, using the ensemble-supplied total order on `D`.
7. **Expand children**: for each candidate error `ei` incident to
   `min_detector` in `active_d2e[min_detector]` (`tesseract.cc:661`), skip if
   `detector_cost_tuples[ei].error_blocked` (this is
   `GetForbiddenErrorsByPrecedence` — see §2). For each survivor:
   - incrementally update `next_detector_cost_tuples` for the *previous*
     survivor's detector footprint before moving on (`tesseract.cc:664-671`) —
     an optimization so the `detector_cost_tuples`/`next_detector_cost_tuples`
     deltas are computed once per sibling rather than recomputed from scratch;
   - flip `next_detectors` over `edets[ei]`, updating `next_num_dets` and
     `next_detector_cost_tuples[oei].detectors_count` for every detector's
     other incident errors (`tesseract.cc:680-687`) — this is the
     **incremental residual-syndrome and cost-tuple maintenance** that makes
     each expansion O(local degree) rather than O(all detectors);
   - beam cutoff check: `next_num_dets > max_num_dets` → skip
     (`tesseract.cc:689`);
   - no-revisit check on the child's syndrome (`tesseract.cc:691-693`);
   - recompute the heuristic delta only over the detectors that actually
     changed cost (the error's own detectors, plus `eneighbors[ei]` —
     detectors that share an error with `ei` and so might have had their
     minimum-cost witness change) (`tesseract.cc:695-713`), using a per-call
     `detector_cost_cache` to avoid recomputing `get_detcost` twice for the
     same detector within one expansion (`tesseract.cc:659,697-698,708-709`);
   - `if (next_cost == INF) continue` (`tesseract.cc:715`) — some residual
     detector has no unblocked witness error at all;
   - append a new `ErrorChainNode` to the arena and push the child
     (`tesseract.cc:718-724`);
   - **pqlimit check** (`tesseract.cc:727-733`): if `num_pq_pushed >
     config.pqlimit`, set `low_confidence_flag = true` and return
     immediately, mid-expansion. This is the heralded-failure exit and can
     fire in the middle of iterating a detector's candidate-error list.
8. If the queue empties without ever hitting `num_dets == 0`,
   `low_confidence_flag = true` (`tesseract.cc:737-743`) — this happens only
   if the beam prunes every path to `EXIT` (pqlimit not reached, but the tree
   under the beam constraint is exhausted).

**Runtime/memory character.** `error_chain_arena` grows monotonically with
`num_pq_pushed`, bounded above by `pqlimit` (reserved up front,
`tesseract.cc:516-518`); the priority queue itself is bounded by the same
count minus however many have been popped. Both are **per-shot, per-trial**
allocations — freed (via `.clear()`, `tesseract.cc:512`) and rebuilt at the
top of every `decode_to_errors_with_graph` call, i.e. once per (order, beam)
pair per shot. Runtime is data-dependent and unbounded above the `pqlimit`
early-exit; the low-confidence outcome is exactly the mechanism that turns
"unbounded" into "bounded but sometimes wrong."

---

## 2. Correlation with the paper's Algorithms 1-3

**Algorithm 1 (Pruning predicate `P_T(F, F')`)** — implemented, but as an
*implicit* filter rather than a literal boolean-returning function. The paper
defines two component predicates:

- `GetForbiddenErrorsByPrecedence(F)` — implemented via the *blocking* half of
  `flip_detectors_and_block_errors` (`tesseract.cc:483-486`): every error at
  or before the current one in `active_d2e[min_detector]`'s sort order, for
  every ancestor on the path, is marked blocked. Then the expansion loop only
  considers `active_d2e[min_detector]` and skips blocked errors
  (`tesseract.cc:661-662`). This reproduces the paper's exactness-preserving
  precedence rule — the code never literally builds a `GetForbiddenErrors`
  set object, it just marks a boolean per error, which is equivalent.
- The "`e` incident to `dmin`" restriction — implemented directly:
  `min_detector` is computed once (`tesseract.cc:649-656`) and only
  `active_d2e[min_detector]` is iterated as candidate children
  (`tesseract.cc:661`), so no other detector's incident errors are ever
  considered as children of this node.
- `GetForbiddenErrorsAtMostTwo(F)` — **NOT implemented as a persistent
  filter** the way the paper's Algorithm 1 frames it. The current source has
  no code path that limits an error to at most 2 per detector; grep for
  "at most two" / "AtMostTwo" over `tesseract.{h,cc}` finds nothing. The
  acceleration paper's §3 still lists "At Most Two Errors Per Detector" as an
  existing optional heuristic (`accel-tesseract-2602.02985.txt:176-181`), but
  no config flag, no code branch. **This appears to have been removed from
  current `main`, or was never actually wired into `TesseractConfig`** — worth
  flagging explicitly, since both papers describe it as present. I did not
  find it anywhere in `deps/tesseract/src/`; I did not exhaustively check
  every historical git tag of the vendored dep, so I can't rule out it being
  gated behind something I haven't seen, but the current `TesseractConfig`
  (`tesseract.h:41-63`) has no such field.

**Algorithm 2 (heuristic `h(F)`)** — implemented directly as the sum over
fired detectors of `get_detcost(d, ...)`: at the root, `tesseract.cc:538-541`;
incrementally per-expansion, `tesseract.cc:695-713`. Matches the paper's
`h(F) = Σ_{d∈x} GetDetCost(x, J, d)` with `x = S ⊕ D(F)`.

**Algorithm 3 (`GetDetCost`)** — `get_detcost` (`tesseract.cc:135-166`,
two overloads at `tesseract.h:160-162`). Matches the paper's definition
`min_{e∈E(d), e∉J} w(e)/|x∩D(e)|` structurally, but the current
implementation has been substantially restructured for performance relative
to both the paper's Algorithm 3 pseudocode and the acceleration paper's own
Listing 2 — see §4.

---

## 3. The heuristics, individually

All are fields of `TesseractConfig` (`tesseract.h:41-63`) unless noted.

- **Beam cutoff** (`config.det_beam`, default `DEFAULT_DET_BEAM = 5`,
  `tesseract.h:36,43`). Enforced via `max_num_dets` bookkeeping,
  `tesseract.cc:548-549,561,637,689`. **On by default** (beam 5; the papers'
  benchmarks use 15/20, so the compiled-in default is more aggressive than
  what either paper reports results for).
- **pqlimit / low-confidence heralding** (`config.pqlimit`, default
  `DEFAULT_PQLIMIT = 200000`, `tesseract.h:37,49`). Enforced at
  `tesseract.cc:727-733`. **On by default** at the paper's "short beam" value.
  Setting it to `SIZE_MAX` disables the arena-reservation optimization
  (`tesseract.cc:516-518`) but the check itself (`num_pq_pushed >
  config.pqlimit`) is unconditional — there's no way to fully disable
  heralding short of passing `SIZE_MAX`.
- **Ensemble reordering** (`config.detector_orders`, default one
  `DetectorOrder()` = `Method::Index` unresolved seed 0, `tesseract.h:53`).
  Four `Method`s (`utils.h:54-64`): `BFS` (randomized BFS over the
  error-induced detector graph, `utils.cc:107-153`, graph built by
  `build_detector_graph`, `utils.cc:76-105`), `Index` (either forward or
  reversed detector-ID order, coin-flipped per resolve,
  `utils.cc:196-207` — **not** what the paper describes for its benchmarks,
  which use `Coordinate`), `Coordinate` (paper's method: sample
  `z ~ N(0,1)^t`, sort detectors by inner product with declared coordinates,
  detectors without coordinates go last, `utils.cc:155-194`), and `Literal`
  (caller-supplied permutation, `utils.cc:280-281`, used by
  `load_detector_orders` from a JSON file, `utils.cc:334-361`). **On by
  default** trivially (exactly one order is always resolved), but the
  *ensemble* effect (multiple orders) requires the caller to populate more
  than one entry in `config.detector_orders`; the compiled-in default has
  only one, so out of the box there is no ensembling.
- **Beam climbing** (`config.beam_climbing`, default `false`,
  `tesseract.h:44`). **Off by default.** When on, iterates beam values
  `0..det_beam` paired cyclically with `detector_orders`
  (`tesseract.cc:431-452`).
- **No-revisit detections** (`config.no_revisit_dets`, default `true`,
  `tesseract.h:45`). **On by default.** Implemented via `visited_detectors`
  (an `unordered_map<num_dets, unordered_set<dynamic_bitset<>>>`,
  `tesseract.cc:521`), checked at pop time (`tesseract.cc:602-603`) and again
  before pushing each child (`tesseract.cc:691-693`) — i.e. it's checked
  *twice*, once as a late-arriving redundant-pop filter and once as an
  eager push filter, which is consistent (both operate on the same map) but
  means a node can be constructed, costed, and then discarded at push time
  without ever being queued.
- **Detection penalty** (`config.det_penalty`, default `0`,
  `tesseract.h:54`). Added unconditionally inside `get_detcost`'s return
  (`tesseract.cc:165`), so it's summed once per residual detector per
  heuristic evaluation — matches the paper's `c·|D(F)⊕S|` exactly, since `h(F)`
  sums `get_detcost` once per element of `x = S⊕D(F)`. **Off by default**
  (0).
- **At-most-two-errors-per-detector**: as noted in §2, **not present** in
  current `TesseractConfig` or anywhere in `tesseract.{h,cc}` — despite being
  documented as an existing optional heuristic in the acceleration paper.
  Flagged as a discrepancy between the paper and the vendored `main`.

---

## 4. The four acceleration-paper optimizations — confirmed present, and one goes further

All four are present in current source:

1. **`vector<bool>` → `vector<char>`.** No `std::vector<bool>` appears in
   `tesseract.h`/`tesseract.cc` at all (grep confirms). Blocked-error flags
   live in `DetectorCostTuple::error_blocked` (a `uint32_t` bitfield-ish flag,
   `tesseract.h:78`, not even a `char` — see optimization 2). Syndrome
   patterns (`detectors`, `next_detectors`, `initial_detectors`) are
   `boost::dynamic_bitset<>` (`tesseract.cc:523,551,563` etc.) — folded
   directly into optimization 4's data structure rather than being a
   `vector<char>` intermediate step, i.e. **current main skipped the
   vector<char> stage for syndrome patterns and went straight to
   dynamic_bitset**; the paper describes these as two separate, sequential
   optimizations (1 then 4) applied to the same data.
2. **SoA → AoS (`DetectorCostTuple`).** `struct DetectorCostTuple {
   uint32_t error_blocked; uint32_t detectors_count; }` (`tesseract.h:77-80`)
   — exactly the paper's structure (two `uint32_t` fields), and it's the sole
   type used for `detector_cost_tuples` / `next_detector_cost_tuples`
   throughout `tesseract.cc`. `ErrorCost { likelihood_cost, min_cost }`
   (`tesseract.h:82-85`) is the parallel AoS structure over the
   *precomputed, error-index-indexed* quantities (see optimization 3).
3. **Early-exit via precomputed lower bounds.** `error_costs[ei].min_cost` is
   precomputed once in `initialize_structures`
   (`tesseract.cc:334-339`, and recomputed on the affected subset in
   `update_internal_costs`, `tesseract.cc:264-279`) as
   `likelihood_cost / |symptom.detectors|` — i.e. cost assuming *all* of the
   error's detectors are currently fired, which is the maximum possible
   `detectors_count` and hence the minimum possible per-detector cost, giving
   a valid lower bound. `d2e[d]` is kept **sorted by `error_costs[idx].min_cost`
   ascending** for every detector (`tesseract.cc:341-346`, re-sorted on the
   affected rows in `update_internal_costs`, `tesseract.cc:273-278`), so the
   early-exit inside `get_detcost`'s loop sees candidates in non-decreasing
   lower-bound order.
4. **`boost::dynamic_bitset` + `boost::hash_value` for syndrome hashing.**
   Confirmed: `#include <boost/dynamic_bitset.hpp>` (`tesseract.h:18`),
   `#include <boost/functional/hash.hpp>` (`tesseract.cc:18`), and the
   `std::hash<boost::dynamic_bitset<>>` specialization delegating to
   `boost::hash_value` (`tesseract.cc:85-94`), used implicitly by
   `unordered_set<boost::dynamic_bitset<>>` inside `visited_detectors`
   (`tesseract.cc:521`).

### The division-elimination point (beyond the paper's own Listing 2)

Paper's Listing 2 (still divides in the hot loop, `accel-tesseract-2602.02985.txt`
lines ~382-400):

```cpp
for (size_t ei : d2e[d]) {
  ec = error_costs[ei];
  if (ec.min_cost >= min_cost) break;
  dct = detector_cost_tuples[ei];
  if (!dct.error_blocked) {
    double error_cost = ec.likelihood_cost / dct.detectors_count;
    min_cost = std::min(min_cost, error_cost);
  }
}
return min_cost + config.det_penalty;
```

Current `get_detcost` (`tesseract.cc:140-166`):

```cpp
double min_cost = INF;
uint32_t min_det_cost_det_count = std::numeric_limits<uint32_t>::max();
for (int ei : active_d2e[d]) {
  ec = error_costs[ei];
  if (ec.likelihood_cost * min_det_cost_det_count >=
      min_cost * errors[ei].symptom.detectors.size())
    break;
  dct = detector_cost_tuples[ei];
  if (!dct.error_blocked) {
    error_cost = ec.likelihood_cost;
    if (error_cost * min_det_cost_det_count < min_cost * dct.detectors_count) {
      min_cost = error_cost;
      min_det_cost_det_count = dct.detectors_count;
    }
  }
}
return (min_cost / min_det_cost_det_count) + config.det_penalty;
```

This confirms `NOTES.md`'s claim (lines 273-278) precisely. The transformation:
instead of maintaining `min_cost` as an already-divided ratio
`likelihood_cost / detectors_count`, it carries the winning `(likelihood_cost,
detectors_count)` pair as `(min_cost, min_det_cost_det_count)` — `min_cost`
here is *not yet divided*; it's the numerator of the best ratio seen so far.
Every comparison `a/b < c/d` (with `b, d > 0`) is replaced by the equivalent
cross-multiplication `a*d < c*b`, which is exact for the positive integers
`b, d` involved (`detectors_count` values) and avoids a floating-point divide
per candidate. Both:
- the early-exit test (line 151-152, comparing the current candidate's
  *precomputed* lower-bound numerator `ec.likelihood_cost` against the current
  best via cross-multiplication with `errors[ei].symptom.detectors.size()`,
  the error's *maximum* possible detector count — this is the lower-bound
  quantity from optimization 3, still valid as a break condition since it's
  monotonic in the same sorted order), and
- the update test (line 158, comparing the *actual current* `detectors_count`
  via cross-multiplication)

are done without division. The **single remaining divide is deferred to the
return statement** (`tesseract.cc:165`), so a call to `get_detcost` now does
at most one floating-point division total, regardless of how many candidates
it scans, versus one division per surviving (non-early-exited, non-blocked)
candidate in the paper's Listing 2. This is a valid transformation because
division is monotonic for positive divisors, and `detectors_count` values are
always positive when reached (a detector's candidate errors are, by
definition, incident to it, so `detectors_count ≥ 1` for any error scanned via
`active_d2e[d]`... **caveat**: this holds for the sum inside `get_detcost`
where `d ∈ active_d2e`'s domain, but I did not independently re-derive a proof
that `detectors_count` can never be `0` for an `ei` reached inside this loop
in every call site — it's plausible from the incremental-update logic in
`decode_to_errors_with_graph` that a currently-blocked error could transiently
have `detectors_count == 0`, but such entries are skipped via the
`error_blocked` check before the division would matter for `min_det_cost_det_count`
update. Marking this as inference, not verified by a formal argument or by
running the code.)

---

## 5. What exists in the code but in neither paper

### `sparse_d2e` / `sparsify_*` — a per-shot candidate-error pre-filter

Config fields: `sparsify_errors` (default `false`), `sparsify_base_degree`,
`sparsify_max_degree`, `sparsify_reactivate_limit` (`tesseract.h:57-60`).
State: `sparse_d2e`, `sparse_error_active`, `sparsify_mandatory_errors`,
`sparsify_optional_errors` (`tesseract.h:148-151`).

**Setup** (`initialize_structures`, `tesseract.cc:375-411`), done once at
decoder-construction time (shared, read-only after this point except
`sparse_error_active`/`sparse_d2e` which are rebuilt every shot):

- Every error is partitioned by its detector-symptom degree
  (`errors[ei].symptom.detectors.size()`):
  - degree `≤ sparsify_base_degree` → **mandatory**
    (`sparsify_mandatory_errors`, `tesseract.cc:400-403`) — always active,
    every shot.
  - degree `> sparsify_base_degree` and (`sparsify_max_degree == -1` or
    degree `≤ sparsify_max_degree`) → **optional**
    (`sparsify_optional_errors`, `tesseract.cc:404-407`) — candidates for
    per-shot reactivation.
  - degree `> sparsify_max_degree` (when `sparsify_max_degree ≥ 0`) → dropped
    entirely, forever, from every shot's candidate set. This is unconditional
    exclusion, not per-shot.
- `sparsify_reactivate_limit`, if `-1`, is auto-derived by
  `suggest_sparsify_reactivate_limit_capped` (`tesseract.cc:57-81`), an
  empirical formula: `limit ≈ round(exp((base_degree-2)*ln(4.5) - ln(3) +
  ln(num_detectors)))`, capped at the number of errors. This is a heuristic
  scaling rule with no derivation given in-source (no comment ties it to a
  formula in either paper) — **the magic constants `4.5` and `3` are
  unexplained in the source**; I could not find their derivation anywhere in
  `deps/tesseract/src/` or the papers. Treat as an empirically tuned default,
  not a principled bound.

**Per-shot** (`build_sparse_d2e`, `tesseract.cc:824-888`), called from
`decode_to_errors`/`decode_to_errors` before the search
(`tesseract.cc:420-422,497-499`):

1. All mandatory errors are marked active (`tesseract.cc:834-836`).
2. For every optional error, compute `overlap` = number of its detectors that
   are in this shot's fired set (`tesseract.cc:848-854`). Errors with zero
   overlap are dropped from consideration entirely for this shot
   (`tesseract.cc:855`, the `if (overlap > 0)` guard).
3. Remaining candidates are sorted by `(overlap desc, degree asc,
   likelihood_cost asc, error_index asc)` (`tesseract.cc:861-873`) — highest
   overlap-with-the-actual-syndrome first, then prefer lower-degree, cheaper,
   lower-index as tiebreaks.
4. The top `min(sparsify_reactivate_limit, candidates.size())` are marked
   active (`tesseract.cc:875-877`).
5. `sparse_d2e[d]` is rebuilt for every detector as the subsequence of
   `d2e[d]` (which is cost-sorted, §6) whose errors are active
   (`tesseract.cc:880-887`) — order is preserved from `d2e`, so
   `sparse_d2e[d]` is still cost-sorted, meaning the early-exit logic in
   `get_detcost` (§4) still works correctly over it.

**Is this exact, or heuristic?** It is **explicitly a heuristic, not an
exactness-preserving filter**. The mandatory/optional split by degree, and
especially step 2's outright removal of zero-overlap optional errors and step
4's hard cap at `sparsify_reactivate_limit`, can remove the true optimal error
from candidacy for a given shot: an optional-tier error with zero detectors in
common with the actual fired set by construction can never explain any
residual detector directly, so dropping it is *locally* safe *for that
error's role as a direct explanation of a detector* — but errors don't just
explain detectors independently; a higher-degree error could still be part of
the eventual minimum-cost `F` if some of its detectors get activated
transiently by other chosen errors during the search... **actually, on
inspection**, `sparse_d2e` is fixed for the whole shot (built once before the
search starts, from the *original* fired detections, not updated as the
search hypothesizes intermediate states), so any error excluded from
`sparse_d2e` is unavailable as a candidate for the *entire* A* search on this
shot, regardless of what intermediate residual syndromes arise during search.
This means the search can only construct chains from the surviving error set,
so if the true minimum-cost solution requires an excluded error, sparsify
will produce a suboptimal (or, in edge cases where no valid chain exists at
all within the sparse set, a `low_confidence` or outright failed) result. This
is stated nowhere in the code comments as a formal guarantee either way; the
above is my own analysis of the control flow, not a quoted design comment —
flagging as inference from reading the code, since no comment or paper states
whether this preserves exactness. **Given the mechanics, it clearly does
not preserve exactness in general** — it is a candidate-set pruning heuristic
analogous in spirit to the beam cutoff, but operating on the *error graph*
before the search starts rather than on *search nodes* during the search.
This heuristic does not appear anywhere in either paper; it is
undocumented in prose beyond the CLI help text in `tesseract_main.cc`.

### `src/tesseract_trellis.{h,cc}` — a separate beam-search decoder

`TesseractTrellisDecoder` (`tesseract_trellis.h:65-105`) is a **structurally
different algorithm**, not a variant of the A* decoder sharing code with it
(it doesn't include `tesseract.h`, and defines its own error/config/state
types). Based on reading `tesseract_trellis.h` and the top ~110 lines plus
grep of `tesseract_trellis.cc` (1294 lines total; not traced line-by-line per
the task's scope):

- It processes detectors in a **layered trellis**: `wide_layer_templates`
  (`tesseract_trellis.h:101`), one `TesseractTrellisWideLayerTemplate`
  (`tesseract_trellis.h:42-52`) per layer, each carrying
  `surviving_local_indices`, `current_active_detectors`, and
  `detcost_transition` data — this is a precomputed, DEM-derived structure
  (shared/read-only, built once), distinct from the per-shot search state.
- Per-shot state is a beam of `FixedWideStateEntry<Words>`
  (`tesseract_trellis.cc:66-71`): each entry is a fixed-width bit-packed
  `state_words` array (`std::array<uint64_t, Words>`, `Words` chosen at
  compile time up to `kMaxCompiledWideStateWords = 4`,
  `tesseract_trellis.cc:44`) plus `mass0`, `mass1` (probability mass for
  each of the two observable outcomes), `penalty`, and `score`.
- The per-layer loop is expand → collapse (bucket by state, merging mass for
  identical states — see `FixedWidePairBucket`, `tesseract_trellis.cc:76-78`)
  → truncate to `beam_width` via `keep_top_compiled_states`
  (`tesseract_trellis.cc:873-895`, using `std::nth_element`, i.e. a partial
  sort, not a full sort) → repeat for the next layer
  (`tesseract_trellis.cc:1006-1111` is the per-shot `decode_shot` method of
  the internal kernel).
- Ranking is configurable via `TesseractTrellisRankingMode`
  (`tesseract_trellis.h:29-33`): `MassOnly` (rank by accumulated probability
  mass alone), `FutureDetcostRanked` / `FutureActiveDetcostRanked` (also
  factor in a forward-looking `get_detcost`-style estimate,
  `future_detcost_scale` config knob, `tesseract_trellis.h:58`).

This is **a sum-product / forward-algorithm-style beam search that tracks
probability mass through a fixed layer sequence**, not a best-first search
over an implicit graph with a priority queue and an admissible heuristic. It
has no analog to `pqlimit`/A*-optimality; instead it has a fixed per-layer
`beam_width` (default `1024`, `tesseract_trellis.h:56`) that bounds *both*
runtime and memory deterministically per layer (`O(layers × beam_width ×
local expansion factor)`), which is a materially different runtime-bound
character from the A* decoder's `pqlimit`-bounded-but-otherwise-unbounded
runtime. Neither paper describes this decoder; it is not mentioned by name in
either paper's text.

### `src/multi_pass/` — `dem_decomposition`, `error_correlations`

Neither is a decoding algorithm on its own; both are **preprocessing/analysis
utilities feeding `MultiPassTesseractDecoder`**
(`multi_pass/multi_pass_tesseract_decoder.h`), which is itself a wrapper
around two `TesseractDecoder` instances run over a monolithic DEM split into
two components:

- `dem_decomposition.{h,cc}`: `prepare_two_component_dem`
  (`dem_decomposition.h:40-41`) validates a caller-supplied
  detector→component assignment (each detector labeled 0 or 1) and produces a
  `TwoComponentDem` (`dem_decomposition.h:25-32`) — each component's own
  sub-DEM plus a "decomposed" combined DEM where cross-component errors have
  been split. `decompose_errors_using_detector_assignment`
  (`dem_decomposition.h:50-51`) does the actual per-error splitting, rejecting
  detectorless error groups (ambiguous which component they belong to) and
  requiring existing Stim `^`-separated decomposition groups to already
  respect the assignment.
- `error_correlations.{h,cc}`: `collect_correlation_evidence`
  (`error_correlations.h:52`) walks an already-decomposed `TwoComponentDem`
  and tabulates, per component symptom, how often it co-occurs with error
  mechanisms that touch symptoms in *both* components
  (`CorrelationEvidence::paired_mechanism_probabilities`,
  `error_correlations.h:41-46`, explicitly documented as *not* the true joint
  probability — it "deliberately omits combinations of independent one-sided
  mechanisms"). `derive_reweight_probabilities`
  (`error_correlations.h:54-60`) turns this into a
  paired/marginal probability ratio, explicitly flagged in its own doc
  comment as "not, in general, an exact conditional probability" — a
  correlated-matching-style heuristic reweighting signal.
- `MultiPassTesseractDecoder` uses these to run one component decoder, feed
  its output as `ErrorProbabilityUpdate`s (the reweighting hook in
  `tesseract.h:119-122`, §1 above) into the other component's decoder, and
  repeat for `num_passes` (`multi_pass_tesseract_decoder.h:41`) — effectively
  a two-block Gauss-Seidel-style iterative correlated decoding scheme,
  scheduled either `Static` (both components every pass) or `Causal` (derived
  from a dependency graph) (`multi_pass_tesseract_decoder.h:31-34`).

Neither paper mentions multipass decoding; it does not appear in either PDF's
extracted text (grepped both for "multi_pass", "multipass", "component" —
no matches tied to this feature).

### `src/simplex.{h,cc}` — the integer-program baseline

`SimplexDecoder` (`simplex.h:43-90`) wraps `HighsModel`/`Highs`
(forward-declared, `simplex.h:24-26`) — **HiGHS**, an open-source
linear/mixed-integer programming solver. This is the "Integer Program
decoder" both papers use as the accuracy baseline/ground truth
(`tesseract-2503.10988.txt:225-236`, `accel-tesseract-2602.02985.txt:132-135,
216-218`). `init_ilp()` (`simplex.h:85`) presumably builds the MILP
formulation of the most-likely-error problem directly (one binary variable
per error, parity constraints per detector, minimize total weight) — I did
not read `simplex.cc`'s 420 lines in full, since the task scope only asks to
identify this file, not trace it.

---

## 6. Data structures and their cost

All of the following are built once in `initialize_structures`
(`tesseract.cc:321-412`) and are **read-only, shared across all shots and all
ensemble trials** for the lifetime of a `TesseractDecoder`, with the sole
exception of the rows touched by `apply_error_probability_updates` (rare,
multipass-only) and the sparsify-specific structures which are shot-scoped.

- **`d2e`** (`tesseract.h:147`, `vector<vector<int>>`, size `num_detectors`):
  detector → list of error indices incident to it. Built by iterating every
  error's `symptom.detectors` and appending to each detector's list
  (`tesseract.cc:325-330`). **Sorted per-detector by `error_costs[idx].min_cost`
  ascending** (`tesseract.cc:341-346`), which is what makes the `get_detcost`
  early-exit (§4, optimization 3) valid. Indexed by detector ID; accessed as
  `d2e[d]` — this is the array that `get_detcost` walks, and per the
  acceleration paper is the single largest source of cache misses in the
  whole decoder (the paper's entire optimization effort targets the two
  vectors — blocked flags and detector counts — that get **randomly gathered**
  at the error indices found while walking `d2e[d]`; `d2e[d]` itself is
  accessed sequentially, but what it's used to index into,
  `detector_cost_tuples[ei]` and `error_costs[ei]`, is a **random gather**
  keyed by whatever error indices happen to live in that detector's sorted
  list).
- **`sparse_d2e`** (`tesseract.h:148`): same shape as `d2e`, rebuilt every
  shot as a filtered subsequence (§5). Order-preserving from `d2e`, so
  remains cost-sorted.
- **`edets`** (`tesseract.h:154`, `vector<vector<int>>`, size `num_errors`):
  the reverse map — error index → its detector symptom (a copy of
  `errors[ei].symptom.detectors`, `tesseract.cc:326`). Indexed by error index;
  accessed sequentially inside `flip_detectors_and_block_errors` and the
  expansion loop wherever an error's own detector footprint needs to be
  walked (e.g. `tesseract.cc:488,680,695`).
- **`eneighbors`** (`tesseract.h:153`, `vector<vector<int>>`, size
  `num_errors`): for each error, the set of detectors that are touched by
  *some other* error sharing at least one detector with it, minus its own
  detectors (`tesseract.cc:348-373`, built via per-error `dynamic_bitset`
  unions over `d2e[d]` for each of its own detectors, then bitwise-andnot of
  its own footprint). This is what lets the heuristic-delta recomputation in
  the expansion loop (`tesseract.cc:706-713`) know exactly which *other*
  detectors' `get_detcost` might have changed as a side effect of adding
  error `ei`, without re-scanning every detector in the code. Precomputed
  once (potentially expensive at construction: `O(Σ_e degree(e) × avg d2e
  degree)` bitset unions), then a pure read-only lookup at decode time.
- **`error_costs`** (`tesseract.h:156`, `vector<ErrorCost>`, size
  `num_errors`): `{likelihood_cost, min_cost}` per error, `min_cost` being the
  precomputed lower bound from §4 optimization 3
  (`likelihood_cost / degree`). Indexed by error index; this is the array
  `get_detcost` reads sequentially (per detector's sorted `d2e[d]` list) to
  decide when to early-exit.
- **`detector_cost_tuples`** (a *local variable* inside
  `decode_to_errors_with_graph`, not a decoder field — `vector<DetectorCostTuple>`,
  size `num_errors`, reallocated fresh per popped node,
  `tesseract.cc:564,647,552`): `{error_blocked, detectors_count}` per error,
  for *this* node's reconstructed state. This is the array that's **randomly
  gathered** by error index while walking a detector's `d2e[d]`/`sparse_d2e[d]`
  list inside `get_detcost` — the AoS optimization (§4 opt 2) exists
  specifically because `error_blocked` and `detectors_count` used to live in
  two separate arrays both gathered at the same random index.
- **Sort order on `d2e`**: ascending `error_costs[idx].min_cost`, ties broken
  by index (`tesseract.cc:342-345,274-277`). This sort is what the early-exit
  in `get_detcost` depends on for correctness (an unsorted or wrongly-sorted
  `d2e[d]` would make the `break` on line 153 of `tesseract.cc` unsound,
  since it assumes no later candidate in the list can beat the current best).

**Gather vs. sequential summary**: `d2e[d]` and `edets[ei]` themselves are
walked sequentially (good cache behavior, small vectors, likely a few
detectors/errors each for LDPC-sparse DEMs). What's *looked up* while walking
them — `error_costs[ei]`, `detector_cost_tuples[ei]`, `errors[ei]`,
`edets[ei]` (nested), `eneighbors[ei]` — are keyed by the error/detector index
found in the outer sequential list, and those indices are not contiguous or
predictable, so those accesses are the random-gather component the
acceleration paper's profiling identified as dominant (70-90% of runtime for
bivariate-bicycle codes, per `accel-tesseract-2602.02985.txt:1164-1168`).

---

## 7. Per-instance vs. shared state — summary for the parallelism question

**Shared, read-only after construction** (one copy suffices across many
concurrent decode calls **on different shots**, given separate `Node`/queue/
arena state — but see the caveat below about sparsify and probability
updates):
`d2e`, `edets`, `eneighbors`, `error_costs`, `errors`, `dem_error_to_error`,
`error_to_dem_error`, `config.detector_orders` (post-resolution),
`sparsify_mandatory_errors`, `sparsify_optional_errors`.

**Per-shot, per-(order,beam)-trial mutable state, freed/rebuilt every call**:
`error_chain_arena`, the `pq`, `visited_detectors`, `initial_detectors`,
`detector_cost_tuples` (root and per-node), `sparse_d2e`/`sparse_error_active`
(shot-scoped, rebuilt by `build_sparse_d2e` every `decode_to_errors` call —
**this makes a `TesseractDecoder` instance with `sparsify_errors=true` NOT
safely usable across two shots concurrently on the same instance**, since
`sparse_d2e` is an instance field mutated per-shot, not a local variable, in
contrast to `detector_cost_tuples` which is stack-local per search call).
`predicted_errors_buffer` and `low_confidence_flag` are likewise instance
fields overwritten per shot (`tesseract.h:135-136`), so the whole
`TesseractDecoder` object is **not thread-safe for concurrent decode calls on
one instance**, sparsify or not — parallelism across shots requires one
`TesseractDecoder` instance per concurrent shot (each holding its own copy of
the shared read-only tables, unless the caller explicitly shares one decoder
object across sequential calls only). This matches `NOTES.md`'s framing: the
per-instance mutable state (arena + queue, data-dependent size up to
`pqlimit`) is the expensive-to-replicate part; the shared tables
(`d2e`/`edets`/`error_costs`/etc.) are the cheap-to-share part, but current
code structure ties both to the same C++ object rather than separating them,
so replication in practice means either duplicating the shared tables per
instance or restructuring the code to pass them by reference into a
per-shot-only mutable context.
