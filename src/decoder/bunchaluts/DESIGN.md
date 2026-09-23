# bunchaluts: a bunch-of-LUTs decoder

Design notes, written for an agent picking this up cold. Status: design settled,
implementation barely started. Code lives in `src/decoder/bunchaluts/`.

Wider context: this is a candidate **L1** decoder for runahead decoding — see
`docs/runahead-decoding.md`.

## Goal

Build a lookup-table decoder that scales to high code distance. The target is to
design a LUT such that the probability of encountering a syndrome the table
cannot decode falls below some target logical error rate.

The obstacle is table size: a monolithic LUT over all syndromes is hopeless. The
central idea is to use *a bunch of* LUTs — one for weight-1 errors, one for
weight-2, and so on up to some `W_max` — where each table's size grows with the
error weight it covers. Table contents are stored relative to a **base
detector**, which is what makes the tables independent of code distance.

Current scope is the **toric code**, chosen to avoid spatial boundaries so that
the base detector is translation-invariant.

**Weight is counted in DEM mechanisms, not physical faults.** A DEM mechanism
already encapsulates multiple physical faults.

## The cube

For a weight-k error cluster, the relevant region is the **cube**: the set of
detectors reachable from the base detector by a BFS of depth k in the decoding
graph. It is called a cube because it looks like one on a surface-code decoding
graph.

**Radius is k, not 2k.** A weight-k connected cluster anchored at a base
detector has all of its lit detectors within k graph steps of the base, so
radius k suffices to contain the syndrome. No extra halo radius is needed
because the decode order runs from larger cubes down to smaller ones (see
below), and because the halo condition is enforced separately as a
thickness-1 check rather than by inflating the cube.

Note: `hg::extract_error_cube()` in `hypergraph.tpp` currently documents itself
as taking a BFS depth of `2*error_count`. Per the above, the intended depth is
`error_count`. The docstring is stale.

The cube is a **build-time** construct — the region enumerated over when
generating table entries. It is *not* the match window at decode time. Keep
these separate in the code.

## Decode algorithm

Tile the input syndrome and slide a window across it.

1. Start with the largest LUT (largest cube, highest weight).
2. Scan across the syndrome looking for a match. On a match, decode that window
   and consume those syndrome bits.
3. Once the whole syndrome has been scanned, if lit bits remain, drop to the
   next-largest LUT and repeat.
4. Continue down to the weight-1 LUT. If lit syndrome bits remain when the
   ladder is exhausted, that is a decode **failure**.

### Match predicate: exact, with a thickness-1 halo

A window matches an entry when:

- the lit detectors inside the window are **exactly** the entry's support
  (exact match, not subset match), and
- **no lit detector is adjacent to the matched support** (halo of thickness 1).

The halo thickness is 1, *not* the full cube radius. This was the single most
important correction made during design discussion. Getting it wrong in either
direction breaks the decoder:

- **Subset matching with no halo** lets you commit to a truncated piece of a
  larger cluster.
- **Requiring the entire cube to be quiet** is far too strict. Two independent
  weight-1 errors whose base detectors sit at distance 2 would each be
  invalidated by the other, and neither would decode. At toric d=11, p=1e-3
  (~1.6e4 mechanisms, ~16 firing per shot, ~30 mechanisms per radius-2
  neighborhood) that situation arises in roughly 20% of shots — ten orders of
  magnitude above any usable LER.

With thickness 1, the pair above decodes correctly as two independent weight-1
events, while a genuinely interacting neighbor lights something adjacent to the
support, fails the match, and correctly falls through to the larger LUT. The
predicate draws the line exactly where interaction begins.

### Key consequence: thickness-1 halo == connected-component decomposition

"No lit detector adjacent to the matched support" is equivalent to "the matched
support is a full connected component of the syndrome." The tiling scheme and an
explicit connected-component decomposition are therefore the same algorithm, and
the tiling version gets the component version's properties without needing a
separate clustering pass.

This is what makes the table sizes work:

- The LUTs need **connected patterns only**. No disconnected/"random cluster"
  entries are required.
- Enumeration cost is therefore the lattice-animal count, ~c^k with c ≈ 30, not
  C(mechanisms-in-cube, k). At k=5 that is ~2e7 before symmetry reduction —
  large but tractable, whereas the binomial is hopeless.

The halo choice is what makes connected-only enumeration *valid*, and
connected-only enumeration is what makes large k *reachable*. The two decisions
hold each other up.

## Building the tables

### Order: small to large

Build the weight-1 LUT first, then weight-2, and so on. This is both the
deduplication order and the right order for the probability comparison.

### Deduplication across tables

Local syndrome patterns are degenerate — the same lit set can arise from a
weight-2 error and from a weight-4 error. If both tables hold the key, the
largest-first decode order would pick the less likely explanation.

So: when generating a candidate entry, check the already-built smaller LUTs for
a match. If found, discard it; do not store. Each distinct key ends up in the
lowest-weight LUT that can explain it.

### Deduplication within a table

If two generated errors of the same weight produce the same key, retain the one
with the highest error probability.

### Degeneracy is logically harmless when 2k < d

Two competing explanations of the same key each have weight ≤ k, so their
symmetric difference is a cycle of length ≤ 2k. A homologically nontrivial cycle
requires length ≥ d. So as long as **2k < d**, every explanation of a given key
induces the *same* logical action, and the highest-probability tiebreak above is
just picking a representative — it cannot cost you a logical error.

(An earlier version of this argument used the cube diameter and concluded
4k < d. That is too pessimistic; explanations are weight-bounded, so 2k < d is
the real condition.)

### Frame flips must not be stored in the entry

The syndrome pattern is translation-invariant but the logical action is **not** —
whether an error crosses the logical cut depends on absolute position. The base
detector is not representative in this respect.

Therefore an entry stores the **relative mechanism support** (mechanism offsets
from the base detector), not frame flips. At decode time, map the offsets to
concrete mechanism ids via the translation and XOR their actual `frame_flips`
out of the DEM (already carried in `hg::basic_edata`). One extra indirection per
correction; exact; keeps the table position-independent.

### Enumerate, don't sample

`build_lut` in `builder.h` currently takes an `RNG&&`, suggesting sampled error
configurations. Prefer exhaustive enumeration of connected clusters per weight,
using the DEM probabilities analytically for the degeneracy tiebreak. Coverage
is the entire correctness argument of this decoder, and sampling cannot certify
that a cluster shape is genuinely absent rather than merely unlucky.

## Error budget and scaling

Two failure sources:

1. **Non-coverage** — a connected cluster of weight > `W_max` appears somewhere.
   This dominates. Rate ≈ N·(c·p)^(W_max+1).
2. **Logical ambiguity** — eliminated by the 2k < d condition above.

Rough numbers at p=1e-3, where c·p ≈ 0.03 and ~16 mechanisms fire per shot:

| W_max | residual |
|-------|----------|
| 3     | ~4e-4    |
| 4     | ~1e-5    |
| 5     | ~4e-7    |

So `W_max` ≈ 5 for a ~1e-7 target, which needs d > 10 to satisfy 2k < d. d=11 is
already fine.

**The scaling claim the whole approach rests on:** total table size depends only
on cube size and `W_max`, not on d. Since failure goes like N·(c·p)^(W_max+1)
with N ~ d³, `W_max` only needs to grow like log(N/target)/log(1/(c·p)) —
logarithmic growth in table count against polynomial growth in the code. Worth
confirming or killing early with experiments.

## Important framing

At p=1e-3 with d=11 you expect ~16 mechanisms to fire per shot. You are never
decoding "a weight-k error" — you are decoding a **dilute gas of clusters**.
`W_max` bounds individual *cluster* weight, not total shot weight. The statistic
that governs the design is the cluster-size distribution.

## Open items

- **Temporal boundaries.** The toric code has no spatial boundary, but round 0
  and the final round still have different detector structure, so the base
  detector must sit in an interior round. Clusters whose earliest detector is in
  round 0 have no valid interior anchor. Likely resolution: run as a streaming
  sliding window in time with a commit region, so temporal edges are always
  covered by an adjacent window. Not yet settled.
- **Measurement not yet done.** A Monte Carlo that histograms (a) connected
  component sizes and (b) separations between nearest components, at the target
  p, would pin down `W_max` and the required cube radius empirically. Nothing
  here has been validated against simulation yet.
- **Symmetry reduction.** Toric-code point-group symmetry should shrink the
  tables by roughly 8x, at the cost of un-canonicalizing on lookup. Unexplored.

## Code pointers

- `hypergraph.h` / `.tpp` / `.cpp` — `Hypergraph<VDATA, EDATA, K>` template.
  `K=2` uses a vertex→(vertex,edge) hashmap for adjacency; `K>2` falls back to an
  interaction graph. `hg::basic_vdata` / `hg::basic_edata` are the DEM-backed
  payloads; `basic_edata` carries `error_probability` and `frame_flips`.
- `hg::from_toric_code_dem()` — builds the graph from a `stim::DetectorErrorModel`.
- `hg::extract_error_cube()` — BFS cube extraction. Stale docstring re: radius.
- `builder.h` — stub. Has a syntax error at line 20: `size_t error_max. RNG&&`
  uses a period instead of a comma.
