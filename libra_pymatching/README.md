# Libra with PyMatching assignment

The stable implementation is now part of the normal `qudeclib` library.
It uses `src/decoder/libra_pymatching.cpp`, PyMatching auxiliary assignment,
and exact subset dynamic programming for components with at most eight vertices.
The correlated matching, synthesis, degeneracy, and large-gap fallback are retained.

See [the integration guide](../docs/libra.md) for normal build commands, the
`qudec libra` CLI, corrected rare event sampling, and evaluation targets.
`libra_decoder` and `libra_pymatching_decoder` are compatibility aliases for
`qudeclib`; the historical CMake entry points forward to the normal project.
Only the stable implementation is compiled. Historical measurements remain in
[Libra-PyMatching.txt](../Libra-PyMatching.txt).

## Exact assignment reduction

The vertices in this auxiliary problem are already selected CSS graph edges.
Each vertex must be covered exactly once by an independent single channel
of cost `s_i`, or a correlated pair channel of cost `p_ij`. These are explicit
costs; the earlier Blossom V formulation needed no shortest path preprocessing.

Simply firing every vertex in a PyMatching graph would enforce odd degree,
which also permits degree three. Instead, quantize costs to the original
`1e-6` resolution, set `C = 1 + max(s_i)`, and add:

- A boundary edge of cost `C + s_i` at each vertex.
- A pair edge of cost `2*C + p_ij` for each useful correlated pair.

Fire every vertex. If an odd-degree solution has degree at least three at
vertex `v`, remove two of its pair edges `(v,a)` and `(v,b)`, and toggle the
boundary edges at `a` and `b`. Parities are preserved and the cost decreases
by at least `2*C - s_a - s_b > 0`. Consequently every optimum covers every
vertex once. Its added cost is the constant `n*C`, so the original optimum
is preserved. PyMatching receives doubled integer costs so its collision
times are integral. The helper rejects overflow of its 32-bit edge weights.

All optimal paths in this graph are direct. The compressed matched endpoints
therefore identify the selected channels without path recovery or a
`SearchFlooder`. Parallel pairs retain the cheapest channel; dominated pairs
are omitted. The quadratic clique of private dummy vertices is eliminated.

## Path recovery in the main decoder

PyMatching's Sparse Blossom solver avoids constructing an all-pairs shortest
path metric graph. Its ordinary observable-output decode also avoids recovery
of explicit error paths. Libra's correlated matching and synthesis need the
actual graph edges and hyperedges, so those main candidate decodes still use
PyMatching's edge-output API and its Dijkstra-style path recovery. Replacing
the auxiliary assignment solver does not remove that work.

## Logical-class probabilities

For syndrome `s`, error configuration `e`, detector map `D`, and logical parity
`L(e)`, the ideal objective is proportional to

```text
P(L = l | s) ∝ sum_{e: D(e)=s, L(e)=l} exp(-W(e))
W(e) = sum_{active channels i} log((1-p_i)/p_i).
```

Core Libra synthesizes improved representative configurations. The degeneracy
extension enabled here approximates each class's sum using retained cycles
around those representatives. It does not enumerate every configuration.
The resulting log masses are approximate, unnormalized scores; the
`matching_data.probability` field is also unnormalized. The original gap gate
can return correlated matching's prediction without running the ensemble.
