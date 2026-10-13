# Libra with PyMatching assignment

This build uses `src/decoder/libra_pymatching.cpp` as an alternative
implementation of the existing `dec::Libra` API. It replaces the auxiliary
Blossom V assignment solver with PyMatching. Assignments with at most eight
vertices use an exact subset dynamic program to avoid solver setup overhead.
All other matching synthesis, correlation handling, gap gating, and degeneracy
logic follows the original implementation.

Existing source, build files, and `Libra.txt` remain available. The new results
and algorithm explanation are in [Libra-PyMatching.txt](../Libra-PyMatching.txt).

## Build and test

```sh
cmake -S libra_pymatching -B .agentic_artifacts/libra_pymatching/build \
  -DCMAKE_BUILD_TYPE=Release -DENABLE_MPI=OFF
cmake --build .agentic_artifacts/libra_pymatching/build -j8 \
  --target libra_pymatching_evaluate libra_pymatching_splitting_evaluate \
           libra_pymatching_test libra_pymatching_assignment_test
ctest --test-dir .agentic_artifacts/libra_pymatching/build \
  -R '^libra_pymatching_' --output-on-failure
```

Use a C++23 compiler. Link `libra_pymatching_decoder` instead of `libra_decoder`;
both implement the same class and must not be linked into one executable.
The included project still builds and links its existing Blossom V library
for other decoders. This Libra implementation makes no Blossom V calls.

The evaluators accept the same arguments as the original
[Libra evaluators](../libra/README.md). The new splitting target enables the
corrected sampler's fast failure callback, local swaps, and systematic
resampling:

```sh
.agentic_artifacts/libra_pymatching/build/libra_pymatching_evaluate mc \
  9 64 10000000 20261001 .agentic_artifacts/libra_pymatching/example_mc.txt \
  --max-errors 20

taskset -c 0 .agentic_artifacts/libra_pymatching/build/libra_pymatching_evaluate latency \
  9 1 65536 20261001 .agentic_artifacts/libra_pymatching/example_latency.txt

.agentic_artifacts/libra_pymatching/build/libra_pymatching_splitting_evaluate \
  13 64 3000000 20261001 .agentic_artifacts/libra_pymatching/example_splitting.txt \
  0.8 8 1000 512 4
```

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
