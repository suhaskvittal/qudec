<!-- Contributor: OpenAI Codex. -->
# Seeded matching

`qudec seeded DISTANCE` runs `SeededMatching<D1, D2>`. The first decoder (`D1`)
is uncorrelated PyMatching or weight-aware union find. It reports its matched
detector pairs and the physical correction edges in each pair's path. The
second decoder (`D2`) is correlated PyMatching or Libra.

For correlated PyMatching, the seed's physical correction edges condition the
matching graph; the second decoder then matches the remaining syndrome. No gap
estimate is computed. Without postprocessing, the remaining syndrome is the
original one, and a PyMatching seed reproduces PyMatching's native correlated
two-pass decoder.

For Libra, the seed gap gates the expensive ensemble. A high gap returns the
seed result. A lower gap sends the remaining syndrome to Libra and forces its
ensemble path. The default is a 40 dB **cluster gap**. A PyMatching seed may
instead use the complementary gap, which requires a separate constrained
matching calculation and is substantially slower. Empty syndromes return the
seed result immediately.

`--seed-postprocess` enables pair retention. A seed pair is frozen only if its
selected physical path has no conditioned weight change and no individual
connection from either endpoint to another fired detector or relevant boundary
is strictly cheaper after conditioning. Both endpoints are removed together;
their observable parity is combined with the second decoder's prediction.
The path checks use canonical shortest routes in the original graph, scored
at conditioned weights. This is a heuristic filter and can change the logical
error rate.

Useful options:

| Option | Values / default |
| --- | --- |
| `--seed-decoder` | `pymatching` (default), `uf` |
| `--second-decoder` | `correlated_pymatching` (default), `libra` |
| `--seed-postprocess` | Off by default |
| `--seed-gap` | `cluster` (default), `complementary` with PyMatching seed |
| `--seed-gap-threshold-db` | `40` by default; only used with Libra |
| `--libra-ensemble-size` | `60` by default |
| `--uf-weight-scale` | `1024` by default |

For example:

```sh
qudec seeded 9 --seed-decoder pymatching --second-decoder correlated_pymatching --seed-postprocess
qudec seeded 9 --seed-decoder uf --second-decoder libra --seed-gap cluster --seed-gap-threshold-db 40
```

The SI1000 circuit for this decoder always includes opposing-basis detectors.
Printed statistics include seed and second-decoder call counts, retained pairs,
input and residual detection-event counts, full pipeline latency, and the
second decoder's latency weighted by its call fraction.
