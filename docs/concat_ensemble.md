<!-- Author: OpenAI Codex. -->
# Concatenated MWPM first-stage ensembling

`qudec concat_mwpm` supports `--concat-ensemble K`. The default is `1`, which
preserves the baseline decoder. `K` counts one unperturbed matcher and `K-1` perturbed matchers for each
of the three ignored colors. For example:

```sh
./qudec concat_mwpm 9 --code cc666 --code-capacity -p 0.03 \
    --concat-ensemble 20 --concat-seed 20261024 --seed 42 \
    --mc-max-samples 100000 --mc-stop-at-errors 100000
```

Use `--phenomenological` for ideal check extraction with data and measurement
noise; its default number of rounds is the code distance. Omitting the ideal
noise selectors uses the existing SI1000 circuit model. Ensembling is available
for both cc666 and cc488.

The optional `--concat-boundary-gap` flag also evaluates the four first-stage
parity sectors of the two spatial boundaries for each ignored color. This is
available with ensembling and selects the minimum-weight final correction over
all feasible color, ensemble-member, and boundary-sector candidates. See
[`concat_boundary_gap.md`](concat_boundary_gap.md) for the parity construction.

`--concat-seed` controls the fixed graph ensemble independently of the sampling
`--seed`. Its default is `20261024`, matching the validated exploratory decoder.
It is shared across MPI ranks. The CLI requires positive `K` and a nonnegative
ensemble seed.

## Algorithm

Only the first matching stage changes. The unperturbed graph is retained. For
successive additional members, each compressed first-stage edge probability is
multiplied independently by `exp(N(0,sigma^2))`, alternating `sigma=ln(2)` and
`sigma=ln(4)`. Probabilities are clipped to `[1e-12,0.499999999]`, retaining
positive log-odds weights and the original graph topology. Different colors have
separate fixed random streams. Increasing K retains the same graph prefix.
This uses the probability perturbations described for Libra, without its
matching synthesis or degeneracy calculation.

For K greater than one, equal selected virtual-detector sets are coalesced per
ignored color across members. Each unique input is sent through the original
second-stage matching once. The first representative
retains its detector insertion order, so coalescing does not introduce new
PyMatching tie choices. All unperturbed candidates precede perturbed ones.

Second-stage graphs and final correction scores use the original model
probabilities. The minimum original negative-log correction weight wins; equal scores keep
the first candidate. Every correction must reproduce the complete spacetime
syndrome. Ensembling does not flatten measurement detectors or change flag
weighting. Circuit-level flag observations select representatives and update
one persistent matching ensemble in place. Graph objects are retained across
flag patterns; memory scales with K and the union of possible graph edges.
The existing seeded perturbation sequence is preserved as active projections
change. See [in-place updates](concat_inplace.md).

## C++ API

`dec::ConcatMWPM::config_type` adds `ensemble_size` (default 1) and
`ensemble_seed` (default 20261024). `dec::concat::CircuitDecoder` and the matching
core accept `EnsembleConfig`; the latter also offers `spread_scale`, mainly for
controlled tests. `spread_scale=0` makes all graph members identical. It has no
separate CLI option.

Inspected candidates expose `member` and canonical `virtual_detectors` fields.
The matching core reports first-stage calls, unique second-stage calls,
coalesced duplicates and infeasible candidates for its last decode.
Instances own mutable PyMatching workspaces and should be used by one worker.

The capacity and phenomenological benchmarks motivating K=20 are exploratory
results, not an exhaustive correction guarantee or a claim that larger ensembles
always improve logical accuracy. They are retained in `.agentic_artifacts/concat_ensemble/`.

With `config_type.debug=true` or CLI `--debug`, `dec::Result::debug_info` reports
first-stage matching attempts (including coalesced duplicates), second-stage
paths, candidate weights, the winning color/member/sector, and matching edge
weights; flag metadata is omitted. See [decoder_debug.md](decoder_debug.md). Debugging defaults to off.
