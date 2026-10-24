<!-- Author: OpenAI Codex. -->
# Concatenated MWPM boundary-parity ensemble

Enable the first-stage boundary-parity extension with:

```sh
qudec concat_mwpm 7 --code cc666 -cap -p 0.01 --concat-boundary-gap \
    --mc-max-samples 100000 --mc-stop-at-errors 1000
qudec concat_mwpm 7 --code cc666 -pheno -p 0.01 --concat-boundary-gap \
    --mc-max-samples 100000 --mc-stop-at-errors 1000
```

`--concat-boundary-gap` replaces `--concat-permissions`. The old option is not
an alias. Ordinary three-color concatenated MWPM remains the default. The new
option is accepted only for `concat_mwpm` and applies to Monte Carlo and
rare-event sampling. Despite its name, it enumerates boundary-parity candidates;
it does **not** compute a numerical gap or use a confidence score for selection.

## First-stage parity sectors

For each ignored color, known spatial data hyperedges are completed with their
missing-color physical boundaries before projection. The two retained-color
boundaries are explicit parity vertices. On `!R` these are G and B. All physical
edges keep their probabilities, including the positive-weight G-boundary to
B-boundary edge supplied by a corner fault. No artificial zero-weight bridge is
added.

The decoder runs the four boundary syndromes: even/even, odd/even, even/odd,
and odd/odd. A fired boundary vertex requires an odd number of selected edges
incident on it; an unfired vertex requires an even number. Even allows zero,
two, four, etc.; odd allows one, three, five, etc. These are parity constraints,
not forbidden/allowed endpoints or capacity-one slots. For a pure spatial
code-capacity or native phenomenological model, only two sectors per color
have feasible total parity, normally yielding six candidates.

The first-stage solution is a selected-edge T-join. If converted to physical
matching paths, split at spatial boundary vertices: multiple paths can end
there, but a physical detector-to-detector path must not pass through a boundary
and return to the bulk. This changes neither selected edges nor parity. The
public `restricted_edges` field contains selected individual edges, not a
pairing of physical detector endpoints.

## Original second stage

Each feasible first-stage solution feeds the original second-stage graph.
Only an event's other-color **real detector support** creates a virtual check.
Purely boundary-to-boundary projections remain unobserved in stage two. Thus
first-stage parity generates alternatives; it is not imposed again on the final
correction. In particular, spatial corner faults are not pinned off by an added
virtual boundary check.

Temporal measurement edges keep their original detector times in both stages.
Own-color measurement faults remain ordinary two-detector edges in the second
stage. There is no time flattening.

The returned correction minimizes the sum of original event log-odds weights
among feasible candidates. Ties use R/G/B color order, then parity masks
0/1/2/3. A disconnected component with odd syndrome and no implicit sink makes
a candidate infeasible in either stage; it is skipped. Equal-weight bulk-path
alternatives within a sector are not enumerated. Opt-in event compression
preserves input event order for reproducible PyMatching tie choices; default
ordinary matching retains its previous sorted-event order.

## Spatial provenance and circuit-level limits

A single-detector DEM event does not identify whether its missing endpoint is
spatial or temporal. A logical-observable label is never used to guess this.
The low-level `dec::concat::Event::boundary_hint` supports:

- `BoundaryHint::Spatial`: an event has at most one real detector of each color;
  missing colors are supplied by explicit spatial boundary vertices in stage one.
- `BoundaryHint::Temporal`: projected singleton access remains an implicit,
  colorless sink. Same-color measurement pairs remain time edges.
- `BoundaryHint::Unknown`: conservatively retain the original implicit sink
  access instead of inventing a spatial boundary color.

Equivalent events with conflicting hints fall back to `Unknown`. Native
code-capacity and phenomenological generators tag data and readout errors as
`qudec_data` and `qudec_measurement`; Stim preserves these tags in the DEM.
The adapter uses these exact hints, including for corner data faults and
readout singletons. See [ideal noise models](ideal_noise_models.md).

For untagged circuit-level DEMs, the adapter recognizes same-round multicolor
support as spatial and same-color multi-detector support as temporal. Singletons
and cross-round mixed events remain conservative. These unknown/temporal sinks
are not assigned a color and can relax the global parity relationship; the
number of feasible sectors can differ from the ideal six. Consequently the
circuit adapter need not reproduce the candidates of a model with complete
physical spatial-boundary provenance. Existing colorless flag conditioning,
DEM-derived flag-readout probabilities, and flag-context caching are unchanged.

`Decoder::ambiguous_boundary_events()` and
`CircuitResult::ambiguous_boundary_events` count active merged events whose
projected singleton has unknown boundary provenance. `CircuitResult::candidates`
contains all feasible candidates. `Result::boundary_mask` records **first-stage
parity**, not final correction boundary parity: bit 0 fires boundary color
`(color+1)%3`; bit 1 fires `(color+2)%3`. In `restricted_edges`, a boundary vertex
has ID `model.colors.size()+boundary_color`; real detectors keep their IDs.
`CircuitResult::colors` holds the best candidate per color; a color with none
has infinite weight.

## C++ API

```cpp
dec::ConcatMWPM decoder(dem, {.boundary_gap = true});

// Known spatial data event; provenance is independent of its logical label.
dec::concat::Event event{
    p, detector_ids, observable_mask, dec::concat::BoundaryHint::Spatial};

// Existing uniform-restricted-weight parameter, then the boundary ensemble.
dec::concat::Decoder core(model, false, true);
auto candidates = core.decode_candidates(syndrome);
```

Each instance owns mutable PyMatching workspaces. Use separate instances for
concurrent sampler workers.
