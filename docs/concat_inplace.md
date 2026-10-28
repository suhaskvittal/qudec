<!-- Author: OpenAI Codex. -->
# Persistent ConcatMWPM graph updates

ConcatMWPM now updates its existing PyMatching graphs when the relevant flag
syndrome changes. This applies to both the conditional flag policy and
`--concat-flag-renormalization`, including ensembles.
No additional CLI flag is required.

Each decoder instance compiles the union of its possible check-support classes
and projected edges once. It retains `3*K` first-stage matching graphs and three
second-stage graphs. A flag change selects the class representatives, updates
their probabilities and logical frames, aggregates projected probabilities,
and writes the resulting integer weights to both PyMatching backend graphs.

Disabled events are physically absent from the active adjacency lists. They
are not represented by large finite weights. Matching graph objects and their
detector-node arrays remain allocated; adjacency arrays reuse their capacities.
Connected-component feasibility information is recomputed only when active
endpoints change. Normalization is recalculated using the existing PyMatching
conversion rule, and ephemeral node/queue state is reset before updates.

Projection membership is precomputed. Active virtual detector IDs stay compact
and retain the original ordering. Ensemble perturbation multipliers are reused
while the number of active projected edges stays unchanged; when it changes,
the original seeded normal-draw sequence is regenerated. This preserves the
existing ensemble prefixes and probability clamps.

The decoder stores only its current flag state. It no longer caches full graph
sets for whole flag patterns. Repeated consecutive flag states reuse their
current weights. `flag_cache_capacity` and the corresponding positional C++
constructor argument remain accepted for source compatibility; they no longer
control graph storage. Instances remain mutable and require one instance per
sampler worker.

Diagnostic incidence neighborhoods are populated only through `inspect()` or
debug decoding. Ordinary decoding does not compute them.

`concat_inplace_test` compares persistent updates against freshly constructed
graphs, including disabled edges, logical-frame changes, unit probabilities,
boundary sectors, and ensembles. It checks exact paths and virtual syndromes as
well as weights and logical predictions. End-to-end reference comparisons and
benchmarks are recorded in `.agentic_artifacts/concat_inplace/`.
