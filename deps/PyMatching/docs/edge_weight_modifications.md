<!-- Author: OpenAI Codex. -->
# Scoped C++ edge-weight modifications

`pm::EdgeWeightModifications` provides temporary absolute weights on an existing
`pm::Mwpm` object, without rebuilding its topology or scanning all graph edges:

```cpp
pm::EdgeWeightModifications bp_weights(mwpm);
const auto edge = bp_weights.edge_id(u, v);  // v == -1 means boundary
bp_weights.modify_edge_weight(edge, integer_weight);
// Run one or many completed decoding calls with these weights.
{
    pm::EdgeWeightModifications ensemble_weights(mwpm);
    ensemble_weights.modify_edge_weight(u, v, perturbed_integer_weight);
    // Decode; matching growth and path extraction see the same edge weight.
}  // Restore the enclosing BP weights.
bp_weights.clear_edge_weight_modifications();  // Restore original weights.
```

A handle represents one undirected edge. Registration canonicalizes endpoint
order and caches pointers to both adjacency directions in both the matching and
path-search graphs. A boundary has one direction in each graph. Registration does
local neighbor lookup; subsequent updates by ID are constant time. Each scope
saves the original values only on its first modification of that edge. Clearing
is proportional to the number of touched edges, retains registered handles, and
is also performed by the destructor, including exception exits. Repeated sets
replace weights rather than accumulate deltas; both increases and decreases are
supported. Nested scopes must be destroyed/cleared in reverse order.

Weights use **existing graph integer units**, not probabilities. They must be
even, nonnegative, and at most `INT32_MAX-1`. Even weights preserve integer
collision times. The existing graph normalization scale remains fixed; callers
converting likelihoods must use that scale and choose adequate headroom. This API
does not recompute a global maximum/scale and cannot exactly reproduce an eager
update that chose a different scale for each perturbation. It does not alter
endpoints, availability, observable labels, or correlation rules. Numeric zero
means zero-cost, not disabled. Public overrides require a nonnegative graph and
an existing path-search graph.

The owning matcher must outlive all scopes. Graph topology, adjacency storage,
and the matcher location must remain unchanged while handles exist. Weights may
change only between completed solves, never in a callback while search events
are active. No callback branch or table lookup is added to edge traversal; the
API writes the already allocated weights directly. An unmodified edge has no
per-shot cost. This is a sparse eager update of selected edges, rather than lazy
on-traversal random-number generation.

## Correlated matching

The existing correlated decode entry points still perform two solves. They first
decode to physical edges, then apply the precomputed implied-weight rules through
this generic API, and decode again. The correlation rule retains the previous
**lower-only minimum** behavior even though the general API permits increases.
An internal scope restores only that correlated pass's changes, preserving any
outer BP or ensemble overrides. Conditional probability rules themselves remain
those precomputed when loading the graph; overriding base weights does not infer
new correlation probabilities.

Both observable decoding APIs and `decode_detection_events_to_edges_with_edge_correlations`
use the same scope. Original graph-specific `reweight_for_edges`/`undo_reweights`
remain available for compatibility and baseline tests. Signed preprocessing is
unchanged for ordinary matching. The DEM correlation loader already rejects
probabilities greater than one half; legacy signed-graph correlation bookkeeping
is retained as a fallback for manually constructed signed matchers.
