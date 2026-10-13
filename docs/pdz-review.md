# PDZ implementation and review notes

Contributor: OpenAI Codex.

## Implemented

The DEM and GZIP constructors create an `AffinityTable` and delegate to the main
constructor taking the table by value. Convenience overloads preserve calls that
omit configuration, avoiding GCC's restriction on default arguments using nested
configuration types with member initializers.

`apply()` removes accepted assignment endpoints from the supplied syndrome and
returns the accepted matching data and its observable parity. XOR this parity with
the subsequent decoder's prediction for the residual syndrome. Matching input must
cover every fired detector exactly once; any number of pairs may use boundary ID
`-1`. Calls validate the matching before modifying the syndrome. PDZ instances use
mutable scratch space and must not be shared between concurrent calls.

Basis classification uses the SI1000 generator's checkerboard coordinates and
first detector as the memory-basis reference. Without those coordinates, graphlike
DEM components identify logical-bearing checks. Observable-free DEMs without
coordinates treat all checks as one basis. Models with neither suitable coordinates
nor graphlike decompositions are rejected.

Opposing-basis scoring (`ignore_opposing_basis_detectors=false`) remains deferred
at the user's request. The constructor rejects this setting explicitly.

## Straightforward fixes

- Corrected undefined identifiers, syntax errors and missing includes.
- Stopped Hypergraph template macros from leaking into later headers.
- Cleared all scratch structures between calls, including after failed validation.
- Corrected the reversed event-list classification and opposing-basis condition.
- Excluded a matched pair's own endpoints from its competing events.
- Compared each pair to other events in that pair's basis.
- Kept the shared boundary out of maps keyed by unique detector assignments.
- Accumulated both matching data and returned observable correction.
- Rejected negative score margins instead of terminating the process.
- Fixed flat distance indexing, boundary storage, table symmetry, packed entry
  updates, quantization and saturation of long/disconnected paths.
- Checked GZIP reads/writes and used a fixed-width four-byte detector count.
- Added the preprocessor sources and ZLIB linkage to the library build.

## Decisions needing review

1. **No competitors accepts a pair regardless of its distance.** The user changed
   both competition helpers to `min` and initialized their scores to the maximum
   `value_type` (127). This fixes nearest-competitor selection and the previous
   rejection of isolated pairs. With no competitor, even a long-distance or
   disconnected pair passes the default threshold. This is the current policy;
   whether to add an absolute match-distance limit remains a separate decision.
2. **Distance omits error probabilities.** The user selected five-bit values
   (32 buckets). Values 0..30 are exact hop counts; 31 saturates longer and
   disconnected paths. The default threshold of two means two hops. Error
   probabilities do not participate in table construction, unlike
   `measure_affinity()`'s probability-based paths. Storage currently uses one
   byte per five-bit value, leaving three unused bits.
3. **Correlated faults can create cross-basis shortcuts.** The hypergraph includes
   the complete support of a fault; the distance construction joins every pair in
   that support. Paths may traverse the opposing basis even when PDZ ignores
   opposing-basis fired events. Whether the table should instead use separated
   graphlike components needs an algorithmic decision.
4. **Table construction remains cubic.** Floyd-Warshall uses O(N^3) work and an
   O(N^2) integer workspace. This becomes costly for large circuit histories.
   Repeated bounded breadth-first searches would be a possible future optimization.
5. **Basis inference is specific to supported DEMs.** Checkerboard coordinates
   assume the supplied SI1000 layout. Graph connectivity without coordinates can
   misclassify disconnected memory-basis components with no logical-bearing fault.
   General layouts would benefit from an explicit detector-basis map.
6. **Old table files must be regenerated.** The repaired packed layout includes
   the boundary, now uses five-bit values, and has a different payload size. The file has
   no version identifier or checksum binding it to a DEM; detector-count equality
   alone cannot detect a table from a different DEM with the same detector count.

## Verification

`tests/pdz_test.cpp` checks packed distances and symmetry, boundary distances,
repeated application, nearby competitors despite distant events, acceptance of
isolated pairs, negative-margin rejection, observable correction, shared boundary assignments, invalid
matching recovery, GZIP round trips and constructor delegation, and application to
X- and Z-memory SI1000 circuits with opposing-basis detectors enabled. These are
functional checks, not a logical error rate validation of the scoring rule.
