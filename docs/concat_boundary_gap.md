# ConcatMWPM spatial boundary parity sectors

`qudec concat_mwpm ... --concat-boundary-gap` evaluates four first-stage
boundary parity sectors for each ignored color. It is compatible with the
first-stage graph ensemble and does not change the second-stage matching.

For a restricted lattice that retains boundaries of colors `a` and `b`, each
candidate fixes the parity at neither boundary, only `a`, only `b`, or both.
PyMatching represents those parities by adding the corresponding boundary
vertices to the first-stage syndrome. Spatial data-fault projections include
the missing boundary vertices, so boundary-to-boundary paths carry the
probability of the physical corner events that support them. Temporal events
do not acquire spatial-boundary endpoints.

The decoder runs the second stage for each feasible first-stage result and
returns the lowest model-weight correction across colors, ensemble members,
and parity sectors. Ties retain the first candidate in iteration order. The
ordinary three-color decoder remains the default when the flag is absent.

Author: OpenAI Codex.
