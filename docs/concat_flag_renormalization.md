<!-- Author: OpenAI Codex. -->
# Flag renormalization for concatenated MWPM

This opt-in implements a DEM-based adaptation of Section 3.4 (Section III D in
[the arXiv version](https://arxiv.org/pdf/1911.00355)) of Chamberland et al.,
*Triangular color codes on trivalent graphs with flag qubits*,
[New Journal of Physics 22, 023019](https://doi.org/10.1088/1367-2630/ab68fd).

```sh
qudec concat_mwpm -d 5 -p 0.0005 --code cc666 --concat-flag-renormalization
qudec_chains concat_mwpm -d 5 -r 5 -k 2 -p 0.0005 --all \
    --concat-flag-renormalization --debug --print-failures 1
```

The corresponding `dec::ConcatMWPM::config_type` fields are
`flag_renormalization=true` and `physical_error_rate=p`. The supplied p is the
nominal physical error rate, rather than a DEM-event probability or flag
marginal. Default decoding retains its existing policy.

## Probability rule

For the observed relevant flag syndrome, count m distinct flagged extraction
circuits. Flags sharing `(face_x, face_y, extraction_round, extraction_basis)`
count once. With the current Bell circuits there is one flag per extraction.
The count spans the decoding window. As in the existing decoder, Z memory uses
X-extraction flags and X memory uses Z-extraction flags.

The adapter groups faults by their complete memory-check support, including
time, and chooses a representative using the existing minimum symmetric flag
difference and probability tie-break. In this mode, every required flag of an
eligible flagged representative must have fired. A class without an eligible
flagged member can use its unflagged representative; otherwise it is omitted.

For a selected supported event, use its raw DEM probability. For an ordinary
unflagged representative, multiply its probability by `p^m`. Unsupported
flag-only classes are omitted, corresponding to infinite matching weight.
Matching and final event scores use `-log(probability)`.

First-stage projection can combine several different full-support events into
one edge. If any are supported by the flags, that projected edge belongs to
`C_f`: its probability is computed from the supported contributions only.
Otherwise its contributors already carry the `p^m` factor. Independent DEM
contributions retain the project's XOR-probability merge; this agrees with the
paper's sum to leading order. Second-stage edges retain their full-event
provenance and the same probability policy.

The graph's `C_f` includes ordinary projected edges from flagged single faults,
as well as same-color hook edges. It is constructed from DEM support, rather
than solely from detector color or spatial adjacency. All matching operations
still use PyMatching.

## Scope

This transfers the probability rule to the existing concatenated decoder and
Bell-flagged circuits. The source paper analyzes different extraction circuits,
so its distance claims do not establish a distance guarantee for this
implementation. The implementation leaves the existing matching stages and
correction selection algorithm in place.

The scale `p^m` is global across the observed relevant extraction circuits in the
window. A local or per-round alternative would be a separate adaptation.
An observed flag with p=0, or numerical underflow of a renormalized probability,
is reported as an error.

## Sparse matching-weight updates

Contributor: OpenAI Codex.

Flag-renormalized decoding automatically applies numerical weights through
PyMatching's `EdgeWeightModifications` API in both matching stages and every
first-stage ensemble member. The existing flag selection, p^m probabilities,
projection rules, fixed seeded perturbations, correction scoring, and integer
normalization remain unchanged. Handles cache both adjacency directions in the
matching and path-search graphs; stable-topology updates touch only changed
integer weights. A changed normalization still requires global requantization.

Changes to enabled edges or compact virtual endpoints use the existing topology
rebuild after clearing handles and retain the legacy backend for that graph
thereafter. This avoids repeated handle registration on volatile hook supports.
This optimization does not change other ConcatMWPM variants or add a CLI flag.
