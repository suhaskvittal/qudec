<!-- Author: OpenAI Codex. -->

# Superdense color-code circuits

Implemented from [Denser Planar Color Codes](https://arxiv.org/html/2609.21376v3),
Section 2.2, Appendix B and Appendix E, and the
[published construction sources](https://doi.org/10.5281/zenodo.22820614).
Source attribution and licensing are in [superdense_source_NOTICE.md](superdense_source_NOTICE.md).

Use `--superdense` with `qudec` or `qudec_chains`, selecting `--code cc666` or
`--code cc488`. Existing Bell-flagged extraction remains the default.

```sh
qudec tesseract -d 7 --code cc666 --superdense -p 0.0005
qudec tesseract -d 7 --code cc488 --superdense --superdense-start-phase 1 -p 0.0005
qudec_chains tesseract -d 7 --code cc488 --superdense -k 2 -n 1000
```

CNOT is the default in the C++ API and the only representation exposed by the
CLI. Native CZ compilation remains available through the C++ API for circuit
construction and reference validation. The two representations have different
physical fault locations and therefore different circuit-level error models.
The selector requires SI1000; the existing capacity/phenomenological generators
remain ideal-check models.

Use Tesseract for circuit-level decoding with this implementation. The current
ConcatMWPM wrapper rejects some superdense fault supports with its existing
"Circuit event does not fit Algorithm 1 for every color" check. This change
does not extend that decoder's supported hyperedge shapes.

The C++ entry points in `circuit_generator.h` are `cc666_superdense` and
`cc488_superdense`. They take a `CircuitNoisePolicy`, memory basis, optional
opposite-basis detectors and component distribution; 4.8.8 additionally accepts
starting phase 0/1. Native compilation can be enabled with the final boolean.
There is no Python dependency at runtime or build time.

## Geometry and extraction

6.6.6 accepts every odd distance at least 3. It has `(3*d*d+1)/4` data
qubits, two ancillas per face, and eight CNOT layers per complete extraction
cycle. Both X and Z checks are obtained in that one cycle, without Bell flags.
The planar contact order is the standard Chromobius superdense construction.

4.8.8 implements the paper's brickwork triangular family at `d=4*m-1`:
3, 7, 11, 15, ... . Other distances are rejected. It has
`(d*d+2*d-1)/2` data qubits and `(d*d+3*d-4)/2` ancillas. Red bulk checks
are squares; green/blue bulk checks are octagons, with truncated boundary
checks. Ancillas are shared between faces. The complete 90-class radius-five
local rule includes the boundaries and corners, with twelve CNOT layers.
Successive cycles alternate forward/reversed layer order without reversing
CNOT direction. Both starting phases are available.

The native 4.8.8 compilation includes the published boundary Bell-polarity
change, alternating data Hadamard frames, local H cancellation and earliest
per-wire scheduling, followed by simultaneous ancilla measurement. Recurring
cycles occupy 23 physical moments at d=3 and 24 at larger supported sizes.
SI1000 noise is injected by the existing policy: p after CZ/CX, p/10 after
one-qubit gates and on idle qubits, 2p reset errors, 5p readout errors, and
additional 2p depolarization on every noncollapsed qubit during a reset/readout
moment, including qubits with simultaneous gates. Component-rate disorder is
applied using the existing shared distribution mechanism.

## Detector and logical parities

The generator derives each phase's incoming and outgoing face-sign parities,
ideal record constraints and logical-sign increments by exact binary CSS Choi
elimination. Algebraic input references do not become physical circuit qubits.
No feedback gates are needed: detector/observable annotations track the frame.
Interior detectors compare the current incoming face sign with the previous
outgoing face sign. First-round detectors use the prepared memory basis;
closure detectors compare the last outgoing sign with final data readout.
Logical annotations include the per-cycle sign increments and final all-data
logical parity. Opposite-basis interior detectors remain optional.

6.6.6 uses repeated complete forward extraction cycles, followed by direct data
readout, with this same instrument-derived closure. The bundle's pinned
Chromobius memory instead ends with an inverted initialization chunk. Thus the
6.6.6 memory boundary circuit is not a byte-for-byte reproduction of that file.

Detector metadata is `(x,y,round,color,basis,is_flag)`. Check colors are
0/1/2 = red/green/blue, and check bases are 0/1 = X/Z.
4.8.8 has `(d-1)/2` independent, colorless, single-cycle routing flags per
round. Their color coordinate is -1 and `is_flag=1`; they are never compared
with flags in previous rounds. These are ideal-zero record constraints, not
extra face checks. A flag constraint in the X (Z) stabilizer sector detects Z
(X) faults, so its basis metadata is Z (X), allowing the decoder's existing
opposite-extraction-basis filtering to retain it for X (Z) memory. These
routing flags do not belong to an X/Z Bell-check extraction block; this
metadata convention is specific to the superdense routing constraint.

Chain output draws the selected extraction geometry, labels checks by spatial
ID and round, and lists independent routing flags with spatial coordinates and
round. Chain stress tests keep the existing flags-OR/checks-XOR convention.

## Validation

`superdense_color_code_circuits_test` covers supported distances through 15,
both memory bases, CNOT/native gates, both 4.8.8 starting phases, one/two/d
rounds and both opposite-basis settings. It verifies physical-layer conflicts,
resource/detector counts, color labels, absolute noiseless parities, noiseless
sampling, deterministic detectors/logicals, noisy DEM generation and absence
of undetected single-fault logical events. Geometry construction also verifies
stabilizer rank, commutation and three-coloring.

For an independent gate/timing comparison against a downloaded release, run
the test executable with the release directory as its sole argument. It checks
all 32 published 4.8.8 ideal physical schedules at d=3,7,11,15, both bases,
phases and gate sets. These tests do not prove arbitrary-distance circuit
distance or a particular decoder's error-correction performance.
