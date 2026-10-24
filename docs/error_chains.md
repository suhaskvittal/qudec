<!-- Author: OpenAI Codex. -->
# Connected DEM error-chain testing

Build the `qudec_chains` CMake target. It shares circuit construction, decoder
configuration and decoder-specific options with `qudec`. Its only required
argument is the decoder name. It supports only `cc666` and `cc488`, defaulting
to `cc666` for every decoder.

```sh
cmake --build <build-directory> --target qudec_chains -j 12

# Defaults: distance 5, length 2, SI1000, d rounds, 10,000 random chains.
<build-directory>/qudec_chains concat_mwpm

# Code capacity: K=20 ensemble, 10,000 connected weight-4 patterns at d=9.
<build-directory>/qudec_chains concat_mwpm --code cc666 -d 9 -k 4 \
    --error-model code_capacity -n 10000 --concat-ensemble 20 --seed 42

# Phenomenological noise and the boundary-gap extension.
<build-directory>/qudec_chains concat_mwpm --code cc666 -d 7 -k 3 \
    --error-model phenomenological --concat-boundary-gap -n 10000

# Exhaust all connected weight-2 combinations in a color-code capacity DEM.
<build-directory>/qudec_chains concat_mwpm --code cc666 -d 5 -k 2 -cap --all

# Bell-flagged 4.8.8 circuit noise, with extra failure records.
<build-directory>/qudec_chains concat_mwpm --code cc488 -d 5 -k 2 \
    --error-model si1000 -n 10000 --print-failures 10
```

Options specific to this executable:

| Option | Default | Meaning |
|---|---|---|
| `-d`, `--distance` | 5 | Code distance |
| `-k`, `--chain-length` | 2 | Number of distinct DEM instructions per chain |
| `-n`, `--samples` | 10000 | Number of random connected chains; sampling uses replacement |
| `--all` | off | Enumerate every k-subset, decoding only connected ones |
| `--error-model` | `si1000` | `si1000`, `code_capacity` (also `code-capacity`), or `phenomenological` |
| `--print-failures` | 1 | Maximum failing chains printed; zero suppresses failure records |
| `--print-chains` | 0 | Print the first N chains, including successful decodes |
| `--seed` | 0 | Random chain seed |
| `-pp`, `--print-progress` | off | Print progress every 10,000 decoded chains |

`--code cc666|cc488` selects the color code; the default is `cc666`.
Surface codes are rejected by this executable. SI1000 and phenomenological noise default
to d rounds; `-r`/`--rounds` overrides that. Code capacity always uses one round.
`-cap`/`--code-capacity` and `-pheno`/`--phenomenological` are also accepted.
Opposing noise-model selectors are rejected. These are Z-memory experiments,
matching `qudec`.

The default physical probability is fixed at 0.001. Selection of faults does not
use their probabilities, but decoder weights and flag conditioning can depend
on p. The shared `-p` option remains available for those comparisons.
Phenomenological measurement probability defaults to p and can be changed with
`--measurement-error-rate`. The existing nonuniform-noise options are available.

A chain is a connected set of **distinct positive-probability error instructions
from the original, flattened, undecomposed DEM**. Two events are adjacent if they
flip a common detector, including a flag detector. Separators within one error
instruction do not create independent events. Repeated targets cancel by XOR.
Detectorless events can appear as length-one chains, but cannot join longer
chains. Shared logical observables alone do not connect events. Syndrome and
logical truth are the XOR of all events selected in a chain. Some decoder graphs
require decomposition; they receive their usual decomposed model while chain
selection still uses the original events.

Random mode chooses a seed event uniformly among events in components large
enough for k faults. It repeatedly chooses a uniformly random event from the
current frontier until k distinct events are selected. Separate draws may repeat
a whole chain. This distribution is **not uniform over connected k-subsets**, and
is not probability-weighted IID sampling. It guarantees termination when a
chain exists, including when the requested sample count exceeds the number of
distinct chains.

Exhaustive mode prints the exact integer `C(N,k)` before decoding starts, where
N is the number of positive-probability original DEM events. It scans every
combination once, tests connectivity, and decodes only connected combinations.
The printed total is the number of combinations scanned, not the unknown number
of connected chains. Counts are printed exactly even beyond 64 bits; runtime
counters reject overflow. Enumerating combinations can be very expensive.
`--samples` has no effect in exhaustive mode. No early stop occurs after a failure.

Output includes chains tested, logical failures, their fraction, and the first
requested number of failure records with original zero-based flattened DEM event
IDs, flipped detector IDs, logical truth, and prediction. **The failure fraction
is a stress-test statistic, not an IID logical error rate.** A logical failure is
recorded and testing continues; construction/decoding exceptions terminate with
a nonzero exit status. Empty/impossible chain lengths are rejected. Normal
completion returns zero even when some chains fail. Decoder statistics follow.

All decoder options retain their `qudec` meaning, including `--concat-ensemble`,
`--concat-seed`, `--concat-boundary-gap`, `--correlated`, Libra, Harmony, seeded
matching and union-find settings; see `--help`. Normal Monte Carlo and rare-event
sampler settings do not control this executable. `--method rare_event` is
rejected. Restriction variants retain their existing cc666/two-round limitation
and reject code capacity. Graphlike decoders can reject color-code hyperedges.
Tesseract receives the undecomposed DEM. Each invocation runs on one MPI rank.

## Event supports and ASCII lattice

Each printed chain includes the detector and observable support of every DEM
event separately, even when detectors shared by events cancel in the combined
syndrome. Example notation: `D85=R(13, 4):Z` is global detector D85, red spatial
check 13, round 4, Z basis. `D86=F(13, 4):X` is a **colorless flag** at spatial
check location 13, round 4, X basis. Check IDs follow the generator's face order,
are zero-based, and are shared across rounds and bases. Round numbers are the
zero-based Stim detector time coordinates; the final readout closure is round
`r` after extraction rounds `0..r-1`. Global detector IDs are never substituted
for spatial check IDs.

The ASCII drawing uses the circuit generator's actual data positions and check
supports for either color code. Dots mark data vertices; edges outline the
check polygons. `[R13]` highlights check 13 when its **combined XOR syndrome**
fires in any round or either basis. This is a union across rounds, so two fired
detectors on the same check in different rounds do not cancel in the drawing.
`R13!` marks a fired colorless flag located at that check; `[R13]!` means both
check and flag detections are present. The host check's color does not give the
flag a color. Unbracketed checks are quiet. A legend below the drawing lists
every fired detector with its `(check-id, round)` pair and basis, separating flags
from check detectors. Events that cancel at a detector remain visible in the
individual event listings, but do not highlight the combined-syndrome drawing.

The drawing shows detection events. A DEM event need not identify a unique
physical data-qubit error, so physical fault locations are not inferred.

```sh
# Inspect the first sampled chain even if the decoder succeeds.
<build-directory>/qudec_chains concat_mwpm -d 7 -k 3 -pheno -n 1 --print-chains 1
```

`--debug` enables decoder-specific decisions for concat_mwpm and restriction
variants. Nonempty `dec::Result::debug_info` is printed after the ASCII drawing
for each reported chain. Combine `--debug --print-chains 1` to inspect a successful
chain. See [decoder_debug.md](decoder_debug.md) for matching-path IDs and candidate
selection.
