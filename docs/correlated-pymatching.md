# Correlated PyMatching

Contributor: OpenAI Codex.

Use the existing PyMatching decoder with `--correlated`:

```sh
mpirun -np 64 ./build/qudec pymatching 9 --correlated \
    --mc-max-samples 1000000 --mc-stop-at-errors 1000000
```

The flag enables PyMatching's existing two-pass correlated MWPM routine.
The first matching selects edges, its correlation rules reweight related
edges, and the second matching supplies the logical prediction. The wrapper
uses the backend directly and does not add ensembling. Without the flag,
PyMatching remains uncorrelated.

SI1000 circuits automatically include opposing-basis detectors in this mode.
For an uncorrelated comparison on the same circuit, use
`--include-opposite-basis-detectors`. Both Monte Carlo and rare-event sampling
are supported. The existing PyMatching latency statistics include both passes
and indicate `correlated=true`.

C++ callers can use `dec::PyMatching(dem, {.correlated=true})`. The DEM must
contain decomposed components for the backend's correlation rules. Correlated
mode currently rejects complementary-gap estimation and matched-endpoint
reporting; those existing paths are implemented for uncorrelated matching.
