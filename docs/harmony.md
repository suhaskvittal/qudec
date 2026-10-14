# Harmony decoder

Contributor: OpenAI Codex.

`qudec harmony` runs independently perturbed two-pass correlated MWPM
members and pools their logical predictions by voting. The implementation
follows the vote-pooling variant of arXiv:2401.12434; it does not implement
the likelihood pooling favored in the paper's surface-code benchmarks.
SI1000 circuits include opposing-basis detectors automatically.

```sh
mpirun -np 64 ./build/qudec harmony 9 -r 9 -p 0.001 \
    --harmony-ensemble-size 20 --harmony-full-ensemble \
    --mc-max-samples 1000000 --mc-stop-at-errors 1000000 --seed 0
```

Options:

- `--harmony-ensemble-size`: total members, default 100, must be positive.
- `--harmony-sentinel-size`: initial members in layered mode. Default 0
  means the smaller of four and the total ensemble size.
- `--harmony-full-ensemble`: run all members on every shot. By default,
  unanimous sentinels stop decoding early; disagreement runs the remaining
  members. Use full mode to measure latency versus total ensemble size.
- `--harmony-seed`: ensemble perturbation seed, separate from the sampling
  `--seed`. Members are deterministic and identical across MPI ranks;
  sampling streams differ across ranks.

Output includes logical error rate, expanded shots, total and mean members
run, and mean/minimum/maximum decoding latency in microseconds. Timings
cover a whole syndrome history, including empty syndromes, and exclude
ensemble construction, sampling, and MPI reductions. All ensemble members
are constructed before sampling. Each rank owns its own ensemble.

For example, compare full ensembles on the same sampling streams:

```sh
for size in 1 4 8 16 32 64 100; do
    mpirun -np 64 ./build/qudec harmony 9 --harmony-ensemble-size "$size" \
        --harmony-full-ensemble --mc-max-samples 100000 \
        --mc-stop-at-errors 100000 --seed 0
done
```

The decoder currently supports one logical observable. Perturbed
probabilities are clamped below one half to fit the matching backend's
nonnegative conditional-weight representation. Vote ties select frame zero.
These details and vote pooling should be considered when comparing against
the paper's reported performance.
