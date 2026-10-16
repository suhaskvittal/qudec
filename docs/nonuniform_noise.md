# Nonuniform component error rates

All SI1000 circuit generators support a fixed, random distribution of component
error rates. This represents device disorder: a component's probability stays
fixed for all circuit rounds and simulation shots. It is sampled when building
the circuit, separately from sampling faults during shots. The sampler and
decoder both use the resulting nonuniform detector error model, so decoder
weights know the actual device rates. Monte Carlo and rare-event sampling already
support unequal DEM channel probabilities.

`ComponentErrorDistribution` in `src/component_error_distribution.h` provides
the common representation. For a component with nominal probability `q`, it uses

```text
Z ~ Normal(0, 1)
q_device = min(channel_cap, q * exp(sigma * Z - sigma^2 / 2))
```

`sigma` is the standard deviation of the **natural logarithm** of the multiplier.
Before capping, the arithmetic mean remains `q`, the median is
`q * exp(-sigma^2/2)`, and the coefficient of variation is
`sqrt(exp(sigma^2)-1)`. A lognormal distribution has unbounded support, so rates
must be capped. `channel_cap` is 1 for scalar Pauli errors, 3/4 for
`DEPOLARIZE1`, and 15/16 for `DEPOLARIZE2`. The depolarizing caps give maximally
mixed channels: Stim cannot convert larger, over-mixing rates into a DEM. The
result is lognormal below the cap with a point mass at the cap, rather than a
pure unbounded lognormal. Saturation reduces the actual arithmetic mean below
`q`, especially at large spread. SI1000 prefactors remain the
existing ones; the distribution multiplies each resulting channel probability.
Zero spread preserves the original circuit exactly, and zero physical noise
remains exactly zero for every spread. Very small probabilities may underflow to
zero in floating point.

A component is identified by the Stim noise channel type, its ordered target
qubit(s), and its instruction tag. One-qubit channels are split per qubit;
`DEPOLARIZE2` is split per ordered qubit pair. The same identity shares its
multiplier across rounds, circuit blocks, and different nominal probabilities.
For example, `DEPOLARIZE1(p/10)` and `DEPOLARIZE1(p)` on the same qubit share a
multiplier. Untagged reset and measurement `X_ERROR` channels on that qubit also
share a multiplier; distinct tags can distinguish component mechanisms where
needed. Different channel types or target pairs have independent draws. This is
a channel/component model, not independent disorder for each space-time gate
occurrence, and it does not impose one common multiplier on every mechanism of
a physical qubit.

The device seed fixes these multipliers independently of instruction ordering,
target batching, repeat-block expansion, sampling seed, and MPI rank. Changing
the nominal physical error rate keeps the same device multipliers, allowing
comparisons across physical error rates. Deterministic SplitMix64 hashing and
Box-Muller sampling avoid standard-library random-distribution differences;
transcendental-library floating-point rounding can still differ by platform.
The circuit transformation supports `X_ERROR`, `Y_ERROR`, `Z_ERROR`,
`DEPOLARIZE1`, and `DEPOLARIZE2`, covering every current SI1000 generator. It
preserves annotations, flags, detector coordinates, observables, and repeat
blocks. At nonzero spread, unsupported noise channels raise an error.

## API and CLI

The optional final argument on `sc_si1000`, `toric_si1000`, `cc666_si1000`, and
`cc488_si1000` is a `ComponentErrorDistribution`. Existing direct calls remain
valid through the default zero-spread argument.

```cpp
ComponentErrorDistribution device{.log_stddev = 0.7, .seed = 123};
auto surface = sc_si1000(5, 5, 0.001, false, true, device);
auto color = cc666_si1000(5, 5, 0.001, false, device);
// Or transform an already-generated circuit:
auto nonuniform = apply_component_error_distribution(uniform_circuit, device);
```

```bash
./build/qudec pymatching 5 -p 0.001 --noise-log-stddev 0.7 --noise-seed 123 \
    --seed 456 --mc-max-samples 100000 --mc-stop-at-errors 100001
```

`--noise-log-stddev` defaults to zero. `--noise-seed` defaults to zero and is a
nonnegative integer; it is independent of the syndrome-sampling `--seed`.
Use the same device seed on all MPI ranks; the CLI does this automatically.

## Five-device evaluation

`scripts/eval_decoder_nonuniform.py` requires only Python's standard library and
a compiled `qudec`. It sweeps distances 3, 5, 7, 9, 11, and 13 by default, using
the same shot schedule as `eval_decoder_uniform.py`: one million shots for
each of distances 3 through 9, ten million at distance 11, and one hundred
million at distance 13. Each distance evaluates exactly five random models.
`--distances` overrides the distance list; `--shots` overrides the number of
shots per model at every distance. Custom distances absent from the default
schedule use one million shots unless `--shots` is supplied. An optional
positional distance preserves the previous single-distance invocation; combining
it with `--distances` is rejected.

The script draws five distinct device seeds from a master seed and reuses them
at every distance, defining five reproducible device families. A shared seed
fixes the multiplier for a shared component identity; lattice qubit numbering
can change between distances, so this is not a claim that smaller devices are
physical subpatches of larger devices. Separate syndrome-sampling streams are
defined per distance and model, and remain stable if the sweep is reordered or
restricted. This fixed-shot design supports a binomial Wilson 95%
confidence interval, including a positive upper bound when no errors occur.
The CLI's error-count stopping condition is set above the total requested shots.

```bash
# Full distance sweep with the default shot schedule:
python3 scripts/eval_decoder_nonuniform.py pymatching -p 0.001 \
    --noise-log-stddev 0.7 --seed 42 --qudec build/qudec

# Smaller custom sweep, MPI, and extra decoder arguments:
python3 scripts/eval_decoder_nonuniform.py libra --distances 3 5 --cores 4 --shots 100000 \
    --qudec build/qudec -- --libra-ensemble-size 20

# Existing single-distance syntax remains supported:
python3 scripts/eval_decoder_nonuniform.py pymatching 5 --shots 100000
```

The optional `--rounds` defaults to distance; `--cores` defaults to one. Extra
decoder arguments after `--` cannot override experiment, noise, sampling, or
output-count settings. The default output directory follows the uniform
evaluator's naming convention, adding the spread and `_nonuni` suffix:
`out/<decoder>_sc_error_rates_p<p:.2e>_std<stddev:g>_nonuni`. For example, the
first command writes to `out/pymatching_sc_error_rates_p1.00e-03_std0.7_nonuni`.
`--output PATH` uses exactly that directory instead. Each distance has a folder
such as `d3/`, with five model folders `model_0/` through `model_4/`. Each model
folder contains a readable `setup.log`, a structured `setup.json` with all
parameters and the command, and the complete decoder execution `output.log`.
The global `results.csv` contains seeds, experiment
parameters, exact sample/error counts, LER, and confidence limits, and is updated
after each completed model across all distances. A five-device mean is printed
for each distance. These means are descriptive; the
per-device intervals describe shot noise, not uncertainty over the distribution
of devices. Reusing an output directory overwrites that evaluation's files.

The regression test `component_error_distribution_test` checks zero-spread
identity, seeded reproducibility, fixed rates across rounds and expanded repeat
blocks, preserved multipliers across physical error rates, empirical lognormal
moments, probability bounds, input validation, and all four SI1000 generators.
