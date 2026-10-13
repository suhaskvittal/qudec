# PDZ + Libra syndrome recording

Contributor: OpenAI Codex.

Use `qudec pdz_libra 9 --record-syndromes`. The default destination is a unique
run directory under `out/`. Override it with `--record-syndromes-dir DIR`.
Each MPI rank writes `rank_N.bin.gz`; existing record files are rejected.
`metadata.json` stores the simulation configuration. Recording includes every
decoded shot, including empty syndromes, without changing decoding decisions.

`--pdz-removal-threshold` defaults to `0.5`: if PDZ removes at least 50% of
the input's fired syndrome bits (both bases), return the original uncorrelated
PyMatching result and skip Libra completely. The PDZ correction is not added
to this returned result. Use `1.0` to disable bypass for nonempty residuals.
Thresholds range from 0 to 1 inclusive. Empty post-PDZ syndromes always bypass
Libra, including shots that started empty: there is nothing left to decode.
Bypass counts are printed in decoder statistics, and Libra timings cover only
shots on which Libra was called.

`libra_called_shots` explicitly counts those calls. The weighted decoding
latency is `(libra_called_shots / sampled_shots) * mean_Libra_latency` in
microseconds per full history, with zero cost assigned to bypassed shots.
This metric includes only Libra's decode time; the initial PyMatching call,
PDZ processing, and recording are outside it. MPI runs compute it from global
shot counts and global timing totals, rather than averaging worker averages.

## Run the complete experiment

```sh
python3 scripts/run_pdz_libra_experiment.py --distance 9 --shots 1000000
```

The runner discovers an existing MPI build, uses 64 ranks by default, records
exactly the requested shot count, then runs the analysis. It creates a unique
subfolder of `out/` containing rank files, configuration, commands, logs, CSV
tables, and reports. It prints a summary with ensemble calls, errors fixed,
and errors introduced by PDZ removal fraction. It does not compile the decoder.

Use `--qudec PATH` to select an executable, `--cores N` to change rank count,
or `--cores 1` for serial execution. `--p` defaults to 0.001 and `--rounds`
defaults to the distance. PDZ and ensemble settings are also exposed; see
`--help`. `--dry-run` checks the configuration and prints commands without
creating a run directory or sampling. The shot budget is global across ranks;
the error-count stop is disabled for this experiment.

Analyze and combine ranks with:

```sh
python3 scripts/analyze_pdz_libra_records.py out/RUN_DIRECTORY
```

The script requires NumPy. Reports and CSV tables are written inside that run
directory. Monte Carlo records represent the sampled distribution. Rare-event
records describe conditional samples and are **not** physical-frequency weighted.

## Binary format, versions 1 and 2

All integers and the IEEE binary64 gap are little endian. The gzip stream starts
with a 40-byte header, followed by fixed-size records until EOF. Bit `d % 8` of
byte `d // 8` represents detector `d`; final-byte padding bits are zero.

| Header offset | Type | Meaning |
|---|---|---|
| 0 | 8 bytes | Magic `PDZLIBR1` |
| 8 | uint32 | Version, 2 for new files; reader also supports 1 |
| 12 | uint32 | Record size, `24 + 2 * ceil(detectors / 8)` |
| 16 | uint64 | Detector count |
| 24 | uint64 | Observable count, 1 for Libra |
| 32 | uint32 | MPI rank |
| 36 | uint32 | Header size, 40 |

| Record offset | Type | Meaning |
|---|---|---|
| 0 | uint32 | Initial syndrome Hamming weight |
| 4 | uint32 | Post-PDZ syndrome Hamming weight, both bases |
| 8 | float64 | Libra's initial unsigned complementary weight gap in dB |
| 16 | uint8 | Ensemble used, 0 or 1 |
| 17 | uint8 | Ensemble used **and** final pipeline frame differs from initial PyMatching |
| 18 | uint8 | Initial uncorrelated PyMatching logical frame |
| 19 | uint8 | Final PDZ + Libra logical frame |
| 20 | uint8 | True logical frame from the sampler |
| 21 | uint8 | Version 2: Libra bypassed, 0 or 1; version 1: reserved zero |
| 22 | 2 bytes | Reserved, zero |
| 24 | packed bits | Original syndrome |
| 24 + ceil(detectors / 8) | packed bits | Post-PDZ syndrome |

The gap is computed on the **post-PDZ** syndrome before ensembling; it is not
the final synthesized gap. In version 1, empty residual syndromes have a
placeholder gap of zero; in version 2 they bypass Libra and have a NaN gap.
They never enter ensembling. The frame-change flag compares against the
initial PyMatching result, rather than a hypothetical Libra decode without
ensembling. A frame change can help or harm; ground truth distinguishes these.

Bypassed shots have a NaN gap, ensemble/change flags zero, and final frame equal
to initial PyMatching. No confidence gap is computed for them. The reader
validates these fields against the configured removal threshold, excludes
bypassed shots from gap means and ensemble-association fits, and includes their
counts in the per-bin CSV tables. Older version-1 files remain readable.

PDZ's removed-bit count is `initial_hw - post_hw`. Dividing it by two matches the
histogram convention; accepted boundary matches remove just one bit. Analysis
also considers the fraction of initial bits removed and compares samples within
initial-Hamming-weight strata. Observed associations do not establish a safe
rule for skipping complementary-gap estimation.
