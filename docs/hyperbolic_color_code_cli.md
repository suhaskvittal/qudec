# Hyperbolic color codes in `qudec`

`qudec` accepts a code name from the catalog below. `-d` / `--distance` is
optional and defaults to 5 for planar codes. For a named hyperbolic code, the
name selects its geometry and distance; a supplied `-d` is ignored.

```sh
qudec concat_mwpm --code hycc_888_64_8_6 -p 0.001 --mc-max-samples 1000
qudec concat_mwpm -d 5 --code hycc_688_384_20_12 -cap -p 0.01
```

With `--rounds -1` (the default), circuit-level and phenomenological runs use
the selected code's distance for the number of rounds. Code capacity always
uses one round. An explicit positive `--rounds` overrides the default except
for code capacity, which requires one round. SI1000 uses Bell-flagged extraction;
`-cap` and `-pheno` use ideal product-check measurements. `--superdense` applies
only to planar color codes.

| Code name | Dataset directory | Default rounds |
|---|---|---:|
| `hycc_888_64_8_6` | `8.8.8_L2_D444-32_n64_k8` | 6 |
| `hycc_888_128_12_8` | `8.8.8_L2_D444-40_n128_k12` | 8 |
| `hycc_888_512_36_12` | `8.8.8_L2_D444-58_n512_k36` | 12 |
| `hycc_888_1600_104_14to16` | `8.8.8_L2_D444-109_n1600_k104` | 14 |
| `hycc_888_1344_88_14` | `8.8.8_L2_D444-96_n1344_k88` | 14 |
| `hycc_888_5376_340_18to20` | `8.8.8_L2_D444-271_n5376_k340` | 18 |
| `hycc_888_24288_1522_22` | `8.8.8_L2_D444-PSL2-23_n24288_k1522` | 22 |
| `hycc_888_31200_1954_22to24` | `8.8.8_L2_D444-PSL2-25_n31200_k1954` | 22 |
| `hycc_688_384_20_12` | `6.8.8_L2_D344-18_n384_k20` | 12 |
| `hycc_688_768_36_14` | `6.8.8_L2_D344-22_n768_k36` | 14 |
| `hycc_688_1536_68_16` | `6.8.8_L2_D344-32_n1536_k68` | 16 |
| `hycc_688_3840_164_20` | `6.8.8_L2_D344-45_n3840_k164` | 20 |
| `hycc_688_10368_436_22` | `6.8.8_L2_D344-79_n10368_k436` | 22 |
| `hycc_688_19584_820_24` | `6.8.8_L2_D344-PGL2-17_n19584_k820` | 24 |
| `hycc_688_48576_2028_28` | `6.8.8_L2_D344-PGL2-23_n48576_k2028` | 28 |

Three dataset entries report a distance interval instead of an exact distance.
Their code names show that interval; default rounds use its proven lower bound.
The catalog is in `main/code_catalog.h`. It contains the CLI name, data
directory, and default round count for each non-planar code.
