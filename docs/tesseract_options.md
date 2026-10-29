<!-- Author: OpenAI Codex. -->

# Tesseract options in qudec

Both `qudec tesseract` and `qudec_chains tesseract` accept the same search
options. Existing circuit, noise, sampling, MPI and component-rate options
remain available. For example:

```sh
qudec tesseract -d 7 --code cc666 --superdense -p .0005 \
    --tesseract-beam 20 --tesseract-beam-climbing \
    --tesseract-pqlimit 1000000 --tesseract-num-det-orders 21 \
    --tesseract-det-order-bfs --tesseract-det-order-index \
    --tesseract-det-order-coordinate --tesseract-det-order-seed 518278944
```

Defaults preserve the previous `TesseractConfig{.dem=dem}` behavior.

| Option | Default | Effect |
|---|---:|---|
| `--tesseract-beam` | 5 | Maximum detector beam; 0–65535, or -1 for the upstream unrestricted sentinel (65535). |
| `--tesseract-beam-climbing` | Off | Cycle beam values from zero through the specified maximum across ordering trials. |
| `--tesseract-pqlimit` | 200000 | Positive priority-queue push limit; -1 means unlimited. |
| `--tesseract-num-det-orders` | 1 | Total number of detector orders. |
| `--tesseract-det-order-seed` | 0 | Fixed ordering seed, separate from the sampling `--seed`. |
| `--tesseract-det-order-bfs` | Off | Include BFS orders. |
| `--tesseract-det-order-index` | Implicit | Include index orders; used automatically when no method flags are selected. |
| `--tesseract-det-order-coordinate` | Off | Include coordinate orders. |
| `--tesseract-allow-revisit-dets` | Off | Disable the default prohibition on revisiting detectors. |
| `--tesseract-no-merge-errors` | Off | Disable merging identical error mechanisms. |
| `--tesseract-det-penalty` | 0 | Finite, nonnegative extra cost per residual detector. |
| `--tesseract-verbose` | Off | Print per-search diagnostics and final selected hyperedges, including original flattened DEM error IDs, probabilities, costs, detectors and observables. |
| `--tesseract-sparsify-errors` | Off | Enable sparsification; requires a positive base degree. |
| `--tesseract-sparsify-base-degree` | -1 | Retain all errors with detector degree at most this value. |
| `--tesseract-sparsify-max-degree` | -1 | Maximum degree of optional errors; -1 means no bound. |
| `--tesseract-sparsify-reactivate-limit` | -1 | Optional errors to reactivate; -1 uses Tesseract's automatic choice. |

When multiple ordering methods are selected, the requested **total** is
distributed round-robin in BFS, index, coordinate order, skipping unselected
methods. The count must cover every selected method. Order j uses seed+j.
Thus 21 orders with all three flags gives seven of each type. Ordering seeds
are shared across MPI ranks; the sampler keeps its separate rank/group streams.

With beam climbing, Tesseract runs `max(beam+1, number_of_orders)` trials,
cycling beam values and order indices. Without climbing, every order uses the
configured beam. The unrestricted beam sentinel is consequently intended for
non-climbing searches rather than a short beam-climbing schedule.

Sparsification parameters require `--tesseract-sparsify-errors`. Its maximum
degree must be at least the base degree, or -1; reactivation limit may be zero,
positive or -1. These settings tune the search rather than changing the
sampling DEM. Limits can produce low-confidence outcomes; `--tesseract-verbose`
prints the underlying search diagnostics.

The effective settings are printed before evaluation. Tesseract-specific flags
are rejected for other decoders, and visualization output is not exposed by
this interface. `--help` lists all options and defaults.

For programmatic diagnostics without verbose search output, call
`dec::Tesseract::set_debug(true)` and read `dec::Result::debug_info` after decoding.
Probabilities refer to retained, possibly merged hyperedges; IDs identify their
representatives in the original flattened DEM. Costs use Tesseract's likelihood
convention and should not be compared directly with BP posterior correction costs.
