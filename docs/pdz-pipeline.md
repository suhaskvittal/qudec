# Configurable PDZ pipeline

Contributor: OpenAI Codex.

`qudec pdz` composes a first decoder, PDZ filtering, and a residual decoder.
The default first decoder is uncorrelated PyMatching and the default residual
decoder is correlated PyMatching:

```sh
mpirun -np 64 ./build/qudec pdz 9 --pdz-initial-decoder pymatching
mpirun -np 64 ./build/qudec pdz 9 --pdz-initial-decoder uf
```

Both commands automatically include opposing-basis detectors. All existing
Monte Carlo and rare-event options apply. Use `--uf-weight-scale` to adjust
UF weight quantization. `--pdz-residual-decoder libra` selects Libra instead
of correlated PyMatching. The existing `pdz_libra` command uses the same
generic composition and retains its specialized histograms and recording.

PDZ uses the first decoder's matched detector pairs and their observable
parities. UF reports paths through its actual peeled correction forest;
branching paths are paired deterministically. This preserves UF's prediction,
but its pairing can differ from MWPM's. Enabling reporting adds path extraction
work to the first decoder's latency. Reporting currently requires UF's merged
channel probabilities to be at most one half, as in our SI1000 model.

The bypass rule is unchanged: if PDZ removes at least
`--pdz-removal-threshold` of the original fired bits (default 0.5), return
the first decoder's prediction. Empty residuals always bypass. Set the
threshold to 1.0 to call the residual decoder on every nonempty residual.
On those calls, XOR the PDZ correction with the residual decoder's prediction.

Output reports the residual call count, bypasses, frame changes, each decoder's
latency, full pipeline mean latency, and residual mean latency weighted by its
call fraction. These are sample averages; rare-event conditional sample
averages are not physical-frequency weighted latency estimates.

## C++ interface

```cpp
#include "decoder/pdz.h"
#include "decoder/pymatching.h"
#include "decoder/uf.h"

dec::UnionFind initial(dem, {.report_matching_data=true});
dec::PyMatching residual(dem, {.correlated=true});
dec::PDZDecoder pipeline(dem, &initial, &residual, {}, 1.0);
auto result = pipeline.decode(syndrome, observables);
```

`PDZDecoder<Initial, Residual>` takes two non-owning decoder pointers. Both
decoders and the DEM must outlive the pipeline. The initial decoder must
report matched endpoints. Each decoder supplies `decode`, `print_stats`,
`mpi_accumulate`, `num_detectors`, and `num_observables`. The template does
not require inheritance or virtual dispatch. Use one instance per worker.
`decode_detailed` also returns the initial prediction, PDZ correction,
remaining syndrome, bypass flag, and final prediction.
