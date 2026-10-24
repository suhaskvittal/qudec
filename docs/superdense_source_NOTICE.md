<!-- Author: OpenAI Codex. -->

The superdense color-code implementation adapts the source release accompanying
Noah Shutty, *Denser Planar Color Codes*, arXiv:2609.21376v3:
<https://doi.org/10.5281/zenodo.22820614>.

The release's `NOTICE.md` and `LICENSE` explicitly license its sources under
Apache-2.0. A copy is preserved in [superdense_source_LICENSE](superdense_source_LICENSE).
The Zenodo metadata separately lists CC-BY-4.0; this adaptation follows the
explicit source-release license.

`src/superdense_488_rules.h` is a lossless C++ transcription of the 90 scheduling
classes in `src/vendor/_triangle488/reference/templates.json`. Geometry,
instrument algebra, Bell-polarity optimization and native scheduling in
`src/superdense_color_code_circuits.cpp` adapt `_triangle488/geometry.py`,
`reference/cx_rule.py`, `reference/ideal_css.py`, `candidate.py`, `native.py`,
and `memory.py` from that release. These sources are attributed to Noah Shutty's
release. The C++ adaptation and integration are authored by OpenAI Codex.

The 6.6.6 layout and CNOT contact schedule adapt
`clorco/color_code/_superdense_planar_color_code_circuits.py`,
Copyright 2023 Google LLC, Apache-2.0, distributed in the same bundle from
Chromobius revision `cfcec5333856f2e7a07f0a4fced01dd9d2d418a3`.
Original upstream: <https://github.com/quantumlib/chromobius>.

Changes include translation to C++, integration with qudec's operation/noise
policies and component distributions, six-coordinate detector metadata,
optional opposite-basis detectors, configurable complete-cycle memory
boundaries, and command-line selection. No upstream decoder is incorporated.
