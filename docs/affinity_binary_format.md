# Affinity analysis binary format

Author: OpenAI GPT-6 (24 September 2026).

This document describes the version 1 output of `analyze_syndrome`. All
multibyte integers and IEEE-754 float64 values are little-endian. A `.bin.xz`
file is an LZMA/XZ-compressed stream of these bytes; `.bin` contains the same
bytes without compression. The file has one header followed by shot records
until end of file. It contains no JSON, shot indices, seed, or shot-count footer.

## Header (23 bytes)

| Field | Bytes | Meaning |
| --- | ---: | --- |
| Magic | 8 | ASCII `QDAFF001` |
| Distance | 1 | Rotated surface-code distance |
| Rounds | 1 | Syndrome rounds |
| Physical error rate | 8 | SI1000 `p`, float64 |
| MPI ranks | 2 | Unsigned integer |
| Detector count | 2 | Number of physical circuit detectors; the boundary ID equals this count |
| Flags | 1 | Bit 0: odd-weight boundary rule; bit 1: opposite-basis detectors included |

The current writer always sets flag bit 0. It constructs Z-memory circuits,
uses an undecomposed DEM for affinity, and gives PyMatching a decomposed DEM.

## Shot record

| Field | Bytes | Meaning |
| --- | ---: | --- |
| Physical Hamming weight `h` | 2 | Number of fired circuit detectors |
| Boundary indicator `b` | 1 | 1 for odd `h`, otherwise 0 |
| PyMatching logical error | 1 | 1 if its prediction differs from sampled truth |
| Detector IDs | `2(h+b)` | Unsigned IDs, physical IDs in ascending order, boundary last if present |
| Effective partner counts | `8(h+b)` | One float64 value for each detector ID, in the same order |

Thus each record uses `4 + 10(h+b)` bytes. For affinity row `i`, the stored
effective partner count is `(sum_j a_ij)^2 / sum_j a_ij^2`, excluding `j=i`.
If the denominator is zero, the stored count is zero. The full pairwise
affinity matrix cannot be reconstructed from these counts.

The reader determines the next record from `h` and `b` and continues until
end of file. XZ checksums detect compressed-stream corruption or truncation;
an uncompressed `.bin` file has no independent completion check.
