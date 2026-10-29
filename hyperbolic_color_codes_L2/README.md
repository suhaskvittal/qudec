# L=2 hyperbolic color codes (6.8.8 and 8.8.8)

These are the check matrices and logical operators for some of the modular hyperbolic color codes in

> O. Higgott *et al.*, "Modular hyperbolic codes" (arXiv preprint).

Each code is built on a closed hyperbolic surface. Start from a quotient `G = △(r,s,t)/Γ` of a triangle group, which gives a 6.8.8 or 8.8.8 tiling. Then fine-grain it with side length L = 2 to get a semi-hyperbolic tiling of hexagons and octagons. There is one qubit per vertex, and one X check and one Z check per face (Methods and Supplementary Information, Sec. "Tables of Modular Hyperbolic Codes").

## Codes

The distances are the ones quoted in the paper. A range `a–b` means the distance is bounded below by `a` and above by `b`, and is not known exactly.

| Directory | Tiling | Symmetry group (paper) | `[[n, k, d]]` | Checks | Orientable | Octagon modules | In Fig. 3a | Logical basis |
|---|---|---|---|---|---|---|---|---|
| `8.8.8_L2_D444-32_n64_k8` | 8.8.8 | `△(4,4,4)_32` | [[64, 8, 6]] | 30 | Yes | 2 | simulated | `symplectic_homology_basis` |
| `8.8.8_L2_D444-40_n128_k12` | 8.8.8 | `△(4,4,4)_40` | [[128, 12, 8]] | 60 | Yes | 4 | simulated | `symplectic_homology_basis` |
| `8.8.8_L2_D444-58_n512_k36` | 8.8.8 | `△(4,4,4)_58` | [[512, 36, 12]] | 240 | Yes | 16 | simulated | `symplectic_homology_basis` |
| `8.8.8_L2_D444-109_n1600_k104` | 8.8.8 | `△(4,4,4)_109` | [[1600, 104, 14–16]] | 750 | Yes | 50 | simulated | `tree_cotree_gram_schmidt` |
| `8.8.8_L2_D444-96_n1344_k88` | 8.8.8 | `△(4,4,4)_96` | [[1344, 88, 14]] | 630 | Yes | 42 |  | `tree_cotree_gram_schmidt` |
| `8.8.8_L2_D444-271_n5376_k340` | 8.8.8 | `△(4,4,4)_271` | [[5376, 340, 18–20]] | 2520 | Yes | 168 | extrapolated | `tree_cotree_gram_schmidt` |
| `8.8.8_L2_D444-PSL2-23_n24288_k1522` | 8.8.8 | `△(4,4,4)\|PSL(2,23)` | [[24288, 1522, 22]] | 11385 | No | 759 |  | `tree_cotree_gram_schmidt` |
| `8.8.8_L2_D444-PSL2-25_n31200_k1954` | 8.8.8 | `△(4,4,4)\|PSL(2,25)` | [[31200, 1954, 22–24]] | 14625 | No | 975 | extrapolated | `tree_cotree_gram_schmidt` |
| `6.8.8_L2_D344-18_n384_k20` | 6.8.8 | `△(3,4,4)_18` | [[384, 20, 12]] | 184 | Yes | 12 |  | `symplectic_homology_basis` |
| `6.8.8_L2_D344-22_n768_k36` | 6.8.8 | `△(3,4,4)_22` | [[768, 36, 14]] | 368 | Yes | 24 |  | `symplectic_homology_basis` |
| `6.8.8_L2_D344-32_n1536_k68` | 6.8.8 | `△(3,4,4)_32` | [[1536, 68, 16]] | 736 | Yes | 48 |  | `tree_cotree_gram_schmidt` |
| `6.8.8_L2_D344-45_n3840_k164` | 6.8.8 | `△(3,4,4)_45` | [[3840, 164, 20]] | 1840 | Yes | 120 |  | `tree_cotree_gram_schmidt` |
| `6.8.8_L2_D344-79_n10368_k436` | 6.8.8 | `△(3,4,4)_79` | [[10368, 436, 22]] | 4968 | Yes | 324 |  | `tree_cotree_gram_schmidt` |
| `6.8.8_L2_D344-PGL2-17_n19584_k820` | 6.8.8 | `△(3,4,4)\|PGL(2,17)` | [[19584, 820, 24]] | 9384 | Yes | 612 |  | `tree_cotree_gram_schmidt` |
| `6.8.8_L2_D344-PGL2-23_n48576_k2028` | 6.8.8 | `△(3,4,4)\|PGL(2,23)` | [[48576, 2028, 28]] | 23276 | Yes | 1518 |  | `tree_cotree_gram_schmidt` |

The first four 8.8.8 codes are the ones simulated in Fig. 3a. The `△(4,4,4)_271` and `△(4,4,4)|PSL(2,25)` codes appear there as extrapolated points. All 8.8.8 codes are from Table S "Hyperbolic color codes on the Pareto frontier for 8.8.8 tilings and L=2". All 6.8.8 codes are from Table S "Modular hyperbolic color codes for side length L=2 with octagon modules". (Fig. 3a labels the 50-module code `[[1600, 100, ≤16]]`, but the correct parameters are `[[1600, 104, 14–16]]`.)

## Files in each directory

| File | Contents |
|---|---|
| `H.mtx` | Check matrix `H` (rows = checks/faces, columns = qubits), MatrixMarket coordinate format. **`Hx = Hz = H`.** |
| `H.txt` | Same as `H.mtx`, one check per line: the space-separated, 0-indexed qubits in that check. |
| `check_colors.txt` | One line per check (same order as `H`): the colour of that face, `0`, `1` or `2`. Every qubit is in exactly one check of each colour, and checks of the same colour never overlap. |
| `Lx.mtx` | X logical operators (rows = logical qubits, columns = qubits), MatrixMarket format. |
| `Lz.mtx` | Z logical operators, MatrixMarket format. |
| `code_info.json` | `n`, `k`, `d` (with lower and upper bounds), number of checks, rank of `H`, Euler characteristic, orientability, histogram of check weights by colour, and which logical basis construction was used. |

Conventions:
- Qubits and checks are **0-indexed** in `H.txt`. MatrixMarket files are 1-indexed by definition, and the loaders below convert automatically.
- Checks have weight 6 (hexagons) or 8 (octagons), and every qubit is in exactly 3 checks.
- `H` has `n − k/2` rows but rank `(n − k)/2`. There are exactly two redundancies: for each colour `c`, the checks of colour `c` sum to the all-ones vector.
- **Logical basis:** `Lx @ Lz.T = I (mod 2)` and `H @ Lx.T = H @ Lz.T = 0 (mod 2)`. Together with `rank(H) = (n − k)/2`, this means the rows form a complete basis of `k` logical qubits. The `logical_basis` field in `code_info.json` says which of two constructions was used:
  - `symplectic_homology_basis` (the smallest codes): first pick a symplectic basis of short homologically non-trivial cycles (`α_i`, `β_i`). Restricting each cycle to edges of one colour gives a string operator. Each pair of cycles carries 4 logical qubits, whose X and Z logicals are such single-colour strings on `α_i` and `β_i`. These are low weight.
  - `tree_cotree_gram_schmidt` (larger codes, and the non-orientable ones): homology generators come from a tree-cotree decomposition. Each generator is restricted to colours 0 and 1, and the resulting strings are put in symplectic form by Gram–Schmidt. Here logical qubits `i` and `k/2 + i` have their X and Z supports swapped. Gram–Schmidt adds strings together, so some of these logicals are much heavier than the distance.

  Neither construction is the symmetric automorphism-based basis described in the paper.
- The qubit labelling is arbitrary. In particular, it does not encode the module layout.

## Loading

### Python with SciPy (MatrixMarket)
```python
import numpy as np
import scipy.io
import scipy.sparse as sp

d = "8.8.8_L2_D444-109_n1600_k104"
H  = sp.csr_matrix(scipy.io.mmread(f"{d}/H.mtx"),  dtype=np.uint8)
Lx = sp.csr_matrix(scipy.io.mmread(f"{d}/Lx.mtx"), dtype=np.uint8)
Lz = sp.csr_matrix(scipy.io.mmread(f"{d}/Lz.mtx"), dtype=np.uint8)
colors = np.loadtxt(f"{d}/check_colors.txt", dtype=int)

# Sanity checks
assert not ((H @ H.T).toarray() % 2).any()
assert not ((H @ Lz.T).toarray() % 2).any()
assert ((Lx @ Lz.T).toarray() % 2 == np.eye(Lx.shape[0])).all()
```

### Plain Python (no dependencies)
```python
def load_rows(path):
    with open(path) as f:
        return [list(map(int, line.split())) for line in f]

def load_mtx_rows(path):
    """Read a MatrixMarket coordinate file into a list of rows (0-indexed column lists)."""
    with open(path) as f:
        lines = [l for l in f if not l.startswith("%")]
    nrows, ncols, nnz = map(int, lines[0].split())
    rows = [[] for _ in range(nrows)]
    for l in lines[1:]:
        i, j, _ = map(int, l.split())
        rows[i - 1].append(j - 1)
    return [sorted(r) for r in rows]

checks = load_rows(f"{d}/H.txt")    # checks[i] = qubits in check i
x_logicals = load_mtx_rows(f"{d}/Lx.mtx")
z_logicals = load_mtx_rows(f"{d}/Lz.mtx")
colors = [int(c) for c in open(f"{d}/check_colors.txt")]

# Dense numpy matrix, if wanted:
import numpy as np
n = 1 + max(q for row in checks for q in row)
Hd = np.zeros((len(checks), n), dtype=np.uint8)
for i, row in enumerate(checks):
    Hd[i, row] = 1
```

### Other languages
- **Julia:** `using MatrixMarket; H = MatrixMarket.mmread("H.mtx")`
- **MATLAB:** `mmread.m` from the NIST Matrix Market site.
- **C/C++:** the `.mtx` format is a header, then a `rows cols nnz` line, then one 1-indexed `row col value` triple per line. `H.txt` is even simpler to parse.

### Stim (optional)
To build stabilizers or a code-capacity check in Stim, each row of `H` gives an X stabilizer (`stim.PauliString` of `X` on those qubits) and a Z stabilizer (`Z` on the same qubits).
