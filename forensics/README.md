# forensics/: the authors' original matrices and method fingerprints

The write-up is `report/alpoge_methods_forensics.md`.

## Data (`originals/`)

| file | content |
|---|---|
| `matrix_H{N}_v{v}_{func}_b{b}.txt.gz` | The authors' unnormalised matrices as supplied by the user (header line, then one `+`/`-` row per line), for N = 668, 716, 1244, 1388, 1436, 1772, 1916, 1964. Gzipped; the loaders read `.txt` or `.txt.gz`. |
| `seq_v223_g.txt`, `seq_v283_g.txt` | The authors' record sequences for N = 892 and 1132: A, B, C, D as row 0 of the first block row. |

Not supplied: N = 1676 (v = 419) and 1948 (v = 487). For those, use the repo cores in `data/`.

`.npy` caches are written next to the files on first use and are git-ignored.

## Scripts

| script | what it does |
|---|---|
| `load_originals.py` | Verifies HHᵀ = NI for every original. Compares each with the normalised hoa64 copy in `data/`, and the 668 one with the formal-conjectures matrix. |
| `gbuilder.py` | Infers the `g` (Goethals–Seidel) block layout from H(1244), checks H(1964), and rebuilds 892/1132 from the record files. |
| `census_originals.py` | Multiplier stabiliser, orbit count and fixed point of every original, in natural coordinates (report §3). |
| `halfshift_b1.py` | Half-shift template for the b = 1 cores: σ(y) = x(y)x(y+h), sign changes, null p, arcs and exceptions (report §4.1). |
| `template_parity.py` | Exhaustive parity check: four exact arcs are impossible, and a lone flip must sit at an arc centre (report §4.3). |
| `template_cost.py` | Gaussian cost of the two half-problems, and its scaling with v (report §4.4). |
| `signed_orbits.py` | Signed multiplier relations beyond the stabiliser. Finds the Legendre(11) twist of v = 359 and nothing elsewhere (report §5). |
| `decode668.py`, `blocks668.py` | Decode the Lean-formalised H(668) (google-deepmind/formal-conjectures commit 91ebb29, `FormalConjectures/Wikipedia/Hadamard668.lean`; fetch the file into this directory first) and split it into blocks. |
| `ff_seed_test.py` | Legendre / Sidelnikov seed tests on all repo cores (negative). |
| `pr_structure.py` | Z₂×Z_h P/R split statistics of the repo b = 1 cores. |

"Natural coordinates" means the Goethals–Seidel reversal R has been undone for blocks B, C and D:
- cyclic cores: f(y) = g(−1−y);
- group-developed cores: f(y) = g(−y).
