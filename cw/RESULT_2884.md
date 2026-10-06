# Hadamard matrix of order 2884 = 4·7·103 (v = 721)

Date: 2026-10-06. Status: **constructed and verified exactly** (H Hᵀ = 2884·I, integer arithmetic) [M].

## Construction

This is the Cooper–Wallis theorem: T-matrices of order t together with Williamson matrices of order w give a Hadamard
matrix of order 4tw.

- **t = 103.** The T-matrices are four disjoint ternary circulants whose supports partition Z₁₀₃, with
  Σ PAF = 103·δ. They are taken from the order-2060 construction by L. Schneider
  ([schneiderlo/hadamard-2060](https://github.com/schneiderlo/hadamard-2060), MIT;
  added to SageMath in sagemath/sage#42876, 2026-09-30).
  - That construction is t = 103 with w = 5, giving 2060 (v = 515).
  - The data is a 35-symbol orbit seed under the multiplier subgroup {1, 46, 56} ≤ (Z/103)*.
- **w = 7.** The Williamson quadruple was found by brute force in `cw_build.py`; Williamson matrices of order 7 are classical:

  ```
  A = (1, 1, 1,-1,-1, 1, 1)   B = (1, 1, 1,-1,-1, 1, 1)
  C = (1, 1,-1, 1, 1,-1, 1)   D = (1,-1, 1,-1,-1, 1,-1)
  ```

- **Since 721 = 7·103, the existing T-matrices of order 103 immediately give order 2884.** v = 721 was on the open list
  used in `report/bge3_construction.md` (4v between 2000 and 3000).

## Verification

`scripts/cw_build.py` is a separate implementation that does not use Sage.

**What it checks:**
- **T-matrices:** values in {0, ±1}, supports partition Z_t, and Σ PAF = t·δ.
- **Williamson quadruple:** Σ PAF = 4w·δ.
- **Assembly:** X, Y, Z, W are the Cooper–Wallis sums of T_i ⊗ (Williamson); they are placed in the Goethals–Seidel
  array with R the reversal on Z_t × Z_w; finally H Hᵀ = nI is computed exactly.

**Control:** with w = 5 the same code rebuilds an order-2060 matrix and verifies it.

**Results:**
- `python3 scripts/cw_build.py 7` gives `HHt_eq_nI: true` (about 1 min).
- The stored matrix is `cw/H2884_packed.npz` (np.packbits of H < 0, row-major). It was reloaded from disk and
  re-verified (H Hᵀ = 2884·I). SHA-256 of the packed bits: `a2ce2b665e44379645a1ba815950d32722bbce1a9cfea3174389dc3b2a510513`.

## Caveats

- **Novelty:** checked by web search on 2026-10-06, which found no prior claim of order 2884. The open-order list
  this rests on is the repo's (Djoković's list for v < 750); I have not re-derived it against the primary source.
- **Credit:** the essential ingredient (T-matrices of order 103) is Schneider's. The step here is only that 721 is
  also a multiple of 103.
