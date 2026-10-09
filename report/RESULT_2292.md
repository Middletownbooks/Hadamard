# Hadamard matrix of order 2292 = 4·573 (v = 573 = 3·191)

Date: 2026-10-09. Status: **constructed and verified exactly** (H Hᵀ = 2292·I) [M].

v = 573 was on the open list (orders 2000–3000) used throughout this repo (`report/unattempted_2000_3000.md`).

## Files

| file | content |
|---|---|
| `data/core_v573_b1.txt` | the core: four ±1 sequences of length 572 on cyclic Z₅₇₂, Σ_k PAF_k = −4 off 0, row sums (2, 0, 0, 0) |
| `data/core_v573_b1_paley_FOUND.txt` | the same core in Z₅₂ × Z₁₁ coordinates (entry (g, z) at index 11g + z), as found |
| `data/H2292_packed.npz` | the full matrix: `bits` = `np.packbits(H > 0, axis=1)`, N = 2292. SHA-256 of the packed rows: `fa6f72090debcb07bab9a3cea5f7bcd257f17b67495d39f7770cbcc0e1eeafb8` |
| `scripts/verify_h2292.py` | independent verifier and builder |

## Construction

**Architecture: bordered Goethals–Seidel (b = 1).**
- Four ±1 circulants A, B, C, D of order m = 572, with Σ_k PAF_k(g) = −4 for every g ≠ 0 and row sums (2, 0, 0, 0).
- They go into the Goethals–Seidel array with the back-identity R.
- A border of four rows and four columns is added, S ⊗ 1ᵀ and Tᵀ ⊗ 1 with S, T ∈ H(4), and the 4 × 4 corner is forced
  as −S Σᵀ T / 4, where Σ is the matrix of block row sums.

**The core lies in the Paley(11) + multiplier slice** (`report/paley_design_answer.md` §7.5). This is the slice the
published 419 core occupies, transferred to 573.
- Z₅₇₂ = Z₅₂ × Z₁₁.
- On every Z₁₁-coset, each block is α·δ + u·χ₁₁ + w·ε. That is, Legendre(11)-twisted (u = ±1) or constant (w = ±1) on
  Z₁₁*, with a separate value α at z = 0.
- All labels are invariant under g ↦ 29g on Z₅₂ (order 3), i.e. 20 orbits per block.
- The gadget has K = 21 constant cosets, (6, 6, 2, 7) per block.

**How it was found:**
1. **Gadget design and stage 1:** done by the theory session (`theory/paley_slicepipe.py guided1 573 11 52 29 21 …`).
   It designed cyclotomic w-gadgets ranked by predicted yield, then ran an exact join for α giving the Z₁₁-coset sums
   y. Its output is `theory/logs/paley_slicepipe_stage1_v573_mu29_K21_s0.npz`, and this core uses y-index 20.
2. **Stage 2:** done by `scripts/paley_s2join.py`. It is an exact join for the twisted signs u, with:
   - a Z₄-level class split, which keeps 2^23.7 and 2^27.7 of the pairs on the two sides;
   - a prime-residue hashed 2+2 join.

   It found 2 solutions for this y in 474 s on one core (`logs/s2join_573_paley_slicepipe_stage1_v573_mu29_K21_s0.jsonl`).
   This was the 141st y processed of the 392 recorded, on three cores.
3. **Completion:** the bordered GS array, done by `theory/paley_b1pipe.finish`. That check confirmed H Hᵀ = N·I.

## Verification [M]

`scripts/verify_h2292.py` shares no code with the search pipeline. It:
1. reads the Z₅₂ × Z₁₁ core and maps it to cyclic Z₅₇₂ by CRT;
2. checks Σ PAF = 2288 at 0 and −4 at every other lag (by FFT, rounded, all values integer), with row sums (2, 0, 0, 0);
3. builds the bordered GS matrix from circulants, searching all H(4) pairs S, T for a ±1 corner;
4. checks H Hᵀ = 2292·I exactly. float64 BLAS is exact here, since every entry is an integer of absolute value at most
   2292.

A second process reloaded `data/H2292_packed.npz` and confirmed:
- every entry is ±1;
- H Hᵀ = Hᵀ H = 2292·I;
- an int64 spot check on 200 random row pairs agrees.

Corner found: [[−1, −1, −1, −1], [−1, −1, 1, 1], [−1, 1, −1, 1], [1, −1, −1, 1]].

Reproduce: `python3 scripts/verify_h2292.py`, which takes about 2 minutes and rewrites `data/H2292_packed.npz`.
