# H(3676), v = 919: a b = 0 core in a 3-cycle-twisted multiplier slice

Date: 2026-10-16. Branch `claude/focused-planck-j1645f`.
Scripts: `scripts/b0_open_screen.py`, `scripts/b0_twist3_search.py`, `scripts/verify_h3676.py`, `scripts/gs_twist.py`.
Logs: `logs/b919/`. Bundle: `pub/H3676/`.

**Tags.** [M] exact computation (script named); [P] proved here; [E] Gaussian estimate (estimator and calibration
named).

---

## 0. Bottom line

1. **H(3676) is constructed and verified [M].** v = 919 was on the open list (`report/prompt_alpoge_design.md`). It is
   also listed as unknown in `data/sage_unknown_orders_le4000.json`.
   - **Core.** Four ±1 circulants A, B, C, D of length 919 with Σ_k PAF_k = 3676·δ and row sums (1, 35, 35, 35).
     Files: `data/core_v919_b0_twist51_0.txt` (bundle core) and `data/core_v919_b0_twist51_1.txt` (a second,
     distinct core).
   - **Matrix.** The Goethals–Seidel array of the four circulants (`g` layout, back-identity R).
   - **Verification.**
     - `scripts/verify_h3676.py`: integer PAF, the slice symmetries, and H Hᵀ = 3676·I in int64. Both cores pass
       (`logs/b919/verify_h3676_core{0,1}.json`).
     - The standalone `pub/H3676/verify.py` checks the sequences and the published normalized matrix
       (`logs/b919/pub_verify_h3676.txt`).
     - SHA-256 of the packed rows (`data/H3676_packed.npz`):
       `cfed1db980e75b703e982469364cf41d74944dce5ef1e879f3099ebbdb27bc34`.
2. **The slice [P]** (§1).
   - U = ⟨8⟩ is the subgroup of order 51 of F₉₁₉^*.
   - A is U-invariant: constant on the 18 U-orbits of Z₉₁₉^*.
   - B is invariant under U³ = ⟨8³⟩ (order 17, 54 orbits).
   - C(x) = B(8⁻¹x) and D(x) = B(8⁻²x).
   - So the generator fixes A and cycles B → C → D → B: the 3-cycle twist of `scripts/gs_twist.py`.
   - In the estimator's terms: 19 orbits, 9 equations, log₂E = 24.4, n = 49.6 bits [E].
     **This is the cheapest slice found for any open order**: below T(131) (n = 51.0, `report/tmat_routes.md` §9).
3. **Why only 919 [P].** The block sums satisfy s_A² + 3s_B² = 4v, with s_A ≡ ±1 (mod d) and s_B ≡ ±1 (mod d/3).
   - At v = 919 the identity is 4·919 = 1 + 3·35². Since 35 ≡ 1 (mod 17), d/3 = 17 is allowed, so d = 51.
   - Among the 17 prime open orders, 4v − 1 = 3s² holds only for v = 919.
4. **Search [M]** (§2).
   - Method: a plain 2-list search. The complete A-table (48,620 sequences with row sum 1) is matched against random
     B with row sum 35. Each B is evaluated by orbit-level quadratic forms, about 3·10⁵ B per second on 3 cores.
   - Result: two distinct cores after about 2^27 B, in about 8 minutes on 3 cores.
   - Prediction: 2^28.35 B per core, so the estimate was about 1.3 bits low. This is in line with its validation on
     the same family (§3).
5. **How it was missed.**
   - The b = 0 cyclic screens stopped at 2000–3000 (`logs/nbits.jsonl`, `report/unattempted_2000_3000.md`).
   - The twisted-GS estimator `gs_twist.py` was written for 571 and never run above 3000.
   - The 3000–4000 work went to b = 1 sweeps on Z_{v−1}, b ≥ 3, the staged Paley pipelines and T-matrices.
   - The objection that stopped the 571 twist (a lopsided 2-list problem, 2^54.6 B per hit) does not apply here. The
     slice is small: 2^28 B per hit.

## 1. The slice [P]

- **Twisted multiplier.** Let u = 8, of order 51 mod 919 (51 | 918 = 2·3³·17). The ansatz is X_{π(i)}(u·x) = X_i(x),
  with π the 3-cycle (B C D) and A fixed. That gives:
  - A(ux) = A(x);
  - C(ux) = B(x), D(ux) = C(x) and B(ux) = D(x).

  Hence B(u³x) = B(x), C(x) = B(u⁻¹x) and D(x) = B(u⁻²x).
- **Equations.** PAF_C(c) = PAF_B(u⁻¹c) and PAF_D(c) = PAF_B(u⁻²c). So ΣPAF = 4v·δ becomes, for each of the
  918/(2·51) = 9 classes of ⟨U, −1⟩ on lags,

  PAF_A(c) + PAF_B(c) + PAF_B(u⁻¹c) + PAF_B(u⁻²c) = 0.

  Both sides are constant on classes: PAF_A is U-invariant, and the B-sum is invariant because B is U³-invariant.
- **Row sums.** Write s_A = A(0) + 51·k_A and s_B = B(0) + 17·k_B. Then s_A² + 3s_B² = 3676.
  - The integer solutions are (±1, ±35) and (±53, ±17).
  - ±53 ≢ ±1 (mod 51) and 17 ≢ ±1 (mod 17), so the only block-sum pattern is (±1, ±35, ±35, ±35).
  - Normalised: A(0) = 1 with 9 + and 9 − orbits, and B(0) = 1 with 28 + and 26 − orbits.
- **Freedom.** 2^{4·18+2} = 2^74 labelled configurations. The search space is C(18, 9) = 48,620 A's times
  C(54, 28) ≈ 2^50.8 B's.

## 2. Search [M] (`scripts/b0_twist3_search.py`)

1. Tabulate every A: 48,620 entries, class values at the 9 lag representatives, with a hashed key.
2. Draw B uniformly with its row sum fixed. Compute T_B(c) = PAF_B(c) + PAF_B(u⁻¹c) + PAF_B(u⁻²c) as
   βᵀQ_cβ, where β = (B(0), b₁, …, b₅₄) and Q_c is the symmetrised sum of the three shift-incidence matrices EᵀS_ℓE.
   One BLAS product per batch.
3. Look up −T_B. Every hit is re-checked on all 919 lags, the GS array is built, and H Hᵀ = 4v·I is checked before
   saving.

**Control [M].** v = 61, d = 15 has the same row-sum family (4·61 = 1 + 3·9²).
- The exhaustive version (`--exhaustB`) finds 672 normalised cores, every one verified
  (`logs/b919/control_v61_d15_exhaustive.jsonl`).
- 672 × 4 sign choices = 2,688, which is exactly `gs_twist.brute(61, 15, 3-cycle)`.
- The fast quadratic-form path and the FFT path give the same count.

**Run at 919** (`logs/b919/search_seed2{1,2,3}.jsonl`):
- three workers, seeds 21, 22, 23, batches of 32,768;
- cores at 435 s (seed 21) and 477 s (seed 23), after about 2^27 B in total;
- the workers were stopped after the second core.

## 3. Calibration of the estimate [M]

| case | exact | estimate (`gs_twist.estimate`) |
|---|---|---|
| (37, 9, 3-cycle) | 2^11.8 brute force | 2^10.9 (`report/unattempted_2000_3000.md`, update) |
| (61, 15, 3-cycle) | 2^11.4 brute force; 2,688 by this search | 2^9.6 |
| (919, 51, 3-cycle) | 2 cores in ≈ 2^27 B, so ≈ 2^27 B per core | 2^28.35 B per core, i.e. 1.3 bits low (from two events) |

All three belong to the same row-sum family, 4v = 1 + 3s².

## 4. The b = 0 screen at the prime open orders (`scripts/b0_open_screen.py`, `logs/b919/b0_open_screen.jsonl`)

**Scope.**
- Every odd d with 15–130 orbits, both untwisted and with the 3-cycle twist where 3 | d.
- Exact block-sum feasibility by ILP, then the Gaussian.
- Orders covered: 571, 599, 659, 739, 751, 859, 883, 919, 947 and 991.
- 643, 719, 823, 839, 863, 907 and 971 have no odd d in range.

Feasible slices, cheapest first:

| v | d | twist | orbits | classes | block sums | log₂E | n (bits) |
|---|---|---|---|---|---|---|---|
| **919** | **51** | **3-cycle** | 19 | 9 | (1, 35, 35, 35) | 24.4 | **49.6** (built) |
| 919 | 27 | 3-cycle | 35 | 17 | (−53, −17, −17, −17) | 46.9 | 91.1 |
| 571 | 15 | 3-cycle | 39 | 19 | (31, 21, 21, 21) | 59.4 | 94.6 (the known lopsided lead) |
| 883 | 21 | none | 43 | 21 | (−41, 1, 1, 43) | 60.1 | 111.9 (built: H(3532), `report/b0_883_answer.md`) |
| 919 | 17 | none | 55 | 27 | (1, 35, 35, 35) | 76.5 | 143.5 |
| 859 | 13 | none | 67 | 33 | (−25, 1, 1, 53) | 94.6 | 173.4 |
| 991 | 15 | 3-cycle | 67 | 33 | (61, −9, −9, −9) | 89.6 | 176.4 |
| 859 | 11 | none | 79 | 39 | (−43, 23, 23, 23) | 111.6 | 204.4 |
| 739 | 9 | 3-cycle / none | 83 | 41 | | 120.3 / 121.8 | 209.7 / 210.2 |
| 659 | 7 | none | 95 | 47 | | 143.4 | 236.6 |
| 883 | 9 | 3-cycle / none | 99 | 49 | | 137.3 / 138.8 | 256.7 / 257.2 |
| 919 | 9 | 3-cycle / none | 103 | 51 | | 141.4 / 142.9 | 268.6 / 269.1 |
| 571 | 5 | none | 115 | 57 | | 179.5 | 280.5 |
| 991 | 9 | 3-cycle / none | 111 | 55 | | 149.4 / 150.9 | 292.6 / 293.1 |
| 883 | 7 | none | 127 | 63 | | 178.0 | 330.0 |

**Row-sum infeasible:** 28 (v, d, twist) combinations; the full list is in the log.

**Notes.**
- Untwisted b = 0 costs equal the T-form costs of `report/tmat_routes.md` §10, as §2 there predicts. Examples:
  883/21 gives 111.9 in both, and 919/9 gives 269.1 in both.
- **883/21 (n = 111.9) is the next b = 0 candidate.** Update: built by this route, H(3532)
  (`report/b0_883_answer.md`). It is untwisted, so the flat-pool k-tree of
  `report/REPORT.md` §4–5 applies directly, with four independent lists.
  - That method found two new cores at v = 487 with n = 130 (`report/v487_new_solutions.md`).
  - The flat restriction gained 1.7–1.9 bits per class at 283 and 487. That would put 883/21 near n′ ≈ 75 [E,
    unmeasured].

## Reproduce

```
python3 scripts/b0_open_screen.py 571,599,659,739,751,859,883,919,947,991
python3 scripts/b0_twist3_search.py 61 15 1 -9 OUT --exhaustB          # control: 672 cores
python3 scripts/b0_twist3_search.py 919 51 1 35 OUT --seed 21 --batch 32768
python3 scripts/verify_h3676.py data/core_v919_b0_twist51_0.txt
cd pub/H3676 && python3 verify.py
```
