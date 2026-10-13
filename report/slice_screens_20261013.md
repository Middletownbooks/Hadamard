# Exact slice screens, 2026-10-13: 955, 971, 917, 537 ×59 dead; what is left

Branch `claude/focused-planck-j1645f`. Scripts: `theory/gp_*.py`, `theory/slice_compress_cpsat.py`,
`theory/q971_stage1*.py`. Logs: `logs/gp_screen/`.

Tags: [P] proved here, [M] exact computation (script named), [E] Gaussian estimate.

## 1. Three exact screens

Setting: a b = 1 bordered GS core on G = Z_m (m = v − 1). In a generalized-Paley (GP) slice, Z_F is fixed pointwise
and the slice acts on Z_q by the index-s subgroup (f = (q − 1)/s). The stage-1 coset sums are
y = α + f·σ, with σ a sum of s signs.

1. **Lag 0 + rows + Z₂** (`theory/gp_z2screen.py`) [P, M].
   - Every nontrivial character χ of Z_F gives Σ_k |ŷ_k(χ)|² = T_J(0) + 4q = 4v. The trivial character gives
     Σ_k r_k² = 4, so the rows are (±2, 0, 0, 0).
   - The screen checks, for each lag-0 composition, whether the large letters can be placed over the blocks and their
     two parity halves so that the rows are (2, 0, 0, 0) **and** Σ A_k² = 4v (A = alternating sum).
   - Controls pass: v = 41, 235 (H(940) exists in the slice), 523 (stage-1 solutions known).

2. **α-congruence** (`theory/gp_alpha_screen.py`) [P, M].
   - Write y = β + 2f·j with β = α + f·[s odd].
   - The cross term 2f·Σ(β j′ + j β′) is 2f times an even number. So Σ_k PAF(β_k)(h) ≡ −4q (mod 4f) for h ≠ 0, and
     the rows of β are ≡ r_k (mod 2f).
   - It is a finite check over ±1 sequences on Z_F.

3. **Compressed stage 1 of any slice** (`theory/slice_compress_cpsat.py`) [P, M].
   - Take any multiplier group M = ⟨μ⟩ on Z_m and any Hall factor Z_K. The compression y(b) = Σ_{g≡b (mod m/K)} x(g)
     is M-invariant.
   - Its value at b ranges over the signed sums of the sizes of the M_b-orbits that partition the fibre (M_b = the
     stabiliser of b). These are independent across the M-orbits of Z_{m/K}.
   - It satisfies Σ_k PAF(y_k) = 4m − 4(K − 1) at lag 0 and −4K elsewhere, with rows (2, 0, 0, 0).
   - It is solved by CP-SAT; INFEASIBLE is a proof.
   - **Composite slices are intersections [P].** If μ acts nontrivially on two primes, the slice's compression along
     each prime is the stage 1 of a GP slice (possibly with fewer free letters).
     - 537 ×59 compresses along Z₆₇ to the sextic GP stage 1 on Z₈ (same order-11 subgroup as ×9) with extra
       ×3-invariance. That stage 1 is exactly empty.

## 2. Results

| v | slice | structure | Gaussian [E] | exact result |
|---|---|---|---|---|
| 955 | ×289 | Z₁₈ × quartic Z₅₃ (f = 13) | +61.8 | **empty [P, M]**: rows fail for both lag-0 compositions (below) |
| 971 | ×81 | Z₁₀ × index-8 Z₉₇ (f = 12) | +46.1 | **empty [P, M]**: §3 |
| 917 | ×17 | Z₄ × index-12 Z₂₂₉ | +54.1 | **empty [M]**: the one lag-0 composition fails rows + Z₂ |
| 907 | ×85, ×91 | Z₆ × index 10 and 6 on Z₁₅₁ | +66.6, +32.1 | **empty [M]**: no lag-0 composition (matches CP-SAT INFEASIBLE for ×85) |
| 823 | ×73 | Z₆ × index-8 Z₁₃₇ | +52.4 | **empty [M]**: no lag-0 composition |
| 971 | ×221, ×101 | Z₁₀ × index 6 and 4 on Z₉₇ | +28.4, +9.3 | **empty [M]**: no lag-0 composition |
| 789, 917, 523, 659, 739, 947 | the remaining s ≤ 4 GP slices with E > 0 | | +1.5 … +11.5 | **empty [M]** |
| 537 | ×59, ×15 | fixes Z₂; index 6 on Z₂₆₈ | +18.5, +11.5 | **empty [P, M]**: compressed stage 1 INFEASIBLE (0.2 s) |
| 599 | ×27 | Z₂₆ × Paley(23) | +27.3 | passes rows + Z₂, but already dead by the w-screen (Theorem 5, `report/tractability_after_2292.md`) |
| 523 | ×181 | Z₁₈ × quartic Z₂₉ | +92.1 | alive; production run elsewhere |
| 947 | ×45 | Z₂₂ × cubic Z₄₃ (f = 14) | +30.8 | **empty [P]**: s = 3 is odd, so y = α + 14·(odd) and \|y\| ≥ 13; lag 0 needs Σy² = 3616 over 88 entries (minimum 88·169). An earlier run of `gp_z2screen` passed it because `q523_stage1.compositions` assumes ±1 is a letter; fixed in `gp_z2screen.py` |
| 971 | ×79 | −1 on Z₅; index 16 on Z₉₇ | +17.7 | **empty [M]**: compressed stage 1 INFEASIBLE (31 s) |
| 947 | ×49 | Paley(11) × order 3 on Z₈₆ (the H(2292) shape) | +37.3 | **compression along Z₄₃ FEASIBLE** (verified y in `logs/gp_screen/c947_49_K43_y.json`); along Z₁₁ UNKNOWN at 1200 s; w-screen alive. Paley split [E]: stage 1 +33 … +46, stage 2 −19 … −22 per y (573: −9.6), core +14 … +24 (573: +23) |
| 947, 859 | ×21, ×289 | composite | +28.3, +27.4 | compressions UNKNOWN at 600 s |
| 739 | ×215 | fixes Z₂; −1 on Z₉, index 8 on Z₄₁ | +61.9 | survives; compressed stage 1 along Z₄₁ (Z₁₈, 10 orbit letters per block, 18 letters) UNKNOWN after 900 s and 1200 s; along Z₉ (Z₈₂) UNKNOWN after 1200 s |
| 647 | ×353 | fixes Z₂; index 4 on Z₁₇, index 6 on Z₁₉ | +36.2 | survives; both compressions (Z₃₈, Z₃₄) UNKNOWN after 1200 s |
| 971 | ×469 | −1 on Z₅, index 8 on Z₉₇ | +36.2 | compressed stage 1 on Z₁₀ UNKNOWN after 600 s |

Logs: `logs/gp_screen/composite_compress_cpsat.jsonl`, `composite_batch2.jsonl`, `c739_215_K41_seed7.json`, `c947_49_*.json`.

**Proof for 955 [P].**
- The letters are ±1, ±25 = ±(26 − 1), ±27 = ±(26 + 1), ±51 = ±(52 − 1) and ±53 = ±(52 + 1). A signed big letter
  is 26K ± 1.
- A block with j big letters has row sum 26·K_tot + e + u, with |e| ≤ j and |u| ≤ 18 − j. For |r| ≤ 2 this forces
  K_tot = 0.
- Lag 0 (Σ y² = 3608 over 72 entries) leaves two compositions: (25, 27, 27, 27, 27) and (27, 53).
  - In the first, every big letter carries ±1 unit of 26, so each block needs an even number of them. The total of 5
    is odd.
  - In the second, the units are 1 and 2, and no signed subset per block sums to 0.

So no stage-1 y exists. The Gaussian had given +61.8; CP-SAT had returned UNKNOWN at 1200 s.

## 3. 971 ×81 [P, M]

- Rows + Z₂ leave one lag-0 composition, {23³, 25³, 1³⁴}. So y = α + 24j, with j ∈ {0, ±1} on 6 cells and
  |y| = 25 exactly where αj = +1 (three cells).
- PAF(α + 24j) = PAF(α) + 24C + 576·PAF(j), and C is even. So mod 48, Σ_k PAF(α_k) ≡ −388 ≡ −4. Since
  |Σ PAF(α)| ≤ 40, the α_k form an exact b = 1 core on Z₁₀ (rows exactly (2, 0, 0, 0), hence Σj = 0 per block).
- `theory/q971_stage1.py` runs over every α-core up to per-block shift, reversal and negation (99 classes) and over
  every j: **no solution**. The planted self-test recovers 20/20. Log: `logs/gp_screen/q971_stage1_exact.json`.
- Independent check: `theory/q971_stage1_cpsat.py` (CP-SAT on (α, j) with no reduction): **INFEASIBLE** in 300 s.

## 4. Lessons

- Gaussian slice estimates of +46 to +67 sat on exactly empty slices five times today, after 739 ×127 (+57) and the
  537 sextic slice (+7). Every candidate must pass the exact screens above before any investment.
- Rows + Z₂ and lag 0 cost seconds. The α-congruence mod 4f is the 523 session's "mod-8 XOR" condition in general
  form, and it explains why naive stage-1 joins fail.
