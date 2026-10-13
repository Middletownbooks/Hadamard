# Next design: exact status of the live slices, a new F_{p²} family, and the best remaining route (647 ×183)

Date: 2026-10-13. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_next_design.md`.
Scripts: `theory/nd_*.py`. Logs: `logs/nd/`.

**Tags.** [M] exact computation (script named); [P] proved here (proof in the text); [E] Gaussian estimate (the
estimator is named, and so is its calibration); [L] literature recollection, unchecked (arXiv and OEIS are blocked).

---

## 0. Bottom line

1. **(A) No new Hadamard matrix this session.**

2. **(B) Not reached.** No route is validated end to end at ≤ 2^45 per core.
   - **Best costed route: 647 ×183, i.e. H(2588) (§5).**
     - Both of its compressed stage 1s are exactly non-empty [M].
     - Once a stage-1 solution y_B on Z₃₄ is fixed, the rest of the core is a problem shaped like 523's stage 2, and
       slightly smaller. Its constraint budget, against 523's:

       | | 647 ×183 (given y_B) | 523 |
       |---|---|---|
       | all constraints | 126.4–127.4 bits | 146.2 bits |
       | S-level | 72.6–72.9 bits | 83.3 bits |
       | middle level | 22.8–22.9 bits | 31.8 bits |
       | D-deficit D_def | 25.4–26.7 bits | 26.1 bits |
       | completions per y_B | +15.4 … +18.8 | +34 per y |

     - The 523 route ported to this slice projects to about 2^41–2^44 operations per core [E]. That is not yet
       calibrated for this shape, and the solver is not built.
     - First computation: §5.4.
   - **Other routes, costed and out of reach [E]:**
     - 971 ×469: ≥ 2^60 (§3);
     - 947 ×49: ≥ 2^65 (§6);
     - 883 on Z₁₈ × F₄₉: no route (§4).

3. **(C) Exact obstructions [P, M]** (§1, §2, §7).
   - **739 ×215 is exactly empty**: its stage 1 has no solution. This kills the most valuable candidate below 3000.
   - Also exactly empty:
     - 947 ×21, ×47 and ×23;
     - 823 ×59;
     - 907 ×155;
     - 933 ×37 (lag 0);
     - 883 on Z₁₈ × F₄₉: the combined slices (7,4), (7,20), (17,4), (7,2) and (7,10);
     - 751 on Z₃₀ × F₂₅: the combined slices (17,2) and (17,14).
   - So every slice of 947 except ×49 is dead.

4. **New structure [P, M].**
   - **An exact compressed stage-1 solver** (`theory/nd_s1_sd.py`), built from two ingredients:
     - the α-congruence mod 4f₀;
     - the exact S/D split on Z_F = Z₂ × Z_n.

     It decides in seconds what CP-SAT left UNKNOWN for 900–3600 s. Validation:
     - it matches an independent brute force on 5 controls, including a count of 1,750,702;
     - planted tests recover 20/20, 10/10 and 10/10.
   - **Uniform cyclotomy over F_{p²}** (§4). On G = Z_F × F_q with q = p² and s | p + 1, the Gauss periods are
     integers. Stage 2 then splits exactly into s integer "class systems" on Z_F, coupled only through per-position
     sums. The splitting was checked numerically to 10⁻¹⁵.
     - This opens noncyclic groups for 883, 751 and the F₉ orders. Their sweeps are in §4.3.
   - **Calibrated tools for any cyclic slice:**
     - a stage-2 estimator (`theory/nd_stage2_est.py`): within 0.7 bit of the 523 record; on the 573 record it
       predicted 0.6 solutions where 2 were found;
     - a level-budget tool with exact h_D (`theory/nd_levels.py`): it reproduces the 523 S, D, Z₃, Z₁ and h_D
       figures.

5. **(D) Ranked candidates with their next computations: §8.**

---

## 1. Exact obstructions found this session

| v | slice | structure | Gaussian [E] | exact result | how |
|---|---|---|---|---|---|
| **739** | **×215** | −1 on Z₉ (y symmetric on Z₁₈); index 8 on Z₄₁ (f = 5) | +61.9 whole, +11.2 stage 1 | **stage 1 empty** | §2: α-congruence mod 20 leaves 21 α-quadruples (canonical), S/D split gives 12 S-solutions, 0 completions [M] |
| 947 | ×21 | −1 on Z₁₁ (y symmetric on Z₂₂); order 7 on Z₄₃ (f = 7) | +28.3 whole, −8.1 stage 1 | **empty** | `nd_s1_sd.py 946 21 43`: 40 α-quadruples, 0 S-solutions [M]; planted 10/10 |
| 947 | ×47 | Z₂₂ compression: ×5 on Z₁₁ (QR), order 7 on Z₄₃ | +2.0 | **empty** | 4 α-quadruples, 32 S-solutions, 0 completions [M] |
| 947 | ×23 | Z₂₂ × Paley(43) | +1.7 | **empty** | [P] below |
| 823 | ×59 | Z₆ (−1 on Z₃) × index 8 on Z₁₃₇ (f = 17) | +7.9 | **empty** | no α-quadruple mod 68 [M] |
| 907 | ×155 | Z₆ (−1 on Z₃) × index 10 on Z₁₅₁ (f = 15) | +13.6 | **empty** | no α-quadruple mod 60 [M] |
| 933 | ×37 | Z₄ × index 8 on Z₂₃₃ (f = 29) | +30.6 | **empty** | [P] lag 0: letters ±1 + 58t; 57² − 1 = 3248 > 2800 − 16, so only ±1 fits, but Σy² must be 2800 over 16 entries |
| 883 | Z₁₈ × F₄₉, (7,4) and (7,20) | ×7 on Z₁₈ (order 3) × index 4 on F₄₉* | +64.3 whole, +22.6 stage 1 | **stage 1 empty** | 565,348 α-quadruples, 863,136 S-solutions, 0 completions (1344 s) [M] |
| 883 | Z₁₈ × F₄₉, (17,4), (7,2), (7,10) | −1 or ×7 on Z₁₈ with index 4 or 2 on F₄₉* | +37.4, +20.9, +20.8 | **stage 1 empty** | 144 / 182 / 182 α-quadruples, 0 S-solutions [M] |
| 751 | Z₃₀ × F₂₅, (17,2), (17,14) | ×17 (order 4) on Z₃₀ with index 2 | +12.6 | **stage 1 empty** | 317,548 α-quadruples, 0 S-solutions [M] |
| 883 | Z₁₈ × F₄₉, (1,3) | GP, index 3 | +26.8 | **empty** | [P] no lag-0 composition (letters ±15, ±17, ±47, ±49 only) |
| 883 | Z₁₈ × F₄₉, (1,2) | Paley(49) | +4.0 | **empty** | [P] lag 0: 3264 is not 2208a + 2400b |
| 751 | Z₃₀ × F₂₅, (1,2) | Paley(25) | +25.8 | **empty** | [P] lag 0 (`nd_sym_s1est.py 30 25 2 --nosym`: no quadruple) |

**Proof for 947 ×23 [P].**
1. The slice is Z₂₂ × Paley(43), with p = 43 and n′ = 22. The letters are α + 42w with w ∈ {0, ±1}, i.e.
   {±1, ±41, ±43}.
2. Lag 0 needs Σy² = 3616 over 88 entries, so 1680a + 1848b = 3528. Hence a = b = 1: exactly one |41| and one |43|,
   i.e. exactly K = 2 cosets with w ≠ 0.
3. Lemma 4 of `report/paley_design_answer.md` gives, for every h ≠ 0, 1764·s(h) = −172 − a(h) − 42c(h), with
   |a| ≤ 88 and |42c| ≤ 84K = 168. So s(h) = 0.
4. Then Σ_k t_k² = K + Σ_{h≠0} s(h) = 2.
   - If both w ≠ 0 cosets lie in one block, that block has t ∈ {0, ±2}, so Σt_k² ∈ {0, 4}.
   - If they lie in different blocks, each such block has |row| ≥ 41 − 21 = 20 > 2.

   Either way there is a contradiction. ∎

---

## 2. 739 ×215 (Q1): stage 1 is exactly empty [P, M]

### 2.1 Reduction [P]

- 215 ≡ −1 (mod 18) and has order 5 mod 41. So M = ⟨−1 on Z₁₈⟩ × ⟨order 5 on Z₄₁⟩, and stage 1 is:
  - y_k symmetric on Z₁₈, with letters y = α + 10t (α = ±1, |t| ≤ 4);
  - Σ_k PAF(y_k) = (2792, −164, …);
  - rows (2, 0, 0, 0).
- **α-congruence.** For symmetric α and t the two cross sums agree, so

      PAF(y) = PAF(α) + 20·X + 100·PAF(t),   X(h) = Σ_g α(g) t(g + h).

  Hence Σ_k PAF(α_k)(h) ≡ −164 (mod 20) for h ≠ 0, and the rows of α are ≡ R_k (mod 10).
  - α_k ranges over the 1024 symmetric ±1 sequences.
  - Up to per-block shift by 9 and negation, with row-0 blocks sorted, **21 α-quadruples** survive.
- **S/D split.** Write Z₁₈ = Z₂ × Z₉ and put σ(b) = (y(b) + y(b+9))/2, δ(b) = (−1)^b (y(b) − y(b+9))/2. Both are
  symmetric on Z₉, and stage 1 is equivalent to
  - Σ_k PAF₉(σ_k) = (657, −82, −82, −82, −82), with σ-rows (1, 0, 0, 0);
  - Σ_k PAF₉(δ_k) = (739, 0, 0, 0, 0).
- Given α, σ involves only u = t(b) + t(b+9) and δ involves only w = t(b) − t(b+9), with u ≡ w (mod 2) and
  |u| + |w| ≤ 8.
  - The per-block u-lists have about 150 entries, so both joins are trivial.
- **Search** (`theory/nd_sym_stage1.py`, reimplemented generally in `theory/nd_s1_sd.py`): 21 α-quadruples →
  **12 S-solutions → 0 D-completions** (1 s). The general code, with different canonicalisation, gives 255 →
  102 → 0.

### 2.2 Validation [M]

- **Independent brute force** (`theory/nd_sym_bruteforce.py`): no α, no S/D; a direct 2+2 join of all
  configurations by energy class. It gives the same canonical solution counts on:
  - F = 6, q = 17: 774;
  - F = 6, q = 41: 0;
  - F = 10, q = 17: **1,750,702**;
  - F = 10, q = 41: 0.

  Log: `logs/nd/sym_stage1_crosscheck.log`.
- **Planted, full 739 size:** 20/20 realistic plants recovered (`logs/nd/s1_739_planted_real.log`). The plants are
  drawn uniformly at the real lag-0 energy and rows and use letters up to 21. The general solver recovers 10/10.
- **Model consistency:** the stage-1 Gaussian (`theory/nd_sym_s1est.py`) gives +11.2 bits for ordered quadruples.
  The symmetry group has order ≈ 2^14.2, so about 0.1 orbit is expected; emptiness is unsurprising.
- **CP-SAT** on the full model with all redundant constraints (`theory/nd_739_cpsat.py`, 3600 s, 2 workers):
  **UNKNOWN**. CP-SAT gives no independent certificate.

**Q2 (stage 2 for 739 ×215) is therefore moot.**

---

## 3. 971 ×469: stage 1 non-empty; stage 2 rich but large [M, E]

- **Stage 1 [M]** (`nd_s1_sd.py 970 469 97`): 8 α-quadruples, 66 S-solutions, **12 verified stage-1 solutions**
  (a few orbits under the Z₁₀ multipliers). Log: `logs/nd/s1_971_mu469_exact.json`.
- **Stage 2 per y [E]** (`theory/nd_stage2_est.py`, calibration in §9):
  - H₂ = 265.5 bits of freedom, C₂ = 228.7 bits (rank 47, Mahalanobis 0.36);
  - **+36.8 completions per y** for all 12 y; the target is inside the sampled lattice.
- **Level budget [E]** (`theory/nd_levels.py`):

  | level | bits |
  |---|---|
  | all | 228.7 |
  | S (Z₅ characters) | 130.9 |
  | D | 130.9 |
  | orders {1, 2} | 94.1 |
  | trivial Z₁₀ character | 54.5 |

- **Why it is out of reach.**
  - The S/D split is symmetric. The S-level carries 131 bits on Z₅, which has no intermediate subgroup, against
    about 132 bits of S-freedom. So about 2^1 S-solutions exist per y, and each costs a 4-list join of roughly
    2^{131/2}.
  - Every hierarchy tried leaves ≥ 2^60 per core [E].
  - Sub-slices with more symmetry are dead: 971 ×81, ×221, ×79 and ×101 are all empty.

---

## 4. A new family: generalised-Paley slices over F_{p²} with integer Gauss periods [P, M, E]

### 4.1 The decomposition [P, checked M]

**Setting.**
- G = Z_F × (F_q, +), with q = p². Z_F is fixed pointwise, and the index-s subgroup H ≤ F_q* acts by multiplication.
- Per block and g ∈ Z_F the slice has α(g) at z = 0 and c_i(g) on the class C_i = γ^i H.
- Stage 1 (the F_q-compression) is identical to the prime case.

**Uniform cyclotomy.**
- If −1 is a power of p mod s (here s | p + 1), the periods η_i(a) = Σ_{z ∈ C_i} ψ_a(z) equal A on one class and B
  on the rest, with B = (−1 − p)/s and A = B + p (Baumert–Mills–Ward [L]).
- Checked numerically (`theory/nd_fq.py check`):

  | q | s | periods |
  |---|---|---|
  | 49 | 8 | (6, −1⁷) |
  | 49 | 4 | (5, −2³) |
  | 25 | 6 | (4, −1⁵) |
  | 25 | 3 | (3, −2²) |
  | 9 | 4 | (2, −1³) |
  | 121 | 12, 6, 4 | integer |
  | 169 | 14, 7 | integer |

  For s ∤ p + 1 the periods are not integers (e.g. 49 with s = 3 or 6).

**Hence**, with σ = Σ_i c_i = (y − α)/f:

    x̂(g; ψ_a) = α(g) + B σ(g) + p c_{i*(a)}(g),

so **stage 2 ⟺ for every class i: Σ_k PAF_{Z_F}(w_{k,i}) = 4v·δ**, with w_{k,i} = α_k + Bσ_k + p c_{k,i}.
- The s class systems share **one** quadratic map Φ(c) = R (the same R for every class).
- They are coupled only through the per-position sums Σ_i c_i(g) = σ(g).
- Verified to 4·10⁻¹⁵ against full-group DFTs (`nd_fq.check_decomposition`, five (F, p, s) cases).

### 4.2 Which open orders have such a group

- v − 1 = F·p² with s | p + 1, s even, letters including ±1 (so that lag 0 can be met):
  - **883 = 18·7² + 1** (Z₁₈ × F₄₉, s = 4, 8): the closest analogue of the 523 slice;
  - **751 = 30·5² + 1** (Z₃₀ × F₂₅, s = 6);
  - 739, 919, 955 and 991 (F₉, s = 4, with F = 82, 102, 106, 110).
- Every other open order has no square factor of this kind [P, arithmetic].

### 4.3 Sweeps and stage splits [E]

`theory/nd_fq_sweep.py` runs all subgroups of Z_F* × F_q* with 16–90 orbits; `theory/nd_fq_split.py` gives the exact
lag-0 count × the lattice Gaussian. Logs: `logs/nd/fqsweep/`.

| order | slice (u on Z_F, e on F_q) | orbits | whole slice | stage 1 | exact status |
|---|---|---|---|---|---|
| 883 | (7,4) / (7,20) | 82 | +64.3 | +22.6 | **empty** [M] |
| 883 | (1,4) GP s = 4 | 90 | +52.9 | +27.0 (≈1 orbit) | stage 1 undecided (18 orbits per block; too many α-quadruples) |
| 883 | (1,8) GP s = 8 | 162 | +147.7 | **+76.5** | **stage 1 non-empty**: 6 verified y in 10 min (`q523_ymitm.py`, `logs/nd/y883_s8/`) |
| 883 | (17,4), (7,2), (7,10) | 82, 46, 46 | +37.4, +20.9, +20.8 | −0.4, +2.1, – | **empty** [M] |
| 751 | (7,4), (7,2), (7,14) | 84, 72, 72 | +53.4, +15.0, +15.2 | +41.7, +22.9, +22.9 | undecided (too many α-quadruples) |
| 751 | (17,8) | 90 | +24.0 | +42.8 | **stage 1 non-empty** [M]; stage 2 per y **−15** with Mahalanobis 9 (`nd_fq_split.stage2_estimate`) |
| 751 | (17,2), (17,14) | 70 | +12.6 | – | **empty** [M] |
| 751 | (11,2) / (19,2) | 80 / 78 | +39.3 / +36.9 | +24.8 / +21.8 | undecided (α lists too large) |
| 739 (F₉) | (37,4) / (23,0) | 90 | +66.5 / +62.1 | – | not screened |
| 955 (F₉) | (13,0) | 90 | +65.9 | – | not screened |

### 4.4 Stage 2 of 883 Z₁₈ × F₄₉ (s = 8): rich in estimate, no tractable solver [E]

On a real y (composition 16×|11|, 8×|13|, 48×|1|), with H₂ = 433.6 bits, the two estimates disagree:
- **Joint pattern model** (`nd_fq.levels_classes`): C₂ = 374.2 bits (rank 79), so **+59.3 completions per y**.
- **Independent class systems** (`nd_fq.class_bits`): X = 52.0 bits per class system, so |S_Φ| ≈ 2^20, the coupling
  costs 142.4 bits, and **+17.8 completions per y**.

Neither model yields a cheap solver:
- **S/D split:** h_D = 90.8 and C(D|S) = 160.1, so D_def ≈ 69 bits.
- **Coupling:** an 8-list problem with lists of 2^20 against at least 101 bits needs about 2^62 (extended k-tree).

**Calibration at v = 55** (Z₆ × F₉, s = 4; `theory/nd_fq_control.py`, `logs/nd/fq_control_v55.log`):
- **Exact count.** Six real stage-1 y were completed **exactly** by an exhaustive 2+2 join over all per-block
  pattern configurations, and every one has **0** completions.
  - The pattern model had predicted 2^11 per y (target at the centre of the sampled distribution); the independent
    model 2^−15.5.
  - A planted target is found by the same join (18 hits), so the join is sound.
- **Explanation: a new exact obstruction [P, M] (`nd_fq.dual_row_screen`).**
  - At a real character χ of Z_F (χ = 1, or χ = (−1)^g), each class needs
    Σ_k (χ̂(β_k) + p·T_{k,i})² = 4v, where T_{k,i} = χ̂(c_{k,i}) is a ±1 sum with F's parity and |T| ≤ F.
  - Summed over classes, Σ_i T_{k,i} = χ̂(σ_k) for every block.
  - At χ = 1 this is exactly the **other compression** (along Z_F onto F_q), so F_q slices are **doubly staged**,
    like 647 (§5).
  - At v = 55 the χ = 1 condition has 24 admissible per-class quadruples but **0** class assignments, for every y
    tried (40 y). That is the reason for the zeros. The Gaussian cannot see it, because the row values form a sparse
    discrete set.
  - For the six 883 s = 8 y the screen passes, with 1.3·10⁹ (χ = 1) and 2·10⁷ (χ = (−1)^g) assignments. So 883
    survives, but its stage-2 count stays uncertain between the two models.
- **Lesson.** For F_q slices, run `dual_row_screen` before any estimate. The real-character rows act as a discrete
  middle level, the analogue of Lemma 6 for Paley slices.

---

## 5. 647: the doubly-staged slices (Q3) — the best remaining route [M, E]

### 5.1 Structure [P]

- m = 646 = 2·17·19, and Z₂ is fixed by both slices.
- **×353:** 353 ≡ 13 (mod 17), order 4, index 4, f = 4; and 353 ≡ 11 (mod 19), order 3, index 6, f = 3. So
  M = M₁₇ × M₁₉, a product.
- **×183:** 183 ≡ 13 (mod 17), order 4; and 183 ≡ 12 (mod 19), order 6, index 3, f = 6. Here M ≅ Z₁₂ is diagonal,
  not a product. The fibres over b₁₇ ≠ 0 carry index-6 classes; those over b₁₇ = 0 carry index-3 classes; and a
  twist x(13b, 12z) = x(b, z) links positions.
- **Character levels** (ε, χ, ψ) on Z₂ × Z₁₇ × Z₁₉:
  - level 0: χ = ψ = 1 (rows, alternating sums);
  - stage 1A: χ = 1, ψ ≠ 1 (the Z₃₈ compression along Z₁₇);
  - stage 1B: χ ≠ 1, ψ = 1 (the Z₃₄ compression along Z₁₉);
  - cross: χ ≠ 1, ψ ≠ 1.

  1A and 1B meet only at level 0.

### 5.2 Exact stage 1 [M] (`theory/nd_batch_exact.py`, `logs/nd/exact_batch.jsonl`, solutions in `logs/nd/s1sols/`)

| slice | Z₃₄ compression (1B, 10 orbits per block) | Z₃₈ compression (1A) |
|---|---|---|
| ×353 | **non-empty** (35,907 α-quadruples, ≥ 200 solutions) | undecided by α-first (too many α-quadruples); Gaussian +51.6 |
| ×183 | **non-empty** (35,907 α-quadruples, ≥ 200 solutions) | **non-empty** (732,160 α-quadruples, ≥ 200 solutions) |

So neither 647 slice has an obstruction at the compressed level.

**Joint stage 1 of ×183 is non-empty [M]** (`logs/nd/s1_647_joint_z2.log`).
- The two compressions must agree on the Z₂ level: per block, the row and the alternating sum A_k, with
  Σ A_k² = 2588.
- Among the solutions found, the value with rows (0, 0, 0, 2) and |A| = (2, 30, 30, 28) is shared by 45 stage-1A and
  56 stage-1B solutions. So consistent pairs (y_A, y_B) exist.
- What remains is the cross level (χ ≠ 1, ψ ≠ 1): 126.6 − 40.6 ≈ 86 bits given both y's [E].

### 5.3 The cheap ordering: fix y_B, then run a 523-style S/D route [E]

Given y_B, the remainder of the core is all characters with ψ ≠ 1, which includes the stage-1A content. It is a
cyclotomic problem over Z₃₄ = Z₂ × Z₁₇ with Z₁₉-fibres. Level budgets from `theory/nd_levels.py`
(`logs/nd/levels_647.log`; h_D exact) on real y_B:

| slice, y_B | H₂ | C₂ (Mahalanobis) | S (ε = +1) | D | χ = 1 level (= 1A) | middle (ε = 1, χ = 1) | h_D | completions per y_B | D_def |
|---|---|---|---|---|---|---|---|---|---|
| ×183 #0 | 142.0 | 126.6 (1.3) | 72.6 | 73.8 | 40.6 | 22.8 | 27.3 | **+15.4** | **26.7** |
| ×183 #40 | 143.7 | 127.3 (1.5) | 72.9 | 74.3 | 41.6 | 22.9 | 28.2 | +16.3 | 26.3 |
| ×183 #160 | 145.2 | 126.4 (1.2) | 72.6 | 73.6 | 40.5 | 22.9 | 28.3 | **+18.8** | **25.4** |
| ×353 #0 | 162.8 | 142.3 (4.7) | 80.0 | 81.5 | 63.7 | 32.5 | 31.5 | +20.5 | 30.9 |
| ×353 #160 | 157.5 | 138.9 (4.1) | 80.0 | 77.9 | 59.4 | 32.6 | 30.0 | +18.7 | 28.9 |
| *523, y of record* | 179.7 | 146.0 | 83.3 | 83.2 | – | 31.8 (Z₃) | 36.4 | +34.1 | 26.1 |

**Reading.**
- For ×183 the stage 2 given y_B is a smaller copy of 523's stage 2:
  - fewer constraint bits (127 against 146) and a smaller S-level (73 against 83);
  - a smaller middle level (23 against 32);
  - the same D-deficit (≈ 26);
  - healthy Mahalanobis distances (1.2–1.5).
- **Projected cost of the 523 route here [E, not yet calibrated for this shape]:**
  - 2^{D_def} ≈ 2^26 S-solutions per core;
  - each costs 2^{F + 1 − n} pair visits, with F ≈ 72.8 − 22.9 ≈ 50 fine bits; at the 523 side size n ≈ 36, that is
    about 2^15;
  - plus a D-check of the 523 kind;
  - total about **2^41 pair visits + 2^26 D-checks, roughly 2^41–2^44**: tens of core-hours, below the 2^45 target.
- **Caveats.**
  - 523's realised S-step depended on how many parent pairs share a middle value (§3.3 of
    `report/paley_523_answer.md`); that has to be measured here.
  - The middle level is now the trivial Z₁₇-character (one parent value per block, a Z₁₉-function), not a Z₃ split.
- For ×353 the target sits 4–4.7 standard deviations out, so its estimates are less reliable, and D_def ≈ 29–31.

### 5.4 Exact first computation (647 ×183)

1. **Choose y_B.** Among the ≥ 200 exact Z₃₄ solutions (`logs/nd/s1sols/cyc_646_183_19.json`), take the one with the
   lowest D_def (#160: 25.4) by `nd_levels.py 646 183 19 FILE --idx i`.
2. **Port the 523 pipeline to a general cyclic slice given y.**
   - S level: the Z₂ pair sums on Z₁₇ with the Z₁₉-fibre patterns.
   - Parents: the χ = 1 (Z₁₇-trivial) sums per block. Their equations are the 22.9-bit middle level.
   - Fibres: the S-configurations per parent.
   - Join: a Wagner two-level join on the 50 fine bits.
   - D-check: the D options given S per coset pair (h_D ≈ 28 bits in all), checked by a 2+2 join on the D level.
   - The algebra needs no new derivation: the PAF signature of any configuration on Z₆₄₆ is computed directly by FFT
     (as in `nd_stage2_est.py`), so the levels are integer projections of lag sums.
3. **Positive control at a smaller order of the same shape.** Z₂ × Z₁₇ × Z_q′ with the same order-4 action on Z₁₇
   and a GP action on Z_q′:
   - v = 239 (q′ = 7, index 3, f = 2);
   - v = 443 (q′ = 13, index 4, f = 3).

   Both share the Z₃₄ stage-1B structure. Steps: stage 1B by `nd_s1_sd.py`, then the ported pipeline, then bordered
   GS completion and `verify_h2292.py`-style verification.
4. **Measure, as for 523:** S-solutions per middle value against the Gaussian, D-success per S-solution, ms per
   S-solution. Then run 647 on 4 cores.

---

## 6. 947 ×49 (Q4): a hierarchical stage 2 does not close the gap [E]

- **Status.**
  - The Z₄₃ compression (Z₂₂) is non-empty [M], with ≥ 200 solutions.
  - The Paley(11) compression (Z₈₆, the real stage 1) is undecided: 30 orbits per block is too many for α-first.
  - The gadget search (`paley_slicepipe.py guided1 947 11 86 49 32 11`) is in `logs/nd/guided1_947_K32.log`.
- **Stage split** (`logs/nd/split_947_mu49.log`, the estimator of the 573 design):
  - stage 1: +33 … +48 per gadget class;
  - stage 2: **−18.6 … −23.4 per y**, with 106–109 u-bits and C₂ ≈ 128–130 bits;
  - on 573 the same estimator gave −9.6 at K = 21, and 573 was realised at about −7 (1 core in 141 y's);
  - this session's calibrated estimator predicted 0.6 solutions on the 151 exactly-joined 573 y's, and 2 were found.
- **Hierarchy.** Z₈₆ = Z₂ × Z₄₃, and ×49 acts on Z₄₃ with order 3.
  - The S/D split on Z₂ gives an S-level on Z₄₃ with about 64 bits of constraint against about 54 bits of S-freedom.
    That is about 2^−10 S-solutions per y at about 2^27 per S-join.
  - The D level has a similar deficit.
  - Z₄₃ is prime, so there is no middle level.
  - The cubic classes of Z₄₃ under ×6 are Galois-conjugate, so they add no class keys beyond the Z₂ split (as at
    523, §1.3 of `report/paley_523_answer.md`).
- **Cost.**
  - 2^{19 … 23} y's per core, each at ≥ 2^27 (S-join) plus stage-1 generation: **≥ 2^46–2^50 even with free y's**.
  - With the guided gadget pipeline (≈ 10–20 core-seconds per y at 573) it becomes ≥ 2^65 operations.
- **Not a candidate.**

---

## 7. Screens over the remaining slices (Q5) [M, P]

**Cyclic sweeps for 750–1000** (`logs/msweep34/`, regenerated this session; it is gitignored, so the positive rows
are listed here).
- Every slice with log₂E > 0 not dead before this session is now screened exactly:

| v | slice (orbits, log₂E) | result |
|---|---|---|
| 789 | ×29 (20, +9.3) | dead before (lag 0) |
| 823 | ×73 (54, +52.4) dead before; **×59** (36, +7.9) | **empty** [M] |
| 859 | ×289 (90, +27.4) | Z₇₈ and Z₆₆ compressions too large for α-first; Gaussians +74.4 and +87.1 [E]; **open** |
| 907 | ×85, ×91 dead before; **×155** (44, +13.6) | **empty** [M] |
| 917 | ×17 dead before; ×9 (20, +8.7) | dead before (GP s ≤ 4) |
| 933 | **×37** (36, +30.6) | **empty** [P, lag 0] |
| 947 | ×49 alive; ×45 dead before; **×21, ×47, ×23** | **empty** [M, M, P] |
| 955 | ×289 | dead before |
| 971 | ×81, ×221, ×79, ×101 dead before; **×469** | **alive** (§3) |
| 751, 883, 919, 991 | no cyclic slice > 0 in 16–90 orbits | see §4 for noncyclic |
| 839, 863 | sweeps still running at the time of writing | – |

**Below 750.**
- **537:** ×5, ×137, ×17 (+5.3, +4.9, +1.2) compress to Z₈ = Z₂·Z₄, where the S/D split needs another form
  (negacyclic D level). **Not screened.**
- Everything else was settled before.

**Noncyclic groups** (§4): 883, 751 and the F₉ orders.

**b > 1 borders (750–1000):** not studied this session.

---

## 8. Ranked list (D)

| rank | v (order) | slice | exact status | estimate | next computation |
|---|---|---|---|---|---|
| 1 | **647 (2588)** | **×183** | both compressions non-empty [M] | +15 … +19 completions per y_B; 523-type route ≈ 2^41–2^44 per core [E] | §5.4: port the S/D pipeline; control at v = 239 / 443 |
| 2 | 647 | ×353 | 1B non-empty [M]; 1A undecided | +18 … +20 per y_B, but Mahalanobis 4–4.7; D_def ≈ 29–31 | as rank 1, after ×183 |
| 3 | 971 (3884) | ×469 | stage 1: 12 y [M] | +36.8 per y; S-level critical, ≥ 2^60 [E] | none cheap; a middle level for Z₅ would be needed |
| 4 | 883 (3532) | Z₁₈ × F₄₉, s = 8 | stage 1 non-empty (6 y) [M] | stage 2 between +17.8 and +59 per y (models disagree; v = 55 control favours the low end) | exact-count calibration at a larger control; class-system solver |
| 5 | 859 (3436) | ×289 | open (compressions too large for α-first) | compressions +74 / +87; whole +27.4 | CP-SAT or `q523_ymitm`-style search for both compressions, then a level budget |
| 6 | 751 (3004) | Z₃₀ × F₂₅ (7,4), (17,8), (11,2), (19,2) | exact batch pending | stage 1 +22 … +43 | `nd_batch_exact.py` results; stage 2 per real y |
| 7 | 739, 955 (F₉) | (37,4), (23,0); (13,0) | not screened | whole +62 … +66 | `nd_fq_split.py` + exact batch |
| 8 | 947 (3788) | ×49 | Z₂₂ compression non-empty; Z₈₆ open | ≥ 2^65 [E] | not recommended |
| – | 537 | ×5, ×137, ×17 | not screened (F = 8) | ≤ +5 | negacyclic S/D variant |

**Dead this session:** 739 ×215; 947 ×21, ×47, ×23; 823 ×59; 907 ×155; 933 ×37; 883 F₄₉ (7,4), (7,20), (1,3),
(1,2); 751 F₂₅ (1,2).

---

## 9. Tools and calibrations

### 9.1 Exact stage-1 solvers

- `theory/nd_sym_stage1.py` (symmetric, the 739 case) and `theory/nd_s1_sd.py` (any cyclic slice compressed to
  Z_{2n}, n odd; or an explicit noncyclic structure with `--explicit`).
- `theory/nd_batch_exact.py` runs them restart-safely, with a cap of 200 solutions.
- Controls:
  - `theory/nd_sym_bruteforce.py` and `nd_sym_compare.py`: counts agree at F = 6 and 10;
  - `--plant-real`: 20/20 at 739; 10/10 at 947 ×21; 10/10 with the general code at 739.

### 9.2 Estimators

- `theory/nd_stage2_est.py` (stage 2 given y, any cyclic slice):
  - 523: H₂ 179.7 / C₂ 146.3 (record 180.3 / 146.2);
  - 573: 151 exactly-joined y's, 0.6 predicted against 2 found (`logs/nd/calib573_stage2.log`).
- `theory/nd_levels.py` (levels by Z_F-character order, exact h_D): on 523 it gives S 83.65, D 82.86, Z₃ 32.49,
  Z₁ 11.54, h_D 36.43 (record 83.3, 83.2, 31.8, 11.4, 36.44).
- `theory/nd_sym_s1est.py`, `nd_compress_est.py`, `nd_fq_split.py`: compressed stage-1 Gaussian with exact lag-0
  counts. The 523 control gives +48.9.

### 9.3 F_q tools

- `theory/nd_fq.py`: field, period check, decomposition check, class systems, the two stage-2 models.
- `theory/nd_fq_sweep.py`, `nd_fq_split.py`: sweeps and splits.
- `theory/nd_fq_control.py`: exact v = 55 control (stage 1, exhaustive stage 2, estimators, bordered GS).

### 9.4 Calibration warnings (what this session adds)

- The whole-slice Gaussian was off by a large amount again:
  - 739 ×215 (+61.9), 947 ×21 (+28.3) and 883 (7,4) (+64.3) are empty;
  - 883 GP s = 8 gives +147.7 whole against +94–136 from the stage split.
- For F_q slices the joint pattern Gaussian can overcount stage 2 by more than 11 bits (v = 55).
- For cyclic slices the same estimator is within 1–2 bits (523, 573).

---

## 10. Reproduce

```
python3 theory/nd_sym_stage1.py 18 41 8                         # 739 x215: 21 alpha-quadruples, 12 S-solutions, 0 (1 s)
python3 theory/nd_sym_stage1.py 18 41 8 --plant-real 20         # planted 20/20
python3 theory/nd_sym_bruteforce.py 10 17 8; python3 theory/nd_sym_compare.py 10 17 8   # 1,750,702 both ways
python3 theory/nd_s1_sd.py 970 469 97 --out FILE                # 971 x469: 12 stage-1 solutions
python3 theory/nd_s1_sd.py 946 21 43; python3 theory/nd_s1_sd.py 946 21 43 --plant-real 10   # 947 x21 empty; 10/10
python3 theory/nd_batch_exact.py logs/nd/exact_jobs.txt OUT.jsonl                       # the batch of section 1/5/7
python3 theory/nd_stage2_est.py 970 469 97 logs/nd/s1_971_mu469_exact.json             # +36.8 per y
python3 theory/nd_levels.py 646 183 19 logs/nd/s1sols/cyc_646_183_19.json --idx 160 --levels "1,17;2,34;1,2;1"
python3 theory/nd_fq.py check 7 8; python3 -c "import sys; sys.path.insert(0,'theory'); import nd_fq; nd_fq.check_decomposition(18,7,8)"
python3 theory/nd_fq_sweep.py 18 7 --N 3000; python3 theory/nd_fq_split.py 18 7 "1,8"
python3 theory/q523_ymitm.py 18 49 8 --comp 11:16,13:8,1:48 --N 20 --secs 600 --out DIR   # 883 stage-1 y's
python3 theory/nd_fq_control.py 6 3 4 --maxy 6                   # v = 55 exact control
```
