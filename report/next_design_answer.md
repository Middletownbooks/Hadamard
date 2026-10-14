# Next design: 647 ×353 survives the joint stage 1 (cross level ≈ 2^−5 per joint solution [E]); 647 ×183, 739 ×215 and 537 exactly empty; a new F_{p²} family

Date: 2026-10-13. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_next_design.md`.
Scripts: `theory/nd_*.py`. Logs: `logs/nd/`.

**Tags.** [M] exact computation (script named); [P] proved here (proof in the text); [E] Gaussian estimate (the
estimator is named, and so is its calibration); [L] literature recollection, unchecked (arXiv and OEIS are blocked).

---

## 0. Bottom line

1. **(A) No new Hadamard matrix this session.**

2. **(B) Not yet reached, but 647 ×353 is now the live candidate (§5).**
   - **Joint stage 1 is non-empty [M].** Stage 1B (the Z₃₄ compression) has exactly **495 orbits** of solutions
     under the slice's symmetry group. For the first orbit (y_B #0) there are **20,928 joint solutions (y_A, y_B)**:
     a labelling consistent with both compressions that solves both stage-1 systems.
   - **Cross level [E, uncalibrated]:** given a joint solution, H = 79.4 bits of freedom against C = 84.3 bits
     (rank 23, Mahalanobis 0.48). So **log₂E ≈ −4.9 per joint solution**, about 2^9.5 expected cores from y_B #0
     alone if the Gaussian holds. It is stable over 3 sampled joint solutions (`logs/nd/joint/cross_646_353_k0.log`).
   - **Exact cross solver:** written (`theory/nd_cross.py`: per-layer labelling lists, S-join over the distinct σ,
     then the full join), but **not yet validated or run**. It is the next computation (§11).
   - **Caveat.** No same-shape positive control exists below 647 (239 and 443 are empty). The Gaussian has failed on
     discrete levels before (§9.4), but the cross level has no further compression, so the classic failure mode does
     not apply. Calibrate first by an exact S-level count on one joint solution (§11).
   - **Routes out of reach [E]:**
     - 971 ×469: ≥ 2^60. Its dual compression is non-empty (§3), so it survives the joint test.
     - 947 ×49: ≥ 2^65 (§6).
     - 883 on Z₁₈ × F₄₉: no route (§4).

3. **(C) Exact obstructions [P, M]** (§1, §2, §7).
   - **739 ×215 is exactly empty**: its stage 1 has no solution. This kills the most valuable candidate below 3000.
   - **647 ×183 is exactly empty**: all 4 orbits of its stage 1B fail the joint test (§5).
   - **537 is finished**: ×5, ×137 and ×17 die at lag 0, so every positive cyclic slice of 537 is empty (§7).
   - Also exactly empty:
     - 839 ×13 (Gaussian +101.7), 839 ×7, 863 ×3;
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

5. **Why the Gaussians failed so often [M].** Every slice whose multiplier acts on two prime factors has two
   compressions. Given one compression, the other depends on only a handful of discrete values per block.
   - The lattice Gaussian treats that level as continuous and can be off by more than 10 bits.
   - The F_q control at v = 55 (pattern model 2^11, exact 0) and the 647 joint test (implied ~2^16, exact 0) show it.
   - The 523 slice is immune: its dual compression is a sum of 18 signs per class, close to continuous, which is why its
     calibrations held.
   - **Rule:** for such slices, decide the joint stage 1 exactly before costing anything
     (`theory/nd_joint_exhaust.py`). It is cheap once the y_B are reduced by the symmetry group: 647 ×183 took 6 s
     in total.

6. **(D) Ranked candidates with their next computations: §8.**

---

## 1. Exact obstructions found this session

| v | slice | structure | Gaussian [E] | exact result | how |
|---|---|---|---|---|---|
| **739** | **×215** | −1 on Z₉ (y symmetric on Z₁₈); index 8 on Z₄₁ (f = 5) | +61.9 whole, +11.2 stage 1 | **stage 1 empty** | §2: α-congruence mod 20 leaves 21 α-quadruples (canonical), S/D split gives 12 S-solutions, 0 completions [M] |
| **647** | **×183** | Z₂ fixed; order 4 on Z₁₇ × order 6 on Z₁₉ (diagonal) | +24.6 whole; +15 … +19 per y_B | **empty** | stage 1B has 4 orbits (962 canonical α-quadruples, 18 raw solutions, 6 s); each fails the joint test (0 S-solutions) [M] (`nd_joint_exhaust.py 646 183 17 19`) |
| 537 | ×5, ×137, ×17 | Z₈ × index 3 / 3 / 2 on Z₆₇ | +5.3, +4.9, +1.2 | **empty** | lag 0 [M] (`nd_lag0_screen.py`): minimum energies 7072 and 14112 > 1880; for ×17 (letters ±1, ±65, ±67) 1880 is unreachable |
| 947 | ×21 | −1 on Z₁₁ (y symmetric on Z₂₂); order 7 on Z₄₃ (f = 7) | +28.3 whole, −8.1 stage 1 | **empty** | `nd_s1_sd.py 946 21 43`: 40 α-quadruples, 0 S-solutions [M]; planted 10/10 |
| 947 | ×47 | Z₂₂ compression: ×5 on Z₁₁ (QR), order 7 on Z₄₃ | +2.0 | **empty** | 4 α-quadruples, 32 S-solutions, 0 completions [M] |
| 947 | ×23 | Z₂₂ × Paley(43) | +1.7 | **empty** | [P] below |
| 823 | ×59 | Z₆ (−1 on Z₃) × index 8 on Z₁₃₇ (f = 17) | +7.9 | **empty** | no α-quadruple mod 68 [M] |
| 907 | ×155 | Z₆ (−1 on Z₃) × index 10 on Z₁₅₁ (f = 15) | +13.6 | **empty** | no α-quadruple mod 60 [M] |
| 933 | ×37 | Z₄ × index 8 on Z₂₃₃ (f = 29) | +30.6 | **empty** | [P] lag 0: letters ±1 + 58t; 57² − 1 = 3248 > 2800 − 16, so only ±1 fits, but Σy² must be 2800 over 16 entries |
| 839 | ×13, ×7 | Z₂ × index 38 / 22 on Z₄₁₉ (f = 11 / 19) | **+101.7**, +58.0 | **empty** | [P, M] below |
| 863 | ×3 | Z₂ × index 10 on Z₄₃₁ (f = 43) | +26.7 | **empty** | [P, M] below |
| 883 | Z₁₈ × F₄₉, (7,4) and (7,20) | ×7 on Z₁₈ (order 3) × index 4 on F₄₉* | +64.3 whole, +22.6 stage 1 | **stage 1 empty** | 565,348 α-quadruples, 863,136 S-solutions, 0 completions (1344 s) [M] |
| 883 | Z₁₈ × F₄₉, (17,4), (7,2), (7,10) | −1 or ×7 on Z₁₈ with index 4 or 2 on F₄₉* | +37.4, +20.9, +20.8 | **stage 1 empty** | 144 / 182 / 182 α-quadruples, 0 S-solutions [M] |
| 751 | Z₃₀ × F₂₅, (17,2), (17,14) | ×17 (order 4) on Z₃₀ with index 2 | +12.6 | **stage 1 empty** | 317,548 α-quadruples, 0 S-solutions [M] |
| 883 | Z₁₈ × F₄₉, (1,3) | GP, index 3 | +26.8 | **empty** | [P] no lag-0 composition (letters ±15, ±17, ±47, ±49 only) |
| 883 | Z₁₈ × F₄₉, (1,2) | Paley(49) | +4.0 | **empty** | [P] lag 0: 3264 is not 2208a + 2400b |
| 751 | Z₃₀ × F₂₅, (1,2) | Paley(25) | +25.8 | **empty** | [P] lag 0 (`nd_sym_s1est.py 30 25 2 --nosym`: no quadruple) |

**Proof for the Z₂ × Z_q slices of 839 and 863 [P]** (also `nd_s1_sd.py`: 1 α-quadruple, 1 S-solution, 0 completions).
1. With F = 2 the S/D split is trivial. Stage 1 reduces to (y₀ − 1)² + y₁² + y₂² + y₃² = v, where y_k = y_k(0),
   block 0 has row 2, every y ≡ ±1 (mod 2f), and y₀ ≡ 1 (mod 2f) (because 2 − y₀ must also be a letter).
2. Both v are ≡ 7 (mod 8), so y₀ = 1 is impossible: v would be a sum of three squares, which Legendre's theorem
   forbids.
3. Otherwise (y₀ − 1)² ≥ 4f².
   - At 863 (f = 43) and at 839 ×7 (f = 19), 4f² > v.
   - At 839 ×13 (f = 11) the only option is (y₀ − 1)² = 484, which leaves 355. But 355 is not a sum of three squares
     from {1, 441, 529, …}.
4. So no stage-1 solution exists. ∎ The +101.7 Gaussian was meaningless.

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
- **Doubly-staged caveat (§9.4).** The 54.5-bit trivial-character level is the *dual* compression onto Z₉₇.
  - Per block it has X(0) fixed by the rows and eight class values X_j = 2Y_j with Y_j ∈ [−5, 5] and ΣY_j = 0.
  - That is discrete, so the +36.8 may be far too high.
- **The joint test does not kill it [M].**
  - All 12 stage-1 y lie in **one orbit** of the symmetry group (`nd_joint_exhaust.Sym`).
  - Given y_B, the dual compression of each block is a transportation problem: column sums of 2 c-rows and 4 s-rows
    with fixed row sums. It excludes almost nothing (for block 0, only top-4 sums of Y above 18), so the joint test
    is essentially the dual compression itself.
  - **Exhaustive dual join** (`theory/nd_971_ejoin.py`): symmetry-reduced (class shift Z₈, negations, block order),
    energy-split numba hash join (`theory/nd_ejoin.py`), about 2^35.4 operations.
    - Validated: self-test 30/30 against brute force; plants found at E = 32 (3/3) and at E = 48 with forced
      partitioned passes (2/2).
    - The run was stopped at hand-off with energy classes EL ≤ 40 done. It had found **739 solutions, all verified
      by the full cyclic PAF on Z₉₇** (`logs/nd/dual971/`). Some are degenerate (Y₀ = Y₁ = 0).
  - **So 971 ×469 survives the joint stage 1.** Its cost stays ≥ 2^60 [E] (cross level about 174 bits of
    constraint), and it is not a candidate.

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

## 5. 647: the doubly-staged slices (Q3) [M, E]

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

### 5.2 Exact joint stage 1 [M] (`theory/nd_joint_exhaust.py`, logs in `logs/nd/joint/`)

**Method (all exact).**
1. **Symmetry group [P].** It is generated by:
   - multipliers u ∈ Z₆₄₆*, acting on Z₃₄ through Z₃₄*/⟨13⟩ (representatives 1, 3, 7, 9);
   - per-block shifts by the M-fixed element 323 (17 on Z₃₄, 19 on Z₃₈);
   - negation of row-0 blocks;
   - permutations of row-0 blocks.

   Each acts on the whole joint problem, so the joint test depends only on the orbit of y_B.
2. **Canonical stage-1B α-quadruples.** The 35,907 α-quadruples reduce to **962** canonical ones. Every stage-1B
   solution over them is enumerated by the S/D pipeline of `nd_s1_sd.py`, then canonicalised.
3. **Joint test per y_B orbit.**
   - M fixes Z₂, so labels, y_B and y_A split by layer. Per block and layer, every achievable y_A layer comes from a
     vectorised DP.
   - Exact energy filters [P]: Σ over the 8 layers of ‖L‖² = 2(TS₀ + TD₀), plus the per-block S- and D-energy
     bounds.
   - Then the S-join over distinct σ (`theory/nd_ejoin.py`), and a full S+D join over the pairs whose σ occurs in
     some S-solution.
4. **Validation.**
   - Random plants: 10/10 at ×183, 6/6 at 443 ×225.
   - Plants on real y_B at real energy slack: 6/6 at ×353, 4/4 at ×183.
   - `nd_ejoin.py --selftest`: 30/30 against brute force.
   - It reproduces the earlier slow results (y_B #0 and #36 of ×183: 0).

**Results.**

| slice | canonical α-quads | raw stage-1B solutions | y_B orbits | joint test |
|---|---|---|---|---|
| ×183 | 962 | 18 | **4** | **0 for all 4 (0 S-solutions each): slice exactly empty [M]** |
| ×353 | 962 | 1098 | **495** | **y_B #0: 20,928 joint solutions** (928,080 S-solutions, 140 s); orbits #1–#494 not yet run (the run stopped on #1, probably out of memory; see §11) |

The 200 ×183 solutions saved earlier fall into the same 4 orbits. The two y_B that the old code could not decide
(too large) are the other two orbits, and both give 0.

### 5.3 Cross level of 647 ×353 [M lists, E counts] (`theory/nd_cross.py`)

- **Structure [P].** Given a joint solution, per block and layer the labels form a 5 × 7 ±1 matrix Λ: rows indexed
  by {0, A₁..A₄} (the Z₁₇-classes), columns by {0, B₁..B₆} (the Z₁₉-classes). Its weighted row sums are y_B and its
  weighted column sums are y_A. So x₀₀, a and c are fixed, and the 4 × 6 cell block has fixed margins.
- **Per-layer lists, exact.** For one joint solution the layer list sizes are (1200, 584), (740, 1110), (1860, 1110),
  (710, 936). σ = (layer₀ + layer₁)/2 on Z₃₂₃ has about 2^17.7 distinct values per block (212,180 and 227,204 for
  blocks 0 and 1).
- **Remaining equations.** The S/D split of Z₆₄₆ = Z₂ × Z₃₂₃ gives
  Σ_k PAF₃₂₃(σ_k) = (645, −2, …) and Σ_k PAF₃₂₃(δ_k) = (647, 0, …) at 20 lag classes. Of these, 12 + 12 are new
  (cross) constraints, and ‖σ‖² + ‖δ‖² = 323 per block links the two lag-0 equations.
- **Gaussian [E]** (3 joint solutions, N = 3000 per block):

  | quantity | bits |
  |---|---|
  | H | 79.3–79.4 |
  | C (S and D together) | 84.3 |
  | C_S | 49.8 |
  | C_D | 49.4 |

  So **log₂E = −4.9 … −5.0 per joint solution**, with the target at Mahalanobis 0.48.
- **Reading [E].**
  - y_B #0 alone has 20,928 joint solutions. Up to the stabiliser of y_B (multipliers acting only on Z₁₉; at most 6)
    that is about 3,500 inequivalent ones, so about 2^7 inequivalent cores are expected if the Gaussian holds.
  - Exact cost per joint solution: an S-join over four σ lists of about 2^17.7 (about 2^21 S-solutions expected),
    then a cheap D-check. A core needs about 2^5 joint solutions.
  - The S-join work after the energy split is **not yet measured**; the 2+2 bound is about 2^35. Rough cost per core:
    about 2^40 [E], **inside the (B) budget if confirmed.**
- **Calibration status.**
  - No same-shape positive control exists (239 and 443 are exactly empty).
  - The cross level has no further compression (every subgroup compression of Z₂ × Z₁₇ × Z₁₉ is already used by 1A
    and 1B), so the discrete-level failure of §9.4 does not apply directly.
  - The available calibration is an exact S-level count on one joint solution against 2^(H_S − C_S) (§11).

---

## 6. 947 ×49 (Q4): a hierarchical stage 2 does not close the gap [E]

- **Status.**
  - The Z₄₃ compression (Z₂₂) is non-empty [M], with ≥ 200 solutions.
  - The Paley(11) compression (Z₈₆, the real stage 1) is undecided: 30 orbits per block is too many for α-first.
  - The gadget search (`paley_slicepipe.py guided1 947 11 86 49 32 11 60 0 14 3300`, the tool that produced 573's
    y's) ran for its full hour without recording a stage-1 y (`logs/nd/guided1_947_K32.log`). At 573 it took 10–20
    core-seconds per y.
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
| 839 | **×13** (78, +101.7), **×7** (46, +58.0) | **empty** [P, M] |
| 863 | **×3** (22, +26.7) | **empty** [P, M] |

**Below 750.**
- **537:** ×5, ×137 and ×17 (+5.3, +4.9, +1.2) are **empty at lag 0** [M] (`theory/nd_lag0_screen.py`,
  `logs/nd/lag0_537.log`).
  - ×137 is GP with s = 3 and ×5 has M_b = M on the even positions: σ is odd, so |y| ≥ 21 on 32 and 16 entries
    respectively. That gives minimum energies 14112 and 7072, against 1880.
  - ×17 (s = 2) has letters ±1, ±65, ±67, and 1880 is not a sum of 32 of their squares.
  - With ×9, ×59 and ×15 dead before, **every positive cyclic slice of 537 is empty**.
- Everything else was settled before.

**Lag-0 sweep of every positive cyclic slice and every Hall factor** (`logs/nd/lag0_all_positive.jsonl`).
- It reproduces every lag-0 death in §1 and in `report/slice_screens_20261013.md`.
- 523 aside, the slices still alive after all exact screens are: **647 ×353** (§5), 971 ×469 (§3), 859 ×289
  and 947 ×49.

**Noncyclic groups** (§4): 883, 751 and the F₉ orders.

**b > 1 borders (750–1000):** not studied this session.

---

## 8. Ranked list (D)

| rank | v (order) | slice | exact status | estimate | next computation |
|---|---|---|---|---|---|
| 1 | **647 (2588)** | **×353** | stage 1B: 495 orbits [M]; **joint stage 1 non-empty** (y_B #0: 20,928) [M] | cross level log₂E ≈ −4.9 per joint solution (H 79.4, C 84.3, Mahalanobis 0.48) [E, uncalibrated]; ≈ 2^40 per core if confirmed | §11: validate `nd_cross.solve_cross` on plants; measure the S-join; exact S-count calibration; then run over the joint solutions of y_B #0 |
| 2 | 971 (3884) | ×469 | stage 1: 12 y, **one orbit** [M]; dual compression non-empty (739 verified in EL ≤ 40) [M] | +36.8 per y; cross level ≈ 174 bits; ≥ 2^60 [E] | none cheap |
| 3 | 883 (3532) | Z₁₈ × F₄₉, s = 8 | stage 1 non-empty (6 y) [M]; dual-row screen passes | stage 2 between +17.8 and +59 per y (models disagree) | class systems = line compressions (§4.1); the coupling is an 8-sum over one set S_Φ (≈ 2^20) |
| 4 | 859 (3436) | ×289 | open (compressions too large for α-first) | compressions +74 / +87; whole +27.4 | doubly staged (Z₁₁ and Z₁₃): `nd_joint_exhaust.py`-style joint test after a non-α stage 1B |
| 5 | 751 (3004) | Z₃₀ × F₂₅ (7,4), (11,2), (19,2) | undecided (α lists too large) | stage 1 +22 … +42 | exact stage 1 by another method |
| 6 | 739, 955 (F₉) | (37,4), (23,0); (13,0) | not screened | whole +62 … +66 | `nd_fq_split.py` + exact batch |
| 7 | 947 (3788) | ×49 | Z₂₂ compression non-empty; Z₈₆ open | ≥ 2^65 [E] | not recommended |

**Dead this session:**
- cyclic: 739 ×215; **647 ×183**; **537 ×5, ×137, ×17**; 839 ×13, ×7; 863 ×3; 947 ×21, ×47, ×23; 823 ×59;
  907 ×155; 933 ×37;
- F₄₉ slices of 883: (7,4), (7,20), (17,4), (7,2), (7,10), (1,3), (1,2);
- F₂₅ slices of 751: (17,2), (17,14), (1,2); (17,8) has stage 2 at −15 per y (Mahalanobis 9), so it is effectively
  dead [E].

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

### 9.3a Doubly-staged and join tools (added at the end of the session)

- `theory/nd_ejoin.py`: exact 4-list join, energy-split, numba open-addressing hash tables, partitioned passes,
  growing match buffer. `--selftest`: 30/30 against brute force.
- `theory/nd_joint_exhaust.py`: symmetry group, canonical α-quadruples, restart-safe stage-1B enumeration, and the
  joint test (per-layer DP, energy filters, S-join, full join). It has `--plant N` and `JointTest.planted_real`.
- `theory/nd_cross.py`: cross level given (y_A, y_B).
  - `Cross.lists`: exact per-layer labellings, by a row-ordered DP with per-step encoding.
  - `Cross.estimate`: lattice Gaussian of the 24 cross constraints.
  - `solve_cross` / `planted_cross`: exact S-then-full join on orbit representatives. **Written, not yet run.**
- `theory/nd_971_ejoin.py`: exhaustive dual compression of 971 ×469 (symmetry-reduced, restart-safe per energy
  class), with `--plant`.
- `theory/nd_lag0_screen.py`: exact lag-0 attainability of any compressed stage 1.
- Older and slower: `theory/nd_doubly.py` and `nd_647_joint.py`. `Doubly.cross` is too slow for ×353 (Python dicts
  over 24-dimensional states).

### 9.4 Calibration warnings (what this session adds)

- The whole-slice Gaussian was off by a large amount again:
  - 739 ×215 (+61.9), 947 ×21 (+28.3) and 883 (7,4) (+64.3) are empty;
  - 883 GP s = 8 gives +147.7 whole against +94–136 from the stage split.
- For F_q slices the joint pattern Gaussian can overcount stage 2 by more than 11 bits (v = 55).
- For cyclic slices the same estimator is within 1–2 bits (523, 573).
- **Doubly-staged slices.** Given one compression, the Gaussian of the other is meaningless (discrete values per
  block). The exact joint test decides it in seconds once y_B is reduced by symmetry:
  - 647 ×183: Gaussian about 2^16 joint solutions per y_B, exact 0;
  - 647 ×353: exact 20,928 for y_B #0.
- The 647 ×353 cross-level Gaussian (−4.9 per joint solution) has **no calibration yet**.

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
python3 theory/nd_ejoin.py --selftest                            # 30/30 against brute force
python3 theory/nd_joint_exhaust.py 646 183 17 19 --dir logs/nd/joint --workers 4   # 647 x183: 4 y_B orbits, joint 0 (6 s)
python3 theory/nd_joint_exhaust.py 646 183 17 19 --plant 10      # random plants 10/10 (also 442 225 17 13 --plant 6)
python3 theory/nd_joint_exhaust.py 646 353 17 19 --dir logs/nd/joint --workers 1   # 647 x353: 495 orbits; y_B #0 joint 20,928
python3 theory/nd_cross.py 646 353 17 19 logs/nd/joint/joint_646_353_17_19.jsonl --k 0 --n 3 --N 3000   # cross Gaussian -4.9
python3 theory/nd_971_ejoin.py --plant 3 --E 32; python3 theory/nd_971_ejoin.py --dir logs/nd/dual971 --workers 3   # restart-safe
python3 theory/nd_lag0_screen.py 536 5 67 536 137 67 536 17 67   # 537: all three empty at lag 0
```

---

## 11. Hand-off (state at the end of this session)

**Nothing is running.** Every job was stopped by PID. All runners are restart-safe:
- `logs/nd/joint/` holds one JSONL line per α-quadruple (`enum_*`) and per y_B orbit (`joint_*`);
- `logs/nd/dual971/` holds one line per energy class.

The 523 production files (`theory/q523_*`, `logs/q523/`) were not touched.

**647 ×353, in order:**
1. **Validate the exact cross solver.** Run `nd_cross.planted_cross(X, yb, ya, seed)` on a joint solution of y_B #0.
   - Load it as in `nd_cross.py __main__`: `X = Cross(646, 353, 17, 19)`, `ys = enumerate_yB(...)`,
     `rec, sols = X.J.test(ys[0])`, `ya = X.J.ya_full(sols[j]['layers'])`.
   - It must return the plant. `solve_cross` works on orbit representatives (35 per layer) to avoid the 5 GB
     arrays that killed the first attempt.
2. **Measure and calibrate the S-level.**
   - Run `solve_cross(X, yb, ya, s_only=True, log=print)` on one real joint solution. Record the S-join time, the
     work, and the exact S-solution count.
   - Compare that count with 2^(H_S − C_S), where H_S = Σ_k log₂(distinct σ_k) (about 4 × 17.7) and C_S ≈ 49.8 from
     `Cross.estimate`. This is the only available calibration of the cross Gaussian.
3. **If the S-join costs at most about an hour, run the full `solve_cross` over the 20,928 joint solutions of
   y_B #0.**
   - Deduplicate them first under the stabiliser of y_B: multipliers u ≡ ±(M mod 34) acting on Z₁₉, at most 6.
   - About 2^5 joint solutions per core are expected [E].
   - Any core → `nd_doubly.bordered_gs_cyclic` → write an independent `scripts/verify_h2588.py` in the pattern of
     `scripts/verify_h2292.py` → commit at once.
4. **Resume the joint enumeration over orbits #1–#494** (`nd_joint_exhaust.py 646 353 17 19 --dir logs/nd/joint
   --workers 1`).
   - The previous run died silently on #1, probably out of memory: its S-join has more than 10⁶ matches, and the
     buffer now grows up to 2^24.
   - Add a cap on S-solutions or chunk the full join before rerunning. This is needed only if y_B #0 is exhausted
     without a core.

**971 ×469.** Resume the dual run only for completeness (`nd_971_ejoin.py --dir logs/nd/dual971 --workers 3`;
EL ≤ 40 done). It cannot kill the slice: the dual compression is already non-empty.

**Report.** §0, §1, §3, §5, §7, §8 and §9 reflect everything above. The level budgets given y_B (the old §5.3 table)
were removed as superseded. The old ×183 numbers are in git history (commit ac5997d).
