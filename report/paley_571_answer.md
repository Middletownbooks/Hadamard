# v = 571 in the Paley + multiplier slices: the p = 3 re-split, exact compression arithmetic, a hierarchical solver, and the Paley(19) slice

Date: 2026-10-10. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_571_design.md`.
Scripts: `theory/p571_*.py`. Logs: `theory/logs/p571/`, `theory/logs/p571_*.jsonl`.

**Tags.**
- **[M]** verified by exact computation here; the script is named.
- **[P]** proved here; the proof is in the text.
- **[E]** estimate: Gaussian unless stated; §7 uses exact per-lag marginals instead.
  - The Gaussian is calibrated against exact counts (§5.3). At p = 3 it overcounts the stage-2 yield per y by
    2^1.5–2^3.5 (v = 79), and the 571 S-level count by about 2^5.
  - Whole slices it predicts rich can be exactly empty (§4, §7).

---

## 0. Bottom line

1. **(A) No H(2284) was constructed.**

2. **(B) Not reached, but the route is redrawn and measured end to end.**
   - The stage-2 starvation of the p = 3 slices was an accounting artifact (§1).
     - For p = 3 the Paley(3) ansatz is vacuous [P]. The repo's stage 1 fixed a bit ("J′ or twisted") that S1 cannot see.
     - With the correct split (stage 1 = the Z₃-coset sums y; stage 2 = Eisenstein phases Z′ = y·ω^{z*}), the ×27 slice
       moves from 2^−39…2^−53 per y to **2^−3…2^−13 per y [E]**.
   - The slice fixes Z₆ = Z₂ × Z₃ pointwise. This gives an **exact four-level decomposition** (§3) [P, M]. On it I
     built a hierarchical exact solver:
     - stage 1 = X₁ = y₀ + y₁, then Xm = y₀ − y₁;
     - stage 2 = S′ = Z′₀ + Z′₁, then D′ = Z′₀ − Z′₁;
     - bordered GS completion and independent verification.
   - **Positive controls [M]:** it constructs verified **H(172)** (v = 43) and **H(316)** (v = 79) in multiplier slices
     fixing Z₆ (§5), and a planted test recovers a known core.
   - **At 571 [M]:**
     - stage 1 works once the exact compression arithmetic of §4 is built in: 1003 exact S1 solutions in 7 s;
     - the remaining cost is stage 2.
   - **Cost at 571 [E]:** screening all 1651 S1 solutions found, the Gaussian gives an S-level join of 2^44.6–2^47.4
     pairs per y and a yield of 2^−1.1 … 2^−8 per y. That is **2^47.6 per core for the best y, 2^50.9 for the median**
     (§6.1).
   - **Calibration at 571 [M]:** the exact S-level join on a real y finds **1, 0 and 0 S-solutions in three runs where
     the Gaussian predicts 2^3.5 per run**, so the Gaussian is optimistic by about 2^5 (Poisson 95 %: at least 2^2.6).
     A planted test at full 571 size recovers 3 of 3 planted solutions, so the join is not the cause (§5.3).
   - **Calibrated cost: about 2^52–2^56 per core**, 7–11 bits above the 2^45 target. This is out of reach here.
   - The reason is specific: 571's stage 2 has no pointwise-fixed factor left to give the cheap class split that made
     573 work (§6.3).

3. **(C) Exact obstructions [P, M].**
   - **θ₅ criterion (new) [P].** If 10 | v − 1 and the multiplier group maps onto (Z/5)*, a b = 1 core needs
     **v = (5s₀ + 1)² + 25(s₁² + s₂² + s₃²)**.
     - It kills, exactly, slices with large Gaussian estimates: the **211 ×27 analogue** (+70 [E]) and 18 small-order
       slices (§4.1).
     - **571 passes it**, with s₀ = 2 forced.
   - **The Z₁₉-cubes slice at v = 115** (a smaller order with the same Z₁₉ structure) is **exactly empty** [M]:
     exhaustive stage 1 (2304 y's), then the full stage 2. Its Gaussian estimate was +22.6.
   - **For 571 ×27 the exact arithmetic thins the slice but does not empty it** [M]:
     - exactly 19,200 compressed (Z₁₀) profiles satisfy the four rational-character equations;
     - 2304 of them survive the exact lag-0 / lag-95 test (§4.3), carrying weighted mass 2^204.3 of 2^237.7;
     - real S1 solutions exist inside them.
   - So the ×27 slice is alive, and its difficulty is stage 2.
   - **571's second-best slice, ⟨61⟩ = Paley(19) on Z₃₀ (Gaussian +39.9), is mostly empty [P, M] (§7).**
     - Exact lemmas reduce it to 70,696 gadget classes, with (K, N₋) = (6, 1) or (7, 6).
     - **All 10,008 classes with (K, N₋) = (6, 1) have no stage-1 solution, even before any S2 condition.**
       - 9,118 of them (91 %) fail a single lag in a solver-free exact DP.
       - All 10,008 are infeasible in two independent CP-SAT formulations.
       - So no core of this slice has six constant cosets.
     - For (7, 6), 49,570 of the 60,688 classes have no stage-1 solution. The other 11,118 are unresolved.
     - Method: CP-SAT, cross-checked by a second formulation, a solver-free lag DP, planted tests, and a positive
       control on the published 419 core.
     - The exact-marginal estimate for the unresolved classes is about 2^−12 stage-1 solutions in total [E]. The same
       estimator gives 2^18.3 on the (non-empty) published 419 class.

---

## 1. The p = 3 re-split (Q1, Q2) [P, M]

**Proposition 1 [P].**
- For p = 3 the squares mod 3 are {1}. So every ±1 triple (x₀, x₁, x₂) on a Z₃-coset is αδ + uχ + wε:
  - x₁ = x₂ gives the constant type, w = x₁;
  - x₁ = −x₂ gives the twisted type, u = x₁.
- Put y = x₀ + x₁ + x₂ ∈ {±1, ±3} and Z = x₀ + ωx₁ + ω²x₂ ∈ Z[ω].
  - |y| = 3 ⟺ Z = 0.
  - If |y| = 1 there is exactly one minority position z*, and Z = −2y·ω^{z*}.
- The Z₃-Fourier transform of the b = 1 condition is exactly:
  - **S1:** Σ_k PAF(y_k) = T_J (4m − 8 at 0, −12 elsewhere);
  - **S2:** Σ_k PAF_C(Z′_k) = v·δ, where Z′ = Z/(−2) = y·ω^{z*} on |y| = 1 and 0 on |y| = 3.
- x ↔ (y, z*) is a bijection.

**Why the old split starved [P].**
- A J′ coset (α, −α, −α) and a twisted coset (−α, u, −u) have the same y. The repo's stage 1 enumerated
  (gadget J/J′, α), and so fixed the bit "z* = 0 or z* ≠ 0" before S1, which cannot see it.
- Stage 1 was therefore inflated by the number of J′-placements, and stage 2 was restricted to the binary choice
  z* ∈ {1, 2}.
- The correct split:
  - stage 1 = distinct y;
  - stage 2 = z* ∈ Z₃ on every |y| = 1 orbit (ternary).

**Exact lattice of stage 2 [P].**
- Z′ ≡ y (mod 1 − ω), and 3 ≡ 0 (mod (1 − ω)²). So Σ_k PAF_C(Z′_k)(h) ≡ Σ_k PAF(y_k)(h) = T_J(h) ≡ 0 (mod 1 − ω).
- Hence for every S1 solution the S2 target lies in the right residue class, and every lag contributes an index-3 cell.
- The naive estimator returned −∞ on random (non-S1) y for exactly this reason.

**Verified [M]** (`theory/p571_ternary.py selftest`):
- S1, S2 and the bijection hold exactly on the b = 1 cores v = 43 (Z₄₂) and v = 85 (Z₈₄);
- the level split holds on random slice sequences (×27, flip slices, v = 211, 391 skeletons).

**The corrected split for 571 ×27 [E]** (`theory/p571_ternary.py 570 103`, `theory/logs/p571_ternary_split_103.jsonl`):

| | repo split (`logs/split_571_p3_mu27_partial_record.log`) | y \| Z split (this report) |
|---|---|---|
| stage 1 | +37 … +49 per (K, N₋): counts (y, J′-placement) | **+44.0**: distinct y |
| stage 2 per y | −39 … −53: binary u | **−2.6 … −13.1**: ternary z* |
| stage-2 lists per block | 2^{61–68}/4 bits | 3^{n_free,k}: 2^14–2^30 |

## 2. All multiplier slices of Z₅₇₀ [M, E]

- (Z/570)* ≅ Z₂ × Z₄ × Z₁₈ has **81 subgroups** (`theory/p571_msweep.py`, `theory/logs/p571_msweep_list.jsonl`).
- The repo sweep saw only "Paley(p) + cyclic ⟨μ⟩".
- Gaussian log₂E for every subgroup with 16–90 orbits is computed by a general estimator (`theory/p571_slicelib.py`).
  It reproduces the repo numbers: 571 ×27 +40.5 (repo +40.6), 573 ×29 +23.8 (+23.2), 419 ×7 +20.6 (+21.1).
- Results: `theory/logs/p571_msweep_est.jsonl`.

| M (generators mod 570) | order | orbits / block | equations | log₂E | residues mod 3 / 5 / 19 | comment |
|---|---|---|---|---|---|---|
| ⟨103⟩ = the ×27 slice | 12 | 66 | 37 | **+40.5** | 1 / all / cubes | best; fixes Z₆ |
| ⟨61⟩ | 9 | 90 | 45 | +39.9 | 1 / 1 / squares | = full Paley(19) on Z₃₀; K = 6 part exactly empty, K = 7 mostly (§7) |
| ⟨11, 49⟩ | 12 | 84 | 47 | +22.5 | ±1 / ±1 / {1,7,11} | = Paley(5) + ×11; contains the pure flip, so Z′ is real and z* ≡ 0 |
| ⟨83⟩ (new: flip on Z₃) | 12 | 70 | 37 | +18.3 | ±1 / all / {1,7,11} | stage 2 per y −39 … −46 [E] |
| ⟨11, 31⟩ | 12 | 80 | 47 | +2.0 | ±1 / 1 / cubes | |
| ⟨107⟩ | 12 | 64 | 43 | −1.6 | ±1 / all / cubes | |
| all others | | | | < −7 | | |

- No new slice beats ×27.
- Flip slices (elements ≡ −1 mod 3) tie z* across orbit pairs, or force z* = 0. They lose stage-2 freedom faster than
  equations (`theory/logs/p571_ternary_split_{83,107}.jsonl`).

## 3. The exact Z₆-level decomposition and the hierarchical solver [P, M]

**Proposition 3 [P].**
- Let M fix Z₆ pointwise (t ≡ 1 mod 6), and let Z_m = Z₂ × Z₃ × Z_q with q = m/6.
- For each character λ = (ε, η) of Z₆ put c_λ(t) = Σ_{a,b} ε^a η^b x(a, b, t). Then x̂(λ, θ) = ĉ_λ(θ), so the core
  condition is equivalent to four level conditions on Z_q:

| level | component | condition | role |
|---|---|---|---|
| L1 | X₁ = y₀ + y₁ | Σ_k PAF(X₁) = 4v·δ − 24 | S1 (even) |
| L2 | Xm = y₀ − y₁ | Σ_k PAF(Xm) = 4v·δ | S1 (odd) |
| L3 | S′ = Z′₀ + Z′₁ (Eisenstein) | Σ_k PAF_C(S′) = v·δ | S2 (even) |
| L4 | D′ = Z′₀ − Z′₁ | Σ_k PAF_C(D′) = v·δ | S2 (odd) |

- The levels are coupled only through the per-(block, orbit) ±1 symbol x(·, ·, O) ∈ {±1}⁶.
- In the slices used here (×27 at 571, ×27 = (2, 6) at 211) there are 11 orbits on Z_q and 7 nonzero lag classes.
- Two exact links between levels:
  - S′ ≡ D′ (mod 2);
  - S′ ≡ X₁ and D′ ≡ Xm (mod 1 − ω).
- Ignoring the first link overcounts stage-2 constraints by about 10 bits (§6.2).
- Verified [M] (`theory/p571_levels.py check`) on the v = 43 core and on random slice sequences.

**Per-orbit entropies (uniform symbol) [P]:** H(X₁) = 2.33, H(Xm | X₁) = 1.29, H(S′ | y) = 2.04, H(D′ | S′, y) = 0.34
bits. Stage 2 is 85% S-level.

**Solver** (`theory/p571_hsolve.py`):
- **Stage 1:**
  1. Exact compressed profiles (§4).
  2. Per block, the y-configurations of the profile; the four fibre lists factor exactly.
  3. L1 join on X₁ classes.
  4. For each L1 solution, the L2 join over all profile-consistent refinements.
- **Stage 2:**
  1. Per block, the S′-classes: the product over orbit pairs of the distinct S′ = Z′₀ + Z′₁ values.
  2. S-level join.
  3. For each S-solution, the D-level check over the few D′ options per class.
  4. Assembly, then an exact core check by FFT.
- **Completion:** `theory/p571_finish.py`. This is bordered GS with circulants, the forced corner and exact H Hᵀ = N I,
  written independently of `theory/paley_b1pipe.py`.
- **Joins:** numba 2+2 meet in the middle on a random linear 64-bit hash. The residue passes are bucketed by a separate
  additive hash modulo a prime, because `key mod M` is not additive under 64-bit wraparound. Every hit is re-checked
  exactly.
- **Tests [M]:**
  - `theory/p571_test_join.py`: planted quadruples are found;
  - `theory/p571_test_stage2.py`: planted S-classes recover the original v = 43 core.

## 4. Exact arithmetic of the rational characters [P, M]

The Gaussian treats every equation as a continuous coordinate. At the **rational** characters of a slice the
equations are sums of squares over small integer data, and that is where it fails.

### 4.1 The θ₅ criterion [P]

**Theorem 4.**
- Let 10 | m = v − 1, and let a b = 1 core on Z_m be invariant under a multiplier group M that maps onto (Z/5)*.
- Then v = (5s₀ + 1)² + 25(s₁² + s₂² + s₃²) for some integers s_i.

Proof:
1. Let W_k(c) = Σ_{i ≡ c (5)} x_k(i) be the Z₅-compression.
2. M-invariance makes W_k constant on c ≠ 0. So for χ of order 5,
   x̂_k(χ) = W_k(0) − W_k(1) = (5W_k(0) − r_k)/4 ∈ Z, where r_k = W_k(0) + 4W_k(1).
3. W_k(1) is a sum of m/5 (even) ±1's, so it is even. Hence W_k(0) ≡ r_k (mod 8).
4. With r = (±2, 0, 0, 0) (b = 1), the core condition at χ, Σ_k x̂_k(χ)² = 4v, becomes
   4[(5s₀ ± 1)² + 25Σ_{k≥1} s_k²] = 4v. ∎

**Consequences [M]** (`theory/p571_controls.py`, `theory/logs/p571/controls_scan.jsonl`):
- **571:** 571 = 11² + 25·18 with 18 = 3² + 3² + 0² = 4² + 1² + 1². s₀ = 2 is forced (the case 5s₀ + 1 = −14 needs 15 as a sum
  of three squares, but 15 ≡ 7 mod 8). So **571 passes** for every slice mapping onto Z₅*. In the ×27 slice:
  - the Z₅-compression of X₁ is pinned to A₀ = 18;
  - {A₁, A₂, A₃} = {±24, ±24, 0} or {±32, ±8, ±8}.
- **211:** the only candidate s₀ = 1 needs 7 to be a sum of three squares, so **every 211 slice mapping onto Z₅* is empty**,
  including the same-skeleton ×27 analogue (Gaussian +70). The general per-prime test (`theta_p_test`) also kills slices
  at v = 151, 295, 331 and 391 (19 slices in all).
- Random lists at 571 essentially never meet the pinned A-values. This is why the first L1 joins found nothing; stage 1
  must enumerate the profiles first (Đoković–Kotsireas compression).

### 4.2 Z₁₀-compression profiles of S1 [P, M]

- Z₁₉₀ = Z₂ × Z₉₅ has four rational characters under ×27: (±1) ⊗ (θ₀, θ₅).
- Per block they see four odd integers:
  - Y(a, 0) = Σ y(a, ·) over the c = 0 fibre;
  - Y(a, *) = Σ y(a, ·) over one c ≠ 0 fibre.
- The conditions are:
  - (i) rows r_k = Y(0,0) + Y(1,0) + 4(Y(0,*) + Y(1,*)) = (2, 0, 0, 0);
  - (ii) Σ_k[(Y(0,0) + Y(1,0)) − (Y(0,*) + Y(1,*))]² = 4v;
  - (iii) Σ_k[(Y(0,0) − Y(1,0)) + 4(Y(0,*) − Y(1,*))]² = 4v;
  - (iv) Σ_k[(Y(0,0) − Y(1,0)) − (Y(0,*) − Y(1,*))]² = 4v.
- The c = 0 fibre meets only the orbits o₀ and o_{C_i}; the c ≠ 0 fibre meets only o₅ and the six mixed orbits. So the
  reachable values factor exactly.
- `theory/p571_z10profiles.py` enumerates **all** solutions:
  - 571 ×27: **19,200 profiles**, weighted mass 2^209.5 (x-multiplicity weights; all row-sum-consistent configurations:
    2^237.7);
  - 211 ×27: **0 profiles**, consistent with Theorem 4.

### 4.3 The lag-0 / lag-q filter [P, M]

- Inside a profile, the y-values on the small orbits are pinned. For example Y(a, 0) ≡ y(a, o₀) + 6 (mod 12) fixes y(a, o₀).
- S1 at lag 0 forces the J-count Σ_k J_k = 189 (Σy² = 760 + 8J).
- S1 at lag q = 95 forces Σ_k ρ_k = −6, where ρ_k = Σ_t y_k(0, t) y_k(1, t).
- Per block and profile, the reachable (J, ρ) set is computed exactly: group lists, pairing of the a = 0 and a = 1
  configurations, OR-convolution (`Lag0` in `theory/p571_z10profiles.py`, `theory/p571_lag0_run.py`).
- **Result [M]:**
  - **2304 of 19,200 profiles survive**, with weighted mass 2^204.3;
  - the 400 heaviest profiles all fail;
  - the survivors are exactly the profiles feasible mod 16 (independent DP over orbits, `theory/p571_empty_check.py 16`,
    same 19,200 profiles).
- The slice is not empty at this level.

### 4.4 Stage 1 at 571 works [M]

`theory/p571_hsolve.py stage1z 570 103 --profiles theory/logs/p571/z10_lag0_all.pkl`, run on the heaviest surviving
profile:
- L1 lists of 4096 X₁-classes per block give **22,550 L1 solutions**;
- L2 gives **1003 distinct y**, every one an exact S1 solution (`y_to_x_check`), **in 6.6 s** (`theory/logs/p571/y571_zl0.npz`).

With a pairing bias (`--lam 1.5/3`), 334 and 314 more y's have more both-free orbit pairs (mean 24.4 and 25.0, against
21.9).

## 5. Controls [M]

### 5.1 Positive controls (constructions)

`theory/p571_run.py`, run end to end (stage 1, stage 2, bordered GS, independent check):

| v | slice (on Z_m, Z₆ fixed) | stage 1 | stage 2 | result |
|---|---|---|---|---|
| 43 | ⟨13⟩ (±1 on Z₇), R = 4 | 432,384 y | 1.5·10⁶ S-solutions per y, 6·10⁵ cores per y | **H(172)**, `data/core_v43_b1_p571method_c43.txt` |
| 79 | ⟨55⟩ (order 3 on Z₁₃), R = 5 | 1,002,816 y | 3.3·10⁵ S-solutions per y, 11,988 cores per y | **H(316)**, `data/core_v79_b1_p571method_c79.txt` |

- Both cores re-verified by `python3 theory/p571_finish.py <file>`: H Hᵀ = N I. Both are multiplier-invariant.
- SHA-256 of the packed matrices:
  - 8446577b…;
  - dd464e6c….
- **Larger p = 3 controls, not completed.** Two R = 9 slices, v = 247 (⟨133⟩ on Z₂₄₆) and v = 211 (⟨109⟩ on Z₂₁₀),
  each gave about 10⁶ exact S1 solutions in 3–8 min (`theory/logs/p571/run_247_sub.log`, `run_211_m4_sub.log`). Stage 2
  (2·10⁴-class sublists) ran about 45 min with no core, and then the runs were lost to a worker restart. They are not
  counted as controls.

### 5.2 Exactly empty slices (negative controls)

| v | slice | Gaussian log₂E | exact result |
|---|---|---|---|
| 211 | ×27 = (2 mod 5, 6 mod 7), same skeleton as 571 | +70.2 | **empty**: Theorem 4 (0 Z₁₀-profiles) [P, M] |
| 115 | cubes mod 19 (R = 4) | +22.6 | **empty**: exhaustive stage 1 (2304 y, complete X₁ lists) and full stage 2, 0 S-solutions [M] |

### 5.3 Calibration of the stage-2 Gaussian [M]

- v = 79, 12 random S1 solutions, full stage-2 enumeration:
  - predicted per-y yield 2^14.2–2^15.3;
  - exact core counts 2^11.7–2^13.9;
  - so the Gaussian **overcounts by 2^1.5–2^3.5**.
- **At 571 the S-level count on a real y is measured directly** (`theory/p571_s3calib.py
  theory/logs/p571/y571_zl0.npz 0 65536`, `theory/logs/p571/s3calib_571_y0.log`):
  - y = 0 of `y571_zl0.npz`, with random sublists of 2^16 S-classes per block (the full lists are 2^18.7–2^24.9), and
    three seeds;
  - the Gaussian (C₃ = 60.5) predicts 2^3.5 = 11.3 S-solutions per run;
  - the exact join finds **1, 0 and 0**.
  - So the S-level Gaussian is optimistic by about 2^5. With one count in three runs the Poisson 95 % bound is
    still at least 2^2.6.
- **The join is not the cause [M]:** `theory/p571_test_plant571.py` plants a random quadruple in 2^14-class sublists of
  the same y and recovers it in 3 of 3 trials (about 37 s per trial).
- This 571 offset is larger than the v = 79 offset (2^1.5–2^3.5), so all 571 yields below are quoted uncalibrated and
  then shifted by about 5 bits.

## 6. 571 ×27: stage 2 and the cost [M, E]

### 6.1 Stage-2 numbers on real S1 solutions

`theory/p571_s2est.py`, `theory/p571_screen.py`, `theory/logs/p571/s2est_571_zl0.log`, `…/screen_571.log`:
- n_free = 62–66 ternary orbits, so H₂ = 98–104 bits;
- the per-block S-class lists are 2^15.5–2^27.1, Σh₃ = 88–91, and the S-join is 2^44.6–2^47.4 pairs per side;
- the stage-2 yield per y (Z₁₉₀ Gaussian, exact index-3 lattice) is −1.1 … −8, with median −5.5;
- S-level: about 2^28–2^31 S-solutions per y by the Gaussian. D-level survival per S-solution: about 2^−33.
- **Screen of all 1651 distinct y** (`theory/p571_screen.py`, `theory/logs/p571/screen_571.log`):
  - Gaussian work per core is 2^47.6 for the best y and 2^50.9 for the median;
  - the best y's have n_free = 66 (H₂ = 104.6), Σh₃ ≈ 92 and yield −1.1 … −1.8;
  - the 20 best are saved in `theory/logs/p571/screen_best.npz`.
- **Calibrated [E]:** the 571 S-level offset (§5.3) lowers every yield by about 5 bits. The best y then costs about
  2^52.6 per core, and the median about 2^56.
- **The method has a floor [E].** In the Gaussian model, with sublists chosen optimally, the work per core is
  max(C₂ − h_D − Σh₃/2, (C₂ − h_D)/2). Here h_D = H₂ − Σh₃ is the D-level freedom left after the S-classes are fixed.
  - Only orbit pairs with both a = 0 and a = 1 free carry h_D: log₂(9/6) bits if y₀ = y₁, log₂(9/7) if y₀ = −y₁.
  - Biasing stage 1 towards such pairs (`--lam 3`) raises h_D from at most 12.5 to at most 14.3 bits
    (`theory/logs/p571/y571_lam3.npz`).
  - With C₂ ≈ 106 the floor is then about (106 − 14.3)/2 ≈ 2^46 uncalibrated.
  - In the positive-yield regime the 5-bit offset counts only half, which gives about 2^48.4.
  - So selecting y cannot reach 2^45.

### 6.2 An estimator pitfall [P]

- Treating L3 and L4 as independent Gaussians gives C₃ + C₄ ≈ 121 against C₂ ≈ 105 for the same y.
- The difference is the exact parity link P_{D′} ≡ P_{S′} (mod 2), an index-2 sublattice in each of about 11
  coordinates.
- With it, C(L4 | L3) ≈ 44 bits and the two estimators agree.

### 6.3 Why 571 has no 573-style shortcut [P/E]

- At 573 the Z₄ class split cut the stage-2 join from 2^34 to 2^23–2^25 pairs per y. With 141 y's to the first core,
  that is about 2^32 in total, roughly 2^9–2^11 below the 2^{C₂/2} meet-in-the-middle bound (C₂ ≈ 82–86).
- The reason is that Z₄ is pointwise fixed by ×29. Its characters give exact integer equations on compressions with
  13-point fibres and binary entries, which are low-entropy classes with strong constraints.
- At 571, stage 2 has already used Z₃ for the y | Z split, and its only remaining pointwise-fixed factor is Z₂.
  - Its compressions (Eisenstein sums over 95-point fibres), the Z₁₀-compressions (19-point fibres, 6–7 values per
    orbit) and S′ itself are all nearly injective on the S-classes.
  - The rational characters (θ₀, θ₅) give sums of four Eisenstein norms, which take only O(v) values. So class pairs are
    almost all completable.
- The S-level is a dense 4-list problem with about 89 bits of class entropy and about 60 bits of constraints. A meet in
  the middle costs about 2^44.5 per y.

### 6.4 Cost to a first core [E]

- Stage 1 is negligible: about 1000 y per 7 s.
- Stage 2 is about 2^{join} × 2^{−yield}:
  - Gaussian: 2^47.6 for the best screened y;
  - with the 571 S-level calibration (§5.3): **about 2^52–2^56 pair operations per core**.
- The current numba kernel does 1–1.5·10⁷ pairs per second per core (`theory/p571_test_plant571.py`: 2·2^28 pairs in
  about 37 s). At 4 cores, 2^52 pairs is about 2^52 / 6·10⁷ s⁻¹ ≈ 2.4 years. Even uncalibrated, 2^47.6 is about 40 days.
  **So (A) was not attempted.**

## 7. The Paley(19) slice ⟨61⟩: K = 6 is empty, K = 7 reduced to a residue (Q3) [P, M, E]

**What the slice is.**
- 61 ≡ 1 (mod 30) and 61 ≡ 4 (mod 19). So ⟨61⟩ fixes Z₃₀ pointwise and acts on Z₁₉ as the squares.
- Its slice is exactly the Paley(19) ansatz on G′ = Z₃₀: on each Z₁₉-coset g of each block, x is constant on QR and
  on QNR.
- The Gaussian gives +39.9, second only to ×27 (§2). The repo's earlier split had put all its cost in a "rigid"
  stage 1 (`report/paley_design_answer.md` §7.4).

**Notation.**
- A coset is either twisted (x = a at 0, ±u on QR/QNR) or a gadget (x = a at 0, w on all of Z₁₉*).
- A_k(g) = a; W_k(g) = w on gadgets and 0 on twisted cosets.
- y_k = A_k + 18W_k is the Z₁₉-compression. A gadget is J if a = w (y = 19w) and J′ if a = −w (y = 17w); write ε = +1
  for J and −1 for J′.
- **S1:** Σ_k PAF(y_k)(h) = −76 for h ≠ 0, Σ_k Σ_g y_k(g)² = 2208, and row sums r = (2, 0, 0, 0) after reordering and
  negating blocks.
- **S2:** with ẑ_k = (A_k − W_k) + u_k√−19 on Z₃₀ (u = 0 on gadgets), Σ_k PAF_C(ẑ_k) = 4v·δ.

**Lemma 7.1 [P].** The number of gadgets K and the number of J′ gadgets N₋ satisfy (K, N₋) ∈ {(6, 1), (7, 6)}.
- Proof: lag 0 gives n_tw + 289N₋ + 361N₊ = 2208 with n_tw + K = 120, so 4N₋ + 5N₊ = 29.

**Lemma 7.2 [P].** t_k = Σ_g W_k(g) ∈ {−1, 0, 1}.
- Proof: Σ_g A_k(g) = r_k − 18t_k and |Σ_g A_k(g)| ≤ 30.
- So a block with an odd number of gadgets has Σ A_k ∈ {−16, ±18, 20}: only 5–7 minority entries out of 30.

**Lemma 7.3 [P] (gadget autocorrelation and X-sums).**
- Expand PAF(y) = PAF(A) + 18X + 324·PAF(W), where X_k(h) = Σ_{g ∈ supp W_k} W_k(g)(A_k(g+h) + A_k(g−h)) is a sum of
  2K_k terms ±1.
- Then for every h ≠ 0, s(h) := Σ_k PAF(W_k)(h) ∈ {0, −1}, and s(15) = 0.
- Also x_h := Σ_k X_k(h) ∈ [8, 2K] where s(h) = −1, and x_h ∈ [−10, 2] where s(h) = 0.
- Proof: in Σ_k PAF(A_k) + 18x + 324s = −76 we have |Σ_k PAF(A_k)| ≤ 120 and |x| ≤ 2K ≤ 14. A pair at distance 15
  counts twice, so s(15) is even.
- Since Σ_{h≠0} s(h) = Σt_k² − K, exactly K − #odd lags have s = −1, in ± pairs, where #odd counts the blocks with
  odd K_k. Call these the special lags.

**Lemma 7.4 [P] (S2 at the two real characters).**
- At the trivial character: Σ_k (r_k − 19t_k)² + 19Σ_k ρ_k² = 4v, with ρ_k = Σ_g u_k(g) ≡ 30 − K_k (mod 2).
  - With |t_k| ≤ 1 this rules out four odd K_k: Σρ² would be 48 or 40, but four odd squares sum to 4 (mod 8).
- At ψ = (−1)^g: Σ_k Ŵ_k(2Â_k + 17Ŵ_k) = Σ_k û_k², with û_k ≡ 30 − K_k (mod 2) and |û_k| ≤ 30 − K_k (hats taken at
  ψ).
- Proof: at a real character the Gauss sum √−19 is imaginary, so S2 reads Σ_k [β̂_k² + 19û_k²] = 4v with β = A − W.
  Subtract S1 at ψ, which is Σ_k (Â_k + 18Ŵ_k)² = 4v.

**Lemma 7.5 [P] (forced terms).**
- If g and g + h are gadgets of the same block, the two X-terms that join them are fixed: ε(g+h)W(g)W(g+h) and
  ε(g)W(g)W(g+h).
- For K = 6 (one J′ gadget), every special lag has a forced −1. So x_d ≤ 10 there and Σ_k PAF(A_k)(d) ≥ 68.

**Computation [M].**
1. **Gadget classes.** `theory/p571_p19gadgets.py` lists every gadget configuration (W, ε) that satisfies Lemmas
   7.1–7.5 at the gadget level.
   - Symmetries used: per-block translation and reflection, negation of blocks 1–3, permutations of blocks 1–3, and
     the global multipliers of Z₃₀.
   - Result: **70,696 classes** (10,008 with K = 6 and 60,688 with K = 7) in 5.4 min
     (`theory/logs/p571/p19/gadgets.log`).
   - **Cross-check [M]** (`theory/p571_p19gadget_check.py`, `theory/logs/p571/p19/gadget_check.log`): an
     independent brute force reproduces the enumerator's class sets exactly on 6 distributions, 3,456 classes in all.
     It enumerates every signed set with a gadget at 0 in each non-empty block, with no excess pruning and no join,
     applies the conditions directly and canonicalises.
     - (0,2,2,2): 296; (1,2,2,1): 224; (2,2,1,1): 216; (0,3,2,1): 456;
     - (1,2,2,2): 800; (2,2,2,1): 1,464.
2. **Stage 1 for every class (exact).** `theory/p571_p19cpsat.py` builds one OR-tools CP-SAT model per class:
   - Booleans for A on the twisted cosets, with A = εW on gadgets;
   - the row sums;
   - the 15 lag equations exactly (one XOR literal per product of two twisted cosets);
   - Lemma 7.4 at ψ as a domain constraint.

   **Result [M]** (`theory/logs/p571/p19/verdicts_pass1.jsonl.gz`, `status_pass1.log`; 1 worker, 0.3–3 s per class,
   about 2 CPU-hours):
   - **(K, N₋) = (6, 1): all 10,008 classes INFEASIBLE.**
   - **(K, N₋) = (7, 6): 49,570 of 60,688 INFEASIBLE; 11,118 UNKNOWN** (`theory/logs/p571/p19/unresolved.txt`).
     - Fully refuted distributions (K₀, K₁, K₂, K₃): (0,5,1,1), (0,5,2,0), (0,6,1,0), (0,7,0,0), (7,0,0,0).
     - The residue sits where three or four blocks carry 1–2 gadgets (most of all in (1,2,2,2), (1,3,2,1), (2,2,2,1)
       and (2,3,1,1)). There the forced X-terms are J′J′ pairs (+1), so the special lags are not extreme and the
       problem is generic.
   - No class was found feasible.
3. **Validation.**
   - Partially planted instances (all but 24 A-values fixed to a random A, targets computed from it) are solved and
     reproduce the planted lag sums (`p571_p19cpsat.py testpartial`, `testy`).
   - Lag-subset models behave consistently. For example, the class in which K₁ = K₂ = K₃ = 2 at distances 3, 5 and 9
     with one J′ is infeasible from lag 5 alone, and this is checked by hand. Lag 5 forces x₅ = 8 and Σ PAF(A)(5) = 104,
     which allows only one minority arc per block along the 6-cycles of the shift. In block 1 the 6-cycle through the
     J′ gadget at 0 is then mixed, with a single minus at 0, and the other 14 minus entries cannot fill whole 6-cycles.
   - **Second formulation [M]** (`p571_p19crosscheck.py y 400 5`, `theory/logs/p571/p19/crosscheck_y.jsonl`). Model y
     is written directly in y = A + 18W: no X/s algebra and no Lemma 7.4, so it is pure S1. On 400 random refuted
     classes it gives **400/400 INFEASIBLE** (29 s in all). So on this sample, S1 alone has no solution.
   - **Solver-free certificates [M]** (`p571_p19crosscheck.py dp 400`, `theory/p571_p19lagdp.py`).
     - For one lag at a time, a DP along the cycles of g ↦ g + h computes exactly the set of σ_k(h) = PAF(A_k)(h) +
       18X_k(h) per block, over all A_k with the right row sum. It agrees with random sampling of A.
     - **130 of 400 random refuted classes (32.5 %) already fail a single lag.** The others need several lags at once.
   - **Positive control [M]** (`theory/p571_p19control419.py`, `theory/logs/p571/p19/control419.log`). Both
     formulations, generalised to Paley(p) on Z_n, are OPTIMAL on the published 419 core (p = 11, n = 38, K = 16,
     N₋ = 11) with its A fixed.
     - From scratch on that class they stay UNKNOWN after 600 s with 4 workers. CP-SAT refutes rigid classes quickly
       but does not find solutions of generic ones.
     - So UNKNOWN must not be read as "probably infeasible", and nothing below relies on it.
     - Simulated annealing on stage 1 (`theory/p571_p19anneal.py control419 300`) does not find a solution of that
       class either (best cost 32 in 300 s). So it is not used as a probe at 571: a failure there would mean nothing.
   - **Smaller p = 19 orders do not calibrate 571 [M].** The Paley(19) ansatz at v = 115, 191, 267 and 305
     (`theory/p571_p19small.py`, `theory/logs/paley_b1pipe_v*_p19.log`) is empty already at the gadget level. For
     n ≤ 16 the bound behind Lemma 7.3 leaves no room for a special lag. The 571 regime, with K ≥ 4 and special lags
     that can be met, starts only at n ≈ 22–28.

**Conclusion (C) [M, E].**
- **[M] No core in the ⟨61⟩ slice has (K, N₋) = (6, 1). Stage 1 alone is already empty for all 10,008 classes.**
  - **9,118** classes have a solver-free single-lag certificate (`theory/logs/p571/p19/lagdp_k6_all.jsonl`). The lag
    that fails is 1 in 2,874 classes, 2 in 2,436 and 3 in 894; the rest are spread over lags 4–14.
  - The other **890** are infeasible both in model A (with Lemma 7.4) and in model y (pure S1). Model y is run on all
    10,008 classes (`theory/logs/p571/p19/modely_k6_all.jsonl`, at most 0.47 s per class).
- **[M] For (K, N₋) = (7, 6), 49,570 of the 60,688 classes are refuted the same way.**
- **11,118 classes are unresolved after pass 1.** A 4-worker, 10 s pass on a random sample of 200 of them refutes
  154 of 200 (77 %; half of those take under 0.4 s, and the other 46 stay UNKNOWN at 10 s; `crosscheck_pass2.jsonl`).
  So a full second pass (about 3 h on 4 cores) would leave roughly 2,600 classes.
- **[E] The residue is very likely empty, but this is not proved.**
  - Estimator: for each unresolved class, N_A·∏_h f_h, where f_h is the exact fraction of row-sum-consistent A
    meeting lag h (`p571_p19lagdp.py estimate`, `theory/logs/p571/p19/estimate_residue.log`). This treats the lags as
    independent but uses exact marginals.
  - Result: on a random sample of 100 unresolved classes, the median is 2^−39.2 and the maximum 2^−20.0 stage-1
    solutions per class, and no class is impossible at a single lag. That is about 3·10⁻⁴ ≈ 2^−11.8 expected stage-1
    solutions over all 11,118. On 100 refuted classes the median is 2^−57.1, and 28 of them are impossible at a
    single lag.
  - The same estimator gives 2^18.3 on the published 419 class, which does contain solutions. In the repo, 419
    gadgets with predicted yield around 2^−27 never gave a stage-1 solution (400 tries,
    `report/paley_design_answer.md` §7.5).
- **So ⟨61⟩, the second-best slice of 571 (+39.9), is mostly empty.** Exactly: every core in it would need seven
  constant cosets in one of the 11,118 listed configurations.

**First computation to run (⟨61⟩ residue).** Re-run `theory/p571_p19cpsat.py --only theory/logs/p571/p19/unresolved.txt
--workers 4 --time 10`, which removes about three quarters of the residue. Then give the survivors to a dedicated exact solver,
for example one per-class join on the Z₂-level split of Z₃₀ (y₊ on Z₁₅, then y₋). Classes that stay UNKNOWN after
120 s with 4 workers, such as the (1,2,2,2) ones, need such a solver.

## 8. Q3 (the 19 factor) and Q4 (b ≥ 3)

**Q3: cubic or sextic twists.**
- In ×27 the Z₁₉-part of the orbits is already the cubic classes C₀, C₁, C₂. On the mixed orbits ×27 adds exactly the
  Jacobi split (c/5)(d/19) (§3, `theory/p571_levels.py`).
- A "cubic Paley" restriction of stage 2 (z* following the cubic class index) is a sub-slice. On a Z₁₉ coset with values
  αδ + uχ₃ + wε, the PSD at the three cubic character classes is |β + uG(χ₃)ω^{−j}|². These are three distinct values:
  one norm equation plus one complex cross equation.
- So the twist removes variables but no equations [P], and log₂E drops.
- The smallest order with the same Z₁₉ structure is v = 115 (Z₁₁₄ = Z₆ × Z₁₉, with M acting as the cubes on Z₁₉ and
  fixing Z₆). That slice is **exactly empty** [M] (§5.2). It is a smaller analogue, not a sub-slice of 571.
- Plain doubly-quadratic twists were already shown empty or negative (`report/paley_design_answer.md` §5).
- The quadratic Paley structure on the 19 factor alone is the ⟨61⟩ slice. It is treated exactly in §7: its K = 6 part
  is empty and its K = 7 part is reduced to a residue.
- Answer: no.

**Q4: b ≥ 3.**
- v + 1 = 572 = 4·11·13 gives c ∈ {4, 22, 26, 44, 52}. c = 4 is dead in the Paley ansatz; no borders are known for
  c ≥ 10.
- Not pursued here. Theorem 4's argument applies verbatim to bordered cores, with targets (T₀ − T_s) in place of 4v.

## 9. What would change the picture

1. **A faster S-level join.** The current kernel does 1–1.5·10⁷ pairs/s/core. A radix-partitioned sort-merge or a GPU
   kernel at 10⁹ pairs/s would put 2^48 at a few days. But the calibrated cost (2^52–2^56) would still take months to
   years. Speed alone does not close the gap.
2. **Selecting y.** The screen picks the best of 1651 y at 2^47.6 [E, uncalibrated]. The method has a floor of
   (C₂ − h_D)/2 (§6.1), so more D-level freedom h_D is the only lever. That would need S1 solutions with about 30 or
   more both-free orbit pairs and y₀ = y₁ on them; the most seen is 29.
3. **A pointwise-fixed factor for stage 2.** The S/D split works because Z₂ is pointwise fixed.
   - The only slices fixing Z₆ with positive log₂E are ×27 and ⟨61⟩ (§2).
   - ⟨61⟩ fixes all of Z₃₀, but its stage 1 is rigid: what is left of it is the residue of §7.
4. **The joint S1+S2 arithmetic.** The rational-character profiles of stage 2 (sums of four Eisenstein norms at θ₀, θ₅
   for L3 and L4) are exact. They are not strong enough for a class split (§6.3), but they are cheap filters on y:
   - y fixes the residues of the stage-2 fibre sums mod (1 − ω);
   - an S1 solution whose compressed stage-2 norm equations have no solution can be dropped before any join.
   - Not implemented.

**First computation to run (×27):**
1. Generate 10⁴–10⁵ y with `stage1z --lam 3` and screen them with `theory/p571_screen.py` (minutes per 10³ y).
2. Run `p571_s3calib.py` on the best-ranked y to measure its S-level offset.
3. Only if the measured cost per core is at or below about 2^47, run `p571_hsolve.stage2` on that y with full S-lists
   (2^45–2^46 pairs, about a week on 4 cores with the current kernel).

## 10. Reproduce

```
python3 theory/p571_slicelib.py 570 217                  # x27 slice estimate (+40.5); also 572 133 157 (573), 418 45 267 (419)
python3 theory/p571_msweep.py 570 --list                  # all 81 subgroups; without --list: Gaussian log2E
python3 theory/p571_ternary.py selftest                   # Proposition 1 on v = 43, 85 cores and random slices
python3 theory/p571_ternary.py 570 103                    # corrected stage split for x27
python3 theory/p571_levels.py check                       # Proposition 3 (Z6 levels), exact
python3 theory/p571_test_join.py; python3 theory/p571_test_stage2.py
python3 theory/p571_controls.py 40 400 --rmax 9           # theta_p tests and small control slices
python3 theory/p571_run.py 78 55 --L 100000 --stop 3      # H(316) control;  42 13 for H(172)
python3 theory/p571_run.py 114 103 --L 1000000 --stop 1000   # v = 115 cubic slice: exactly empty
python3 theory/p571_z10profiles.py 570 103                # 19,200 Z10 profiles (210 97: none)
python3 theory/p571_lag0_run.py 570 103 --top 20000       # 2304 survive the lag-0/lag-95 test
python3 theory/p571_empty_check.py 16                     # independent mod-16 DP: same 2304
python3 theory/p571_hsolve.py stage1z 570 103 --profiles theory/logs/p571/z10_lag0_all.pkl   # 1003 S1 solutions in 7 s
python3 theory/p571_s2est.py theory/logs/p571/y571_zl0.npz 570 103; python3 theory/p571_screen.py theory/logs/p571/y571_*.npz
python3 theory/p571_s3calib.py theory/logs/p571/y571_zl0.npz 0 65536 --seeds 3   # S-level calibration at 571: 1, 0, 0 vs 2^3.5
python3 theory/p571_test_plant571.py                      # planted S-solutions recovered at full 571 size
python3 theory/p571_p19gadgets.py --out theory/logs/p571/p19/gadget_classes.pkl   # 70,696 Paley(19) gadget classes (5 min)
python3 theory/p571_p19gadget_check.py theory/logs/p571/p19/gadget_classes.pkl     # brute-force check of 7 distributions
python3 theory/p571_p19cpsat.py testpartial; python3 theory/p571_p19cpsat.py testy # encoding tests (planted)
python3 theory/p571_p19control419.py                      # both models accept the published 419 core (p = 11, n = 38)
python3 theory/p571_p19cpsat.py theory/logs/p571/p19/gadget_classes.pkl --time 1 --log theory/logs/p571/p19/cpsat_s1_chunk0.jsonl
python3 theory/p571_p19status.py                          # verdicts per gadget distribution; writes unresolved.txt
python3 theory/p571_p19crosscheck.py y 400; python3 theory/p571_p19crosscheck.py dp 400; python3 theory/p571_p19crosscheck.py pass2 200 10
python3 theory/p571_p19lagdp.py estimate theory/logs/p571/p19/gadget_classes.pkl 'theory/logs/p571/p19/cpsat_s1_chunk*.jsonl' 100
python3 theory/p571_finish.py data/core_v79_b1_p571method_c79.txt
```
