# H(2092) (v = 523) in the quartic-Paley slice: exact stage-2 algebra, a calibrated two-level solver, and its cost

Date: 2026-10-10. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_523_design.md`.
Scripts: `theory/q523_*.py`. Logs: `logs/q523/`. Cores: `data/core_v*_b1_q523_*.txt`.

**Tags.**
- **[M]** verified by exact computation here; the script is named.
- **[P]** proved here; the proof is in the text.
- **[E]** estimate. Unless stated otherwise it is the *lattice-exact* Gaussian of `theory/q523_gauss.py`: the
  probability of hitting the target is the Gaussian density of the integer coordinates of the sampled values in a
  basis of the exact lattice they span (no per-coordinate gcd steps). §2.3 shows this lattice is the exact one.
  Calibrations against exact counts are in §3.5.

---

## 0. Bottom line

<!-- BOTTOM-LINE: filled in at the end of the run -->

1. **(A) H(2092): PENDING** (see §6).

2. **(B) A stage-2 method below 2^45 per core, validated.**
   - **Method** (§3): the Z₂ = Z₁₈/Z₉ split of stage 2 into an S level and a D level (both exact quadruple
     conditions over O_K on Z₉), with
     - a **two-level Wagner tree** for the S level: first merge on the Z₃ compression, done at the level of parents
       (row sums), second merge on the order-9 equations; one side is a union of many parents' fibres;
     - a compiled **D-check** for every S-solution: Z₃ class split on D, then the order-18 equations.
   - **Cost at 523 [E, calibrated]**, for the CP-SAT y of `logs/cyclo/stage1_cpsat_523.json`:
     - S level: 2^40.0 pair operations;
     - D-checks: 2^26.5 S-solutions × about 2^17 operations each = 2^43.5;
     - **total about 2^43.6 per core.**
   - **Calibration at full 523 size [M vs E]:**
     - S-step: 8 runs, 477 S-solutions found against 427 predicted; per-run error −1.5 … +0.7 bits;
     - last level (order-18 characters): measured on real S-solutions, the D-success is 2^−26.5 … 2^−26.8 per
       S-solution, against 2^−26.5 from the model.
   - **Controls [M]:**
     - planted 523 instances: the pipeline recovers the planted solution 3/3, at full size;
     - H(316) (v = 79, Z₆ × quartic Z₁₃), verified by two independent finishers;
     - v = 235 (Z₁₈ × quartic Z₁₃, the same shape as 523): PENDING (§4).
   - **The floor of the route** is the D-deficit D_def = C(D|S) − h_D ≈ 26.5 bits (§3.4). It is robust: it does not
     move under any reweighting of the S-lists, because the S-level energy equation pins the D-freedom of
     S-solutions [M]. It depends on y; §5 screens y's.

3. **(C) Exact results [P, M].**
   - **739 ×127 (Z₁₈ × quartic Z₄₁) is exactly empty [P]**: its stage-1 lag-0 equation has no solution (§7).
     The Gaussian had given +57.
   - The stage-2 lattice of the 523 y contains the target (no parity obstruction), and the sampled lattice used by
     every estimate is the exact one [M] (§2.3).
   - No exact obstruction to the 523 quartic slice was found.

4. **New structure [P, M].**
   - Stage 2 is exactly Σ_k Z_k Z_k^† = 4v·δ in O_K[Z₁₈], with Z = Σ_j (c_j − α) η_j in the Gauss-period basis;
     α is a function of y (§1).
   - The coset alphabet splits into a quadratic coordinate e and a quartic coordinate e₁ ∈ 2·Z[i] with
     complementary supports; no stage-2 equation depends on e alone (Q2, §1.3).
   - Per-level constraint budget (§2): the S level costs 83.3 bits, its order-9 part 51.5 bits on 16.4 bits of fibre
     per block; the D level given S costs 63.0 bits on 36.5 bits of freedom.
   - **Energy pinning [M]:** the Z₃ lag-0 equation forces parents of typical fibre size (no "big-fibre" parents), and
     the S-level lag-0 equation forces the D-freedom of S-solutions to its typical value (§3.2, §3.4).

---

## 1. Exact algebra of stage 2 [P, M]

### 1.1 The O_K formulation

- K is the quartic subfield of Q(ζ₂₉). Its Gauss periods η_j = Σ_{z ∈ C_j} ζ^z (j = 0..3) form a Z-basis of O_K,
  with 1 = −Σ_j η_j.
- On a Z₂₉-coset g the slice has x(g, 0) = α and x(g, z) = c_j on C_j. Its Z₂₉-transform is

      Z(g) = α + Σ_j c_j η_j = Σ_j (c_j − α) η_j,     η-coordinates z_j = c_j − α ∈ {0, ±2}.

- Summing the 2-D PAF against ζ^z gives P(h) = Σ_g conj(Z(g)) Z(g + h), with η-coordinates
  R(h, l) − PAF_x(h, 0), where R(h, l) = PAF_x(h, z ∈ C_l).
- **Proposition 1.** Given stage 1, stage 2 is exactly

      Σ_k P_k(h) = 4v·δ_{h,0}   in O_K, for every h ∈ Z₁₈,

  i.e. η-coordinates (−4v, −4v, −4v, −4v) at h = 0 and 0 elsewhere.
- Complex conjugation is η_j ↦ η_{j+2}, because −1 = 2¹⁴ lies in C₂. So K is a CM field with K⁺ = Q(√29).
- **Verified [M]** (`theory/q523_algebra.py`):
  - the η-basis multiplication table against complex numerics, for q = 5, 13, 17, 29 and 41;
  - P(h) against the 2-D FFT PAF on 20 random slice configurations with the 523 y;
  - on the verified H(164) core (v = 41), Σ_k P_k = 4v·δ exactly and stage 1 holds.
- **Counting.** The equations are the orbits of (λ, t) ∈ Ẑ₁₈ × Z₄ under (λ, t) ~ (λ̄, t + 2): **36**, of which one
  linear combination is implied by stage 1 (Parseval on the cosets). So there are 35 independent stage-2 equations.
  - For q ≡ 1 (mod 8), e.g. 739's q = 41, −1 is a fourth power, K is totally real, and the count is **40** (§7).

### 1.2 The coset alphabet [P, M]

- **α is a function of y [P].** y = α + f·Σ_j c_j with Σ_j c_j even, so y ≡ α (mod 2f). Hence stage 2 lives entirely
  in the c-variables, with Σ_j c_j = (y − α)/f fixed per coset.
- **Alphabet [M]** (`q523_algebra.py`), with e = c₀ − c₁ + c₂ − c₃ (quadratic) and e₁ = (c₀ − c₂) + i(c₁ − c₃)
  (quartic):

| y | options | (e, e₁) |
|---|---|---|
| ±1 | 6 | (±4, 0) or (0, 2(±1 ± i)) |
| 13, 15 | 4 | (−2, ±2) or (+2, ±2i) |
| −13, −15 | 4 | (−2, ±2i) or (+2, ±2) |
| ±27, ±29 | 1 | (0, 0): Z = ∓2 or Z = 0 (the coset is dead for stage 2) |

  - e and e₁ have complementary supports.
  - e₁/2 is a unit of Z[i], or (1 + i) times a unit.
  - The Galois generator σ acts by c_j ↦ c_{j−1}, i.e. e₁ ↦ i·e₁. It maps solutions to solutions, so every core
    comes with 4 Galois conjugates.
- **Per-coset Parseval [M]:** y² + f·Σ_t |Z^{(t)}|² = q².

### 1.3 The quadratic subfield (Q2) [P, E]

- Write Z = R + I with R ∈ K⁺ and I ∈ K⁻. Then R = (α − s/4) + (√29/4)·e and I = ¼(G ē₁ − Ḡ e₁), G a quartic
  Gauss sum. The equation splits into
  - Σ_k (R_k R_k* − I_k I_k*) = 4v·δ in K⁺ (rational part and √29 part),
  - Σ_k (I_k R_k* − R_k I_k*) = 0 in K⁻.
- In lag form the rational part is Σ_k [16·PAF(a) + 29·PAF(e) + 58·(PAF(u) + PAF(w))] = 16·4v·δ, with
  a = α − s/4, u = Re e₁, w = Im e₁. The √29 part contains Re(J̄·Σ_g e₁(g) e₁(g+h)) (J a Jacobi sum). The K⁻ part is
  bilinear in (a, e) × e₁.
- **No equation depends on e alone [P].** The u, w terms of the rational part (PAF(u) + PAF(w)) and of the √29 part
  (ReJ·(PAF(w) − PAF(u)) − ImJ·CC_sym(u, w)) are linearly independent quadratic forms. So no linear combination of
  K⁺ equations is free of e₁, and the K⁻ equations are bilinear in e₁.
- **What the quadratic projection constrains:** congruences only. For example (rational part mod 232):
  Σ_k PAF(e_k)/4 + 2·#{both-nonzero e₁ pairs} ≡ (stage-1 data) mod 4.
  - These are exactly the lattice conditions of §2.3, and they hold for the 523 y [M].
  - So a "stage 2a on Q(√29)" has no filtering power beyond parity.
  - The Paley(29) w-screen does not apply: the K⁺ equation carries the −I I* term.
- **Budget split [E]** (`theory/q523_levels.py`): rational + quadratic parts 77.6 bits, quartic (K⁻) part 84.6 bits,
  all 146.2 bits.
  - The field split is not a data split: every norm mixes R and I, so it gives no class keys.

## 2. Constraint budget by level [E] and the exact lattice [M]

### 2.1 Levels

- The 36 coordinates are indexed by the order d | 18 of the Z₁₈-character and the Galois character r of Gal(K/Q).
- Integer projections are Ramanujan sums (`theory/q523_levels.py`):
  - the S level of the Z₂ split is d ∈ {1, 3, 9} (Z₁₈ → Z₉);
  - the D level is d ∈ {2, 6, 18};
  - the Z₃ level is d ∈ {1, 3}.

### 2.2 Budget for the 523 y (`logs/cyclo/stage1_cpsat_523.json`, 62×|1|, 6×|13|, 4×|15|)

| level | equations | bits | per-block data (exact) |
|---|---|---|---|
| all stage 2 | 35 | **146.2** | 45.4 (H₂ = 180.3 in total) |
| Z₁ | 2 | 11.4 | 9.0 |
| Z₃ (d ∈ {1, 3}) | 6 | 31.8 | 19.7 |
| Z₆ | 12 | 57.5 | 29.8 |
| S = Z₉ (d ∈ {1, 3, 9}) | 18 | 83.3 | 36.1 |
| D (d ∈ {2, 6, 18}) | 18 | 83.2 | 36.1 |
| C(D \| S) | | **63.0** | h_D = H(D \| S) = **36.5** in total |

- **Stage-2 yield for this y: +34.1** (180.3 − 146.2). The prompt's +47.3 was an average over y with per-coordinate
  steps.
- The individual d-levels sum to 192 bits against 146 jointly. They are strongly dependent, so conditional bits must
  be used (`theory/q523_chain.py`).

### 2.3 The exact lattice of stage-2 values [P, M]

- **Proposition 2.** For a base configuration c₀, the lattice Λ spanned by all V(c) − V(c₀) (V = Σ_k P_k) is spanned
  by the single-coset first differences and the within-block cross terms.
  - Proof: V is quadratic in the per-coset choices and blocks do not interact. V(c) − V(c₀) = Σ_i Δ_i + Σ_{i<j} X_ij
    over the changed cosets, and every Δ_i and X_ij is itself a difference of values. ∎
- **Exact obstruction.** If T − V(c₀) ∉ Λ, the y has no stage-2 completion.
- **523 y [M]** (`theory/q523_lattice.py`, 6 s):
  - rank 35, target in Λ: no parity obstruction;
  - the lattice spanned by 400 sampled differences has index 1 in Λ, so the Gaussian estimates use the exact
    lattice.

## 3. The method (Q1) [P, M, E]

### 3.1 S/D split, per-parent S-step, and why it costs 2^45

- The pairs (g, g + 9) of Z₁₈ give S(g) = Z(g) + Z(g + 9) and D(g) = (−1)^g (Z(g) − Z(g + 9)) on Z₉.
- Stage 2 ⟺ Σ_k P(S_k) = T_S and Σ_k P(D_k) = T_D on Z₉ [P].
- Per pair the (c, c′) choices give at most 19 distinct S values. The number of D-options given S is
  m(S) ∈ {1, 2, 6} for |y| = 1 pairs.
- First design (`theory/q523_sdsolve.py`): a Z₃ "parent" (row sums of S per block), then a 2+2 join of its fibres on
  the order-9 equations, then a D-check per S-solution.
  - Each parent gives 2^{4f − C_mid} S-solutions at cost 2^{2f}, with f = 16.3 bits per block (typical fibre) and
    C_mid = C_S − C_Z₃ = 51.5.
  - So each S-solution costs 2^18.8. With 2^26.5 S-solutions needed (§3.4) the total is about 2^45.3.
  - This beats the 571 floor formula (C₂ − h_D)/2 = 54.9 by the Z₃ class split, but only to the edge of the target.

### 3.2 No big-fibre parents: energy pinning [M]

- The per-S-solution cost 2^{C_mid − 2f} would fall with larger fibres. Single rows have fibres up to 2^7.45
  (typical 2^6.06), so up to 2^22 per block.
- But the Z₃ lag-0 equation fixes the total row-sum energy Σ_{k,r} |Φ_k(r)|² = 4v, and large fibres sit at small
  row sums.
- An exact DP over the 12 rows (maximal summed log fibre size at the target energy) gives **66.6 bits against 65.35
  typical**, i.e. +0.3 bits per block. So typical fibres are all one gets.
- The first implementation, which chose top-fibre row sums, found no Z₃ parent at all [M].

### 3.3 The two-level Wagner tree (the fix)

- The Z₃ equations are additive over blocks. For a middle value c, every left parent pair (φ₀, φ₁) with
  Z3sig(φ₀) + Z3sig(φ₁) = c combines with every right parent pair with Z3sig(φ₂) + Z3sig(φ₃) = T₃ − c.
- So a side of the second join is the union of the fibre products of all its parent pairs. It is no longer limited
  to one parent's 2^{2f}.
- **Size [E]:** pair Z₃ partial sums carry about 28 bits. So one middle value admits about 2^{2·36.1 − 28} = 2^44 left
  configuration pairs.
- With sides of 2^n pairs a middle value gives 2^{2n − C_mid} S-solutions. 2^{D_def} of them cost 2^{(C_mid + D_def)/2 + 1}.
- Implementation (`theory/q523_wagner.py`):
  - 64 row sums per row (2^18 Z₃ candidates per block); middle values chosen among realised left sums;
  - parent pairs by linear hashing;
  - fibres of up to 24 configurations per row (Gumbel top-K by multiplicity weight);
  - residue-pass join with a cache-sized table (numba), every hit re-checked exactly.

### 3.4 The D-check and the D-deficit

- Given an S-solution, D runs over the product of the D-options (2^9.3 per block). `theory/q523_dcheck_nb.py`
  (compiled, about 6 ms per S-solution at 523):
  1. enumerate the options;
  2. take the Z₃ class split on D (row sums D₃, the d ∈ {2, 6} equations);
  3. run a radix-partitioned 2+2 join of the D₃ values;
  4. for each D₃-solution, join its group members on the full D level.
- Validated against an independent implementation on 24 planted cases at v = 79 and 523 (identical solution sets).
- **Measured on real 523 S-solutions** (`theory/q523_dfine_calib.py`, 34–54 S-solutions per run):
  - 2^7.5–7.7 D₃-solutions per S-solution and 2^10.4–10.8 D-candidates;
  - order-18 equations 37.2–37.3 bits on the real residuals (Mahalanobis 0.2–0.5);
  - **D-success 2^−26.8, 2^−26.7, 2^−26.5 per S-solution** (model: 2^{h_D − C(D|S)} = 2^−26.5).
- **The deficit is pinned [M].** Reweighting the S-lists by m(S)^β, β = 1 … 4, leaves the D-freedom of the
  S-solutions found at 36.0–36.3 bits. The S-level lag-0 equation forces the S-energy to half the total, which
  balances S = 0 pairs (6 D-options) against S-full pairs (1 option).
  - The unconditional model's lower deficits at β > 1 were averages over all configurations, not over S-solutions.
  - A symmetric type tilt (more S = 0 and S-full pairs together) gives D_def = 26.3–26.8 for γ = 1 … 8 [E].
  - The "S-determined D-lag-0" filter rejects nothing in practice (34/34 passed).

### 3.5 Calibration at full 523 size [M vs E]

| run | lists | predicted S-solutions | found |
|---|---|---|---|
| per-parent β = 1 | 2^13.1–13.75 | 2^1.65 | 5 (2^2.32) |
| per-parent β = 1.5 (2 seeds) | 2^13.75 | 2^4.37, 2^3.66 | 12, 10 |
| per-parent β = 2 (2 seeds) | 2^13.75 | 2^4.49, 2^4.29 | 8, 8 |
| Wagner, 2 + 2 parents, logN 28 (2 runs) | 2^28.5 per side | 2^4.72, 2^3.16 | 34, 14 |
| Wagner, 6 + 6 parents, logN 30 | 2^30.1 per side | 2^8.29 | 386 (2^8.59) |

- Total: 477 found against 427 predicted. Per-run log₂ error −1.5 … +0.7, mean −0.2.
- The D-level is calibrated by the measured D-success above (§3.4).

### 3.6 Cost per core (Q1 answer)

- For a y, the route costs, in pair operations,

      2^{(C_mid + D_def)/2 + 1}  (S level, Wagner)  +  2^{D_def} × t_D  (D-checks),   t_D ≈ 2^17.

- **523 (CP-SAT y): S level 2^40.0, D-checks 2^43.5, total 2^43.6 [E, calibrated]**, below 2^45.
  - Wall clock at the measured speeds: about 2·10^7 S-join probes/s and 6 ms per D-check per core. That is about
    200 core-hours, dominated by the 2^26.5 D-checks.
- **Floor of the route:** 2^{D_def}·t_D. D_def does not respond to list reweighting (§3.4), so the levers are y
  (§5) and the constant t_D.

## 4. Controls [M]

- **Planted, full 523 size** (`theory/q523_planted.py`):
  - a random configuration is planted, the target becomes its P, and the planted rows are forced into the lists;
  - the middle value is its left Z₃ sum, joined with other parent pairs of that middle value;
  - **3/3**: the join returns the planted S-solution, the D-check returns the planted D, and the assembly reproduces
    the planted target exactly.
- **H(316), v = 79** (Z₆ × quartic Z₁₃, slice ×55 on Z₇₈ ... F = 6):
  - stage 1 by CP-SAT (0.5 s), stage 2 by the S/D solver (`theory/q523_sdsolve.py`);
  - core `data/core_v79_b1_q523_sdsolve.txt`;
  - bordered GS verified by `theory/q523_finish.py` (independent finisher) and by `theory/paley_b1pipe.py
    verifycyclic`: H Hᵀ = 316·I. SHA-256 of the packed rows: ccef126a….
- **v = 235** (Z₁₈ × quartic Z₁₃: the 523 shape with q = 13): PENDING.

## 5. Choosing y (Q3) [E, M]

PENDING (y-screen `theory/q523_yscreen.py`).

## 6. The 523 run

PENDING.

## 7. 739 (Q5): the quartic slice is exactly empty [P, M]

- At v = 739 (F = 18, q = 41, f = 10) the letters are y ∈ {±1, ±19, ±21, ±39, ±41}, so y² − 1 ∈ 40·{0, 9, 11, 38, 42}.
- Stage-1 lag 0 requires Σ_k Σ_g y_k(g)² = 4m − 4(q − 1) = 2792 over 72 entries, i.e. 9a + 11b + 38c + 42d = 68.
- **There is no solution [P]:** 68, 68 − 38 = 30 and 68 − 42 = 26 are not in the numerical semigroup ⟨9, 11⟩, and
  38 + 42 > 68.
- Confirmed by `theory/q523_stage1.py --list-comps`, `scripts/cyclo_stage1_exact.py 18 41 4` and brute force
  (`logs/q523/v739_quartic_lag0_empty.log`).
- So the ×127 quartic slice of 739 has no stage-1 solution and is empty. The Gaussian had given +57 (stage 1 +29).
- What would have transferred: K would be totally real there (41 ≡ 1 mod 8), with 40 stage-2 equations instead of 36;
  the S/D split, the Z₃ middle value and the D-check are written for any q (`cneg` handles the real case).

## 8. Alternative splits (Q4)

- "y together with α" is the same as y, since α is a function of y (§1.2).
- "y together with the Z₂ compression of c" is exactly the S level of the route.
- Moving D information earlier makes things worse [E] (`theory/q523_chain.py`):
  - fixing D₃ at the parent level (Z₆ parents) raises the cost to 2^42.2 in the per-parent model;
  - adding D₃ to the S-lists to 2^55;
  - each bit of D-freedom fixed early costs two bits in a join and buys at most one bit of deficit.
- Splitting the S level further: Z₉ has no intermediate subgroup. The noncyclic Z₂ × Z₃ × Z₃ gives four line
  compressions, but later lines keep only 1–6 bits of fibre against about 17 bits of equations each, so the
  per-S-solution cost does not drop below the cyclic 2^18.8 in the per-parent model [E].

## 9. What would change the picture

PENDING.

## 10. Reproduce

```
python3 theory/q523_algebra.py                                   # section 1 checks (eta basis, P vs FFT, v = 41 core, alphabets)
python3 theory/q523_levels.py 18 29 4 logs/cyclo/stage1_cpsat_523.json        # level budget (section 2.2)
python3 theory/q523_lattice.py 18 29 4 logs/cyclo/stage1_cpsat_523.json       # exact lattice (section 2.3)
python3 theory/q523_chain.py 18 29 4 logs/cyclo/stage1_cpsat_523.json         # per-parent chains (section 3.1, 8)
python3 theory/q523_calib.py 18 29 4 logs/cyclo/stage1_cpsat_523.json --K 24 --beta 1.5 --seed 1   # S-step calibration
python3 theory/q523_wagner.py 18 29 4 logs/cyclo/stage1_cpsat_523.json --beta 1 --rows 64 --K 24 --logN 30 --ncent 1 --seed 3
python3 theory/q523_dfine_calib.py 18 29 4 logs/cyclo/stage1_cpsat_523.json --beta 1.5 --logN 28   # D-success per S-solution
python3 theory/q523_dcheck_nb.py test; python3 theory/q523_planted.py 18 29 4 logs/cyclo/stage1_cpsat_523.json --trials 3 --logN 29
python3 theory/q523_yscreen.py 18 29 4 logs/q523/y523/*.json                   # y-screen (section 5)
python3 theory/q523_stage1.py 18 41 4 --list-comps                              # 739: no lag-0 composition (section 7)
python3 theory/q523_prod.py 18 13 4 logs/q523/y_v235_a.json --tag v235a --workers 3 --logN 33   # v = 235 control
python3 theory/q523_finish.py data/core_v79_b1_q523_sdsolve.txt                # completion + exact H H^T = N I
```
