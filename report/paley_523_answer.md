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

<!-- status line updated during the 523 run -->

1. **(A) H(2092): not found yet. The 523 production run is in progress (§6).**
   - It runs on 4 workers with the pipeline that produced H(940) at v = 235.
   - Realised cost: about 75 core-hours per expected core (about 3.6 ms per S-solution since the kernel upgrade of
     §3.3, 2^26.1 S-solutions per core), i.e. about 19 hours on this machine. Cores arrive as a Poisson process.
     Status: §6.

2. **(B) A stage-2 method at about 2^44 per core, validated end to end at the 523 shape.**
   - **Method** (§3). The Z₂ = Z₁₈/Z₉ split of stage 2 gives an S level and a D level, both exact quadruple
     conditions over O_K on Z₉. Then:
     - **S level:** a two-level Wagner tree. The first merge is on the Z₃ compression, at the level of parents (row
       sums). The second merge is on the order-9 equations, over a union of parent-pair fibre products, with a
       radix-partitioned, cache-resident join.
     - **D level:** a compiled, row-factored D-check for every S-solution: the Z₃ class split on D, a bitset 2+2
       join, then the order-18 equations with an exact prefilter on lags 0 and 3.
   - **Cost at 523 for the y of §5** [M speeds × E count]:
     - S-solutions needed: 2^{D_def} = 2^26.1;
     - each costs 2^18.1 pair visits (2^{F + 1 − n}, K = 32, n ≈ 37.5) plus one 1.5–1.8 ms D-check;
     - total: 2^44.2 pair visits and 2^26.1 D-checks. Measured: **about 75 core-hours** with the current join
       kernel, 90 with the previous one (the joins find 1.6 bits more S-solutions than predicted; §3.6).
     - The ideal Wagner cost 2^40.7 is out of reach because a Z₃ middle value has only a few hundred parent pairs
       (§3.3).
   - **H(940), v = 235** (Z₁₈ × quartic Z₁₃, the 523 shape; §4) [M]:
     - found by the production pipeline;
     - **the first core came after 1.07·10^7 S-solutions against 1.04·10^7 predicted**, and the whole round gave
       1 core in 15.6M against 1.5 expected;
     - verified by two independent finishers (H Hᵀ = 940·I).
   - **Calibration at 523 [M vs E]:**
     - D-success per S-solution: 2^−26.5 … −26.8 measured, 2^−26.5 from the model (§3.4); for the production y,
       2^−26.02 measured on 1000 production S-solutions against 2^−26.14 from the model;
     - S-solutions per middle value: 1.6–2 bits above the Gaussian prediction at production size, so the prediction
       is conservative;
     - planted full-size 523 instances: 5/5 recovered;
     - H(316) at v = 79.
   - **Floor of the route.** The D-deficit D_def = C(D|S) − h_D ≈ 26 bits is pinned by the S-level energy equation
     (§3.4) and moves by at most about 0.5 bit with y (§5). The S-step per S-solution is set by how many parent
     pairs share a Z₃ middle value (§3.3).

3. **(C) Exact results [P, M].**
   - **739 ×127 (Z₁₈ × quartic Z₄₁) is exactly empty [P]**: its stage-1 lag-0 equation has no solution (§7).
     The Gaussian had given +57.
   - The 523 stage-2 lattice contains the target (no parity obstruction), and the sampled lattice used by every
     estimate is the exact one [M] (§2.3).
   - **Stage-1 structure [M]:**
     - every block PAF value is 2 mod 4;
     - the mod-8 class vectors of the four blocks must XOR to all-ones;
     - the PSDs at j = 9, 6, 3 are integers summing exactly to 4v.
     These explain why naive joins and annealing fail. A residue-aware 4-list join finds many stage-1 solutions (§5).
   - No exact obstruction to the 523 quartic slice was found.

4. **New structure and answers to Q2–Q4 [P, M, E].**
   - Stage 2 is exactly Σ_k Z_k Z_k^† = 4v·δ in O_K[Z₁₈], with Z = Σ_j (c_j − α) η_j in the Gauss-period basis;
     α is a function of y (§1).
   - **Q2.** The quadratic subfield gives parity congruences only; no stage-2 equation depends on the quadratic
     coordinate alone (§1.3).
   - **Q3.** Within a composition, the cost of a y depends only on the pairing of its letters across (g, g + 9), by
     an exact table of H(D|S). Across compositions D_def stays within 25.7–26.5 (§5).
   - **Q4.** Moving D information earlier costs more than it saves (§8).
   - **Energy pinning [M].** The Z₃ lag-0 equation forces typical parents, and the S-level lag-0 equation forces
     the D-freedom of S-solutions to its typical value (§3.2, §3.4).

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
- With sides of 2^n configuration pairs, a middle value gives 2^{2n − F} S-solutions at a cost of 2^{n+1} pair
  visits. Here F is the second-level constraint of the actual lists: C_mid plus whatever the list choice loses.
- **Ideal cost [E]:** one middle value with n = (F + D_def)/2 gives the 2^{D_def} S-solutions of a core for
  2^{(F + D_def)/2 + 1} pair visits.
- **Measured limits [M]** (`theory/q523_wagner.py`, y of §5):

| lists (row sums per row, fibre cap K per row) | parent pairs per middle value | side 2^n | F (fine bits) | F + 1 − n |
|---|---|---|---|---|
| 64, K = 24 | 90 used at n = 34 (about 256 available) | ≤ 2^35.5 | 53.3–53.5 | 19.3 at n = 35 |
| 128, K = 24 | 321–461 (all available) | 2^35.6–35.9 | 53.7–53.8 | 18.9 |
| 256, K = 24 (lazy parents, 10^7 per block) | 182–426 | 2^35.1–35.8 | 53.7 | 18.9 |
| 64, K = 32 (all pairs) | 164–202 | 2^37.2–37.5 | 54.4–54.5 | **18.1** |
| 64, K = 40 (all pairs) | 121–163 | 2^38.5–38.9 | 55.6–55.7 | 18.0 |
| 64, K = 48 | 2 used at n = 34.5 | 2^34.5 | 57.1–57.3 | – |
| full row fibres (K = 1000) | 1 | 2^35.8 | 60.8–61.3 | – |

  - Pair Z₃ partial sums carry 33–38 bits, not the 28 estimated earlier, so a middle value has only a few hundred
    parent pairs.
  - More row sums do not help (`theory/q523_wagner2.py`). Lazy parents with 256 row sums per row (10^7 parents per
    block) still give at most about 430 pairs, because the extra parents have atypical Z₃ signatures.
  - Larger per-row fibres (K) cost F: the top-K choice by multiplicity weight concentrates the S-signatures near the
    target. Each doubling of K adds about 6 bits of side for 3–4 bits of F.
  - The best balance is K ≈ 32–40 with all available parent pairs: n ≈ 37.5–38.9. The ideal n ≈ 40 is out of reach
    with these parents.
  - **Realised S-step per S-solution: 2^{F + 1 − n} ≈ 2^18.1 pair visits (K = 32).**
- **Kernels** (`theory/q523_wagner.py`, all four return identical hits on every test [M]):
  1. `mjoin`: residue passes into a 2^20-entry hash table (about 28 MB, so mostly DRAM latency);
  2. `mjoin_bs`: an L2 bitset prefilter per pass. It is 2.1–2.3× faster, but every pass visits every left fibre
     element, which costs about 2.4·10^6 visits per S-solution at 523 whatever n is;
  3. `mjoin_rad`: one pass per residue class of about 2^23 pairs. The pairs are scattered by their top key bits
     into 16 partitions, and each partition is joined in cache with two bitsets. About 26 ns per pair visit at 523;
  4. `mjoin_rad2` (production from 2026-10-11 23:24 UTC): `mjoin_rad` with two residue classes per pass and 32
     partitions, so each fibre element is visited once per two classes. It also has no modulo in the inner loop,
     puts the smaller side of each partition in the bitset, and clears the bitset by a sweep.
     - On a full-size 523 centre (sides 2^37.6 and 2^37.0, 64 residue classes, three samples) it returns the same
       hit set as `mjoin_rad` (9,142–9,372 hits) and is 1.46–1.61× faster (`theory/q523_joinbench.py`).
     - Whole runs give the same S-solutions as before: v = 235 (13,323, two centres) and 523 (9,438, K = 32,
       logN 33).
     - Profile of `mjoin_rad` per residue class: pair emission 59% (half of it the per-element residue arithmetic
       and range lookups), partition joins 41%.
- **Other details:**
  - 64 row sums per row (2^18 Z₃ candidates per block); middle values chosen among realised left sums;
  - parent pairs by linear hashing;
  - fibres of up to K configurations per row, K = 32 in production (Gumbel top-K by multiplicity weight);
  - every hit re-checked exactly;
  - each S-join's S-solutions and its D-check progress are checkpointed (`theory/q523_prod.py --resume`).

### 3.4 The D-check and the D-deficit

- Given an S-solution, D runs over the product of the D-options (2^9.3 per block). The D-check finds all D that
  satisfy the order-18 equations:
  1. take the Z₃ class split on D: row sums D₃, the d ∈ {2, 6} equations;
  2. join the D₃ values 2+2;
  3. for each D₃-solution, check its group members on the full D level.
- **Row factorisation (`theory/q523_dcheck2.py`).** The D-options of a block are a product over its 9 pairs, and the
  three Z₃ rows {0,3,6}, {1,4,7}, {2,5,8} split this product.
  - So the D₃ groups are exactly the triples of distinct per-row sums (at most 216 row-configurations per row), and
    no configuration of the block has to be enumerated.
  - The Z₃-level signature of a triple is a sum of per-row self terms (lag 0) and per-row-pair cross terms (lag 1).
  - On Z₉, lags 0 and 3 stay inside a row. So they are sums of per-row-configuration terms, cached per S-solution.
    Every member combination is first checked exactly on lags 0 and 3; only survivors get the full signature.
  - The 2+2 join of the distinct Z₃-level signature keys uses a cache-resident bitset prefilter
    (`bitset_join4` in `theory/q523_dcheck_nb.py`).
- **Validation [M]:** identical output to the first compiled version (`theory/q523_dcheck_nb.py`) and to the Python
  reference (`theory/q523_dcheck.py`):
  - 24 planted cases at v = 235 and 523;
  - 122 real 523 S-solutions (almost all with no D-solution), same D₃-solution counts and solution sets.
- **Speed at 523 [M]** (122 real S-solutions, one core): **1.5 ms per S-solution** (was 4.5–7 ms).
  - Split: blocks 0.21 ms, D₃ join about 1.0 ms, expansion 0.2 ms, full level about 0.1 ms.
  - Under the production load (4 workers, 3000 real S-solutions): 1.79 ms, of which blocks 0.27 ms, D₃ join 1.21 ms
    and expansion plus full level 0.36 ms. The D₃ join enumerates about 2^16.9 key pairs per side (378–589 distinct
    keys per block), which is the size of the 2+2 join itself, so little more is to be had there.
- **Bug found and fixed in this session [M].** The first compiled version returned its whole zero-filled solution
  buffer from two early exits (no D₃ match, oversized block). Every S-solution without a D₃ match then reported 64
  all-zero "solutions", and assembly stopped with an alphabet assertion.
  - Both earlier v = 235 control runs died at their first empty D-check.
  - The calibrations below are unaffected: they used the Python reference or planted cases, which always have a
    D₃ match.
  - The test now includes unplanted S-solutions.
- **Measured on real 523 S-solutions** (`theory/q523_dfine_calib.py`, 34–54 S-solutions per run):
  - 2^7.5–7.7 D₃-solutions per S-solution and 2^10.4–10.8 D-candidates;
  - order-18 equations 37.2–37.3 bits on the real residuals (Mahalanobis 0.2–0.5);
  - **D-success 2^−26.8, 2^−26.7, 2^−26.5 per S-solution** (model: 2^{h_D − C(D|S)} = 2^−26.5).
- **Measured on production S-solutions** (y of §5, K = 32 lists; `theory/q523_dsuccess_real.py` on 1000 S-solutions
  drawn from one production join, saved in `logs/q523/ssol_sample_p523f_w2_c7.npy`):
  - 2^7.84 D₃-solutions and 2^11.21 D-candidates per S-solution;
  - order-18 equations 37.23 bits on the real residuals (Mahalanobis 0.05);
  - **D-success 2^−26.02 per S-solution** against the model's 2^−26.14 (a 500-solution subsample gives 2^−26.06).
    So the production y needs about 6.8·10^7 S-solutions per core.
- **The deficit is pinned [M].**
  - Reweighting the S-lists by m(S)^β, β = 1 … 4, leaves the D-freedom of the S-solutions found at 36.0–36.3 bits.
    The S-level lag-0 equation forces the S-energy to half the total, which balances S = 0 pairs (6 D-options)
    against S-full pairs (1 option).
  - The unconditional model's lower deficits at β > 1 were averages over all configurations, not over S-solutions.
  - A symmetric type tilt (more S = 0 and S-full pairs together) gives D_def = 26.3–26.8 for γ = 1 … 8 [E].
  - The "S-determined D-lag-0" filter rejects nothing in practice (34/34 passed).
- **Ordering S-solutions by D-freedom does not pay [M]** (`theory/q523_dorder.py`, 42 real S-solutions).
  - The D-freedom hd = log₂ #D-configurations is 35.8 ± 0.76 bits.
  - It predicts the D-candidates well (correlation 0.95 with log₂ #candidates).
  - But the top half by hd holds 69% of the candidates and takes 59% of the D-check time, a gain of at most 1.17×.

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

- A core needs 2^{D_def} S-solutions [E, calibrated in §3.4]. Each one costs:
  - **S-step:** 2^{F + 1 − n} pair visits (§3.3; F ≈ 54.4–54.5 and n ≈ 37–38 with K = 32 [M]);
  - **D-check:** 1.5 ms on an idle machine, 1.8 ms with four workers sharing it [M] (§3.4, §6).
- **523, y of §5 (D_def = 26.14, so 2^26.14 ≈ 7.4·10^7 S-solutions per core):**
  - S-step: 2^{26.14 + 18.1} ≈ 2^44.2 pair visits;
  - measured join speed: 22–30 ns per pair visit (2^37.4–2^38.8 visits per K = 32 join, four workers) [M];
  - the K = 32 joins find 1.55–1.75 bits more S-solutions than the Gaussian predicts (§6), so a found S-solution
    costs 2.4–3.0 ms of join [M], not the 7 ms of the prediction;
  - **realised total with `mjoin_rad`: 4.2–4.8 ms per S-solution, about 90 core-hours per expected core** [M rates ×
    E count];
  - with `mjoin_rad2` (from 23:24 UTC on 2026-10-11) the production joins run 1.34–1.79× faster (residue classes
    per second on the same three centres, before the switch and over the rest of each join after it), so about
    1.8 ms of join + 1.8 ms of D-check: **about 3.6 ms per S-solution, 75 core-hours per expected core** [M rates ×
    E count];
  - smaller lists are worse: K = 24 at logN 35 cost 8.7 ms of join + 1.9 ms of D-check per S-solution (p523h, §6).
- **Ideal Wagner cost:** 2^{(F + D_def)/2 + 1} ≈ 2^40.7 pair visits. It needs one middle value with about 5000
  parent pairs per side, and the Z₃ middle level only offers about 500 (§3.3). That is the gap between this route's
  2^44 and its ideal.
- **Floor of the route** with these parents: 2^{D_def} × (2^{F + 1 − n_max} + D-check). D_def does not respond to
  list reweighting (§3.4) and moves by at most about 0.5 bit with y (§5).

## 4. Controls [M]

- **Planted, full 523 size** (`theory/q523_planted.py`):
  - a random configuration is planted, the target becomes its P, and the planted rows are forced into the lists;
  - the middle value is its left Z₃ sum, joined with other parent pairs of that middle value;
  - y0: **3/3**, and the y of §5: **2/2**. In every trial the join returns the planted S-solution, the D-check
    returns the planted D, and the assembly reproduces the planted target exactly.
- **Kernel and D-check equivalence:**
  - the three S-join kernels return identical hit sets (v = 235: 6138 and 178 hits; 523: planted and real runs);
  - the three D-check implementations agree on 24 planted and 122 real S-solutions (§3.4).
- **H(316), v = 79** (Z₆ × quartic Z₁₃, F = 6):
  - stage 1 by CP-SAT (0.5 s), stage 2 by the S/D solver (`theory/q523_sdsolve.py`);
  - core `data/core_v79_b1_q523_sdsolve.txt`;
  - bordered GS verified by `theory/q523_finish.py` (independent finisher) and by `theory/paley_b1pipe.py
    verifycyclic`: H Hᵀ = 316·I. SHA-256 of the packed rows: ccef126a….
- **H(940), v = 235** (Z₁₈ × quartic Z₁₃: the 523 shape with q = 13), found by the production pipeline
  (`theory/q523_prod.py`, the same code as the 523 run) [M]:
  - y `logs/q523/y_v235_a.json`; 3 workers; lists K = 24; one middle value each at logN 35;
  - S-joins: 6.36M, 6.01M and 6.05M S-solutions (2^22.6 each; predicted 2^21.0) in about 1450 s each;
  - D-checks: about 0.5 ms per S-solution;
  - **the first core appeared after about 1.07·10^7 S-solutions checked over the three workers.** Model:
    2^{D_def} = 2^23.32 = 1.04·10^7 per core. This is one Poisson event, so it is consistent, not a precise
    calibration. The whole round checked 15.6M S-solutions: 6.01M and 6.05M on workers 1 and 2 with no core, and
    about 3.56M on worker 0 up to its core. That is 1 core against 1.5 expected. There were 60.5 D₃-solutions per
    S-solution;
  - core `data/core_v235_b1_q523_v235e_w0.txt`; bordered completion `data/H940_v235_q523.npz`;
  - verified by `theory/q523_finish.py` (H Hᵀ = 940·I exactly, SHA-256 of the packed rows c9e84f66…), by
    `theory/paley_b1pipe.py verifycyclic` (`logs/q523/finish_v235.json`), and by `scripts/verify_h2092.py`, which
    shares no code with `theory/` (pattern of `scripts/verify_h2292.py`; same SHA-256). The same script will verify
    the 523 core.
  - Three earlier attempts were lost to bugs, all now fixed:
    - v235a and v235b: the D-check buffer bug of §3.4;
    - v235d: the OOM killer, during the exact re-check of about 2^21 hits per centre. That re-check built Python
      objects of about 2 KB per hit; it is now vectorised and stores 36 bytes per S-solution. All processes of the
      session share a 14.3 GB memory cgroup.
  - The driver now runs independent worker processes, so a killed worker no longer hangs a multiprocessing pool.
  - Model: D_def = 23.3, so about 2^23.3 S-solutions per core. The Wagner join at logN 35 predicts about 2^21 per
    middle value.

## 5. Choosing y (Q3) [E, M]

### 5.1 Stage 1 is abundant, but structured

- **4-list join for stage 1** (`theory/q523_ymitm.py`):
  - random block lists with the lag-0 composition and the rows built in;
  - a middle value on a few lags, a 2+2 hash join on the rest, every hit verified exactly;
  - each y is stored in a canonical form (per-block translation and reversal), so duplicates can be dropped.
- **v = 235:** about 10^5 distinct stage-1 solutions per middle value (N = 2^20 per block) [M].
- **523:** composition (6×13, 4×15) gave 4 new y in about 2 min, and (10×13, 1×15) gave 11 y in 5 min [M].
  - Every y in `logs/q523/y523/` is verified exactly.
  - A naive join fails for the other compositions. Two exact structures explain why [M]:
    1. **Residues mod 8.** Every block PAF value is 2 mod 4, so the target needs the four blocks' mod-8 class vectors
       to XOR to all-ones. For some compositions a block has very few patterns (68 mod-16 patterns for a block with a
       single 13).
    2. **Exact PSD values.** At j = 9, 6, 3 the PSDs are integers (A², with A the alternating sum, and the norms of
       the order-3 and order-6 class sums), and Σ_k PSD_k(j) = 2092 must hold exactly. One failed join had every
       left pair at A_A² + A_B² = 4 while the right pairs needed 1860.
  - Simulated annealing on the lag residuals stalls at E = 16 (one lag off by 4) for the same reason
    (`theory/q523_ysa.py`).

### 5.2 What a y changes [M, E]

- **Within a composition, only the pairing of letters across (g, g + 9) matters.**
  - C(D|S) agrees to ±0.05 bit for all y of a composition (16 y screened, `theory/q523_yscreen.py`).
  - h_D is exact from the pair types. H(D|S) per pair (y(g), y(g+9)):

| pair | H(D \| S) bits |
|---|---|
| (±1, ±1) | 1.097 |
| (±1, ±13 or ±15) | 0.792 |
| (13 or 15, same sign) | 0.75 |
| (13 or 15, opposite signs) | 0.5 |
| any dead letter (±27, ±29) | 0 |

  - Example, composition (10×13, 1×15): h_D is 36.16 or 36.44, so D_def is 26.44 or 26.14.
  - The lag-9 equation (Σ_k Σ_{g<9} y_k(g) y_k(g+9) = −58) limits the number of same-sign big–big pairs. The best
    reachable D_def in this composition is about 25.9 [E].
- **Across compositions [E]** (`theory/q523_ycomp.py`):
  - method: covariance-only Gaussian (Mahalanobis term dropped) on random letter placements; it reproduces C(D|S)
    of real y to ±0.07 bit;
  - result: D_def 25.7–26.5 for all six compositions tested, including dead-letter ones. The best estimate is
    (1×15, 2×29) at about 25.7.
- **Conclusion (Q3).**
  - The y choice moves the cost by at most about 0.5 bit (1.4×). The screen is cheap (5 s per y), but it is not a
    lever.
  - Used here: `logs/q523/y523_best.json` (10×13, 1×15, h_D = 36.44).
    - C₂ 146.2, C_mid 51.2, D_def 26.14, so 2^26.14 S-solutions per core;
    - exact lattice rank 35 with the target inside;
    - planted full-size test 2/2.

## 6. The 523 run

- **y:** `logs/q523/y523_best.json` (§5). D_def = 26.14 (model), so a core needs about 2^26.14 ≈ 7.4·10^7
  S-solutions [E]; measured on 1000 production S-solutions, 2^26.02 ≈ 6.8·10^7 (§3.4).
- **Runs** (`logs/q523/prod/`, all on the same y; the live logs are copied to `logs/q523/snap/` at each status update;
  "excess" is log₂ found − log₂ predicted per middle value):

| tag | lists | middle values | S-solutions found | excess (bits) | ms per S-solution: join + D-check | cores |
|---|---|---|---|---|---|---|
| p523a | K = 24, logN 35, first kernel | stopped for the kernel upgrade | – | – | – | – |
| p523b | full row fibres | stopped: F = 61.3 (predicted 2^11.2 per centre) | – | – | – | – |
| p523c (D-checks finished by p523d) | K = 24, logN 35 | 1 | 114,471 | 1.0 | 9.6 + 1.8 | 0 |
| p523e (D-checks finished by p523g) | K = 32, logN 38, pass 2^24 | 2 | 5,815,230 | 1.64–1.69 | 2.6 + 1.7 | 0 |
| p523f, 3 workers | K = 32, logN 38, pass 2^23 | 14 + running | 46,036,472 | 1.55–1.74 | 2.6 + 1.8 | 0 |
| p523g | K = 32, logN 38 | 1 | 948,056 | 1.61 | 4.2 + 1.8 | 0 |
| p523h | K = 24, logN 35 | 26 | 3,522,576 | 1.03–1.22 | 8.7 + 1.9 | 0 |
| p523i | K = 32, logN 38 (replaces p523h) | running | – | – | – | – |

- **Measured rates [M], 4 workers on 4 cores:**
  - K = 32 joins visit 2^37.4–2^38.8 pairs in 1.1–4.0 hours (22–30 ns per pair visit) and find 2^19.9–2^22.7
    S-solutions each, 1.55–1.75 bits above the Gaussian prediction;
  - D-check: 1.7–1.9 ms per S-solution (1.5 ms on an idle machine, §3.4);
  - so 4.2–4.8 ms per S-solution in all at K = 32, about 90 core-hours per expected core; K = 24 lists cost 10.6 ms;
  - since 23:24 UTC the join kernel `mjoin_rad2` (§3.3) runs 1.34–1.79× faster in production (three centres joined
    partly with each kernel: 1.59 → 2.18, 1.89 → 3.39 and 1.24 → 1.66 residue classes per second): about 3.6 ms per
    S-solution, 75 core-hours per expected core. (The `secs` of an S-join record count only the part of the join run
    by the process that finished it.)
- **Memory fix and restart (2026-10-11 22:50 UTC).**
  - A K = 32 worker held 3.3–3.7 GB, so only three fitted into the 14.3 GB cgroup and the fourth ran K = 24 lists
    (p523h), at 2.4× the cost per S-solution.
  - About 2.5 GB of it was the int64 S-signatures of one block's fibres (about 300 bytes per configuration). They
    were kept in a list while the block's join keys were built, and the last block's list stayed alive through the
    join and the D-checks. The fibres are now consumed one at a time.
  - The new code gives identical S-solutions on a v = 235 test (13,323 of 13,323, same digest). A K = 32 worker now
    needs about 1.5 GB during its join, and all four workers run K = 32 (p523f resumed, p523i replaces p523h).
  - A restarted worker now resumes its own interrupted D-check whatever the age of its progress file, then continues
    after the last centre it joined. This was tested by killing and resuming a v = 235 run in mid D-check.
  - The tally (`logs/q523/tally_523.py`) now also counts the checks a killed process made before its D-check was
    resumed: 0.92M checks, previously missed.
- **Status (2026-10-12 05:35 UTC; tally `python3 logs/q523/tally_523.py`):**
  - 51 middle values have been joined, giving 80.0M S-solutions;
  - **75.5M are D-checked (2^26.17): 1.02 expected cores by the model, 1.11 by the measured D-success. There is no
    core yet.**
  - With no core, the probability of seeing none so far is e^−1.11 = 0.33.
  - With four K = 32 workers and `mjoin_rad2` the run checks about 4M S-solutions per hour, about 17–19 hours per
    expected core. The run continues.

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

- **A larger middle level.** The S-step costs 2^{F + 1 − n} per S-solution. With K = 32 lists n ≈ 37–38, because a
  Z₃ middle value has only a few hundred parent pairs (§3.3).
  - Lazy parents (up to 2^24 per block, `theory/q523_wagner2.py`) give no more pairs than 128 row sums per row
    (182–426 against 321–461): the extra parents have atypical Z₃ signatures.
  - Larger K raises F almost as fast as n (§3.3).
  - A middle level with n ≈ 40 at F ≈ 54 would bring the S-step to about 2^15 pair visits per S-solution, well
    below the D-check. A core would then cost about 2^26.1 D-checks ≈ 37 core-hours instead of about 90.
- **The D-deficit.** It is pinned by the S-level energy equation (§3.4) and moves by at most about 0.5 bit with y
  (§5). The only further lever found is the pairing of big letters across (g, g + 9). It is capped by the lag-9
  equation and is worth about 0.25 bit.
- **The D₃ join.** It is the largest part of the 1.5 ms D-check (§3.4). A faster join helps only once the S-step is
  cheaper.
- **Luck.** A run finds cores as a Poisson process. At the expected cost the chance of having one is 63%; at twice
  the cost, 86%.

## 10. Reproduce

```
python3 theory/q523_algebra.py                                   # section 1 checks (eta basis, P vs FFT, v = 41 core, alphabets)
python3 theory/q523_levels.py 18 29 4 logs/cyclo/stage1_cpsat_523.json        # level budget (section 2.2)
python3 theory/q523_lattice.py 18 29 4 logs/q523/y523_best.json               # exact lattice (section 2.3)
python3 theory/q523_chain.py 18 29 4 logs/cyclo/stage1_cpsat_523.json         # per-parent chains (section 3.1, 8)
python3 theory/q523_calib.py 18 29 4 logs/cyclo/stage1_cpsat_523.json --K 24 --beta 1.5 --seed 1   # S-step calibration
python3 theory/q523_dfine_calib.py 18 29 4 logs/cyclo/stage1_cpsat_523.json --beta 1.5 --logN 28   # D-success per S-solution
python3 theory/q523_dcheck_nb.py test; python3 theory/q523_dcheck2.py test      # D-check implementations agree (section 3.4)
python3 theory/q523_dorder.py 18 29 4 logs/cyclo/stage1_cpsat_523.json --logN 29.5 --seed 12   # D-freedom ordering
python3 theory/q523_planted.py 18 29 4 logs/q523/y523_best.json --trials 2 --logN 28           # planted, full size
python3 theory/q523_ymitm.py 18 29 4 --comp 13:10,15:1,1:61 --N 20 --secs 300 --out DIR       # stage-1 solutions
python3 theory/q523_yscreen.py 18 29 4 logs/q523/y523/*.json                   # y-screen (section 5)
python3 theory/q523_ycomp.py 18 29 4 --real logs/q523/y523_best.json --comps "15:1,29:2,1:69" # composition screen
python3 theory/q523_stage1.py 18 41 4 --list-comps                              # 739: no lag-0 composition (section 7)
python3 theory/q523_prod.py 18 13 4 logs/q523/y_v235_a.json --tag v235e --workers 3 --logN 35 --K 24 --seed 600 --resume   # H(940)
python3 theory/q523_prod.py 18 29 4 logs/q523/y523_best.json --tag p523f --workers 3 --hours 400 --logN 38 --K 32 \
    --rows 64 --beta 1 --seed 6000 --pass 23 --resume                         # the 523 run (and p523i: --workers 1 --seed 9000)
python3 logs/q523/tally_523.py                                                  # S-solutions produced / D-checked / expected cores
python3 theory/q523_joinbench.py --classes 64 --p0 1000                         # mjoin_rad vs mjoin_rad2 on one full-size centre
python3 theory/q523_finish.py data/core_v79_b1_q523_sdsolve.txt                # completion + exact H H^T = N I
python3 scripts/verify_h2092.py data/core_v235_b1_q523_v235e_w0.txt OUT.npz  # independent check (H(940); the 523 core likewise)
```
