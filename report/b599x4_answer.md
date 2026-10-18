# 599, c = 6, ×4 on Z₉ has no core: a coset-Parseval obstruction, the rigidity of the Z₃² variants, and what a stage-2 solver can still buy

Date: 2026-10-18. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_599x4_stage2.md`.
Scripts: `theory/b599x4_*.py`, `theory/m2_*.py`. Logs: `logs/b599x4/`, `logs/m2/`.

**Tags.** [P] proved here (proof in the text); [M] exact computation (script named); [E] estimate (estimator and its
calibration named).

---

## 0. Bottom line

1. **(A) No H(2396), and none can come from this slice: 599, c = 6, q = 11, s = 2, H′ = Z₉, U = ×4 on Z₉ has no
   core [P, M]** (§1, `theory/b599x4_parseval.py`).
   - ×4 fixes 3Z₉ pointwise, and its orbits {1, 4, 7} = 1 + 3Z₉ and {2, 5, 8} = 2 + 3Z₉ are cosets of 3Z₉. So every
     character of order 9 on Z₉ sums to 0 over them, and the 396 characters of G of that kind see only x on
     C = Q × 3Z₉ × Z₁₁ (198 points).
   - Parseval then caps their energy at 3·198²·(8/9) = **104,544 per block**. The core equations need 4v = 2396 at each
     of them: **237,204 per block on average** (948,816 over the four blocks, against at most 418,176).
   - The bound does not involve y: every stage-1 solution of this slice fails stage 2.
   - Checked three ways: an exact cyclotomic identity check over all orbits, 2,000 random invariant x (the maximum
     seen is exactly the bound), and a brute-force control at v = 19 (G = Z₂ × Z₉, b = 1). There, ×4 and ×2 on Z₉
     give 0 cores, as the criterion predicts, against 18,432 for −1 on Z₉ with the same 10 orbits.

2. **Why the estimates promised 2^15–2^19 cores per class [E, M]** (§1.3, §4).
   - The stage-2 budgets of `report/bge3_stage1_answer.md` §4 and §6 come from `theory/ab_typical_levels.py`. It reads
     every lattice Gaussian at the *sample mean* (`ab_gauss.gauss_bits_at_mean`), so it charges each equation as if
     its target sat at the centre of the lift distribution. Here the target is outside the reachable set, and nothing
     in that reading can show it.
   - The symptom was visible: on typical y's the compression along 3Z₉ looks extremely constraint-dense (C_M ≈ 48 of
     C_S ≈ 65.6 bits, dense by +17; nested-level formula 2^17.4 per S-solution against 2^26.4 plain;
     `logs/b599x4/levelprobe_seed*.log`). Its complement, the order-9 level, "costs" only 18 bits at the mean,
     because it sees just 18 values per block; its true target is unreachable.
   - **Calibration of the at-mean reading on the two calibration slices** (real y's, same lift samples, true target
     against mean; `theory/m2_levelcal.py`):

     | slice | C₂ true / at mean | C_S true / at mean | D_def true / at mean | Mahalanobis |
     |---|---|---|---|---|
     | 523 (production y) | 146.15 / 146.09 | 83.59 / 83.52 | **26.13 / 26.14** | 0.36 |
     | 347 (published y) | 117.92 / 111.95 | 69.24 / 65.51 | **16.14 / 13.90** | 2.88 |

     The model's "calibration at 347 and 523" used the true-target values on real y's (16.2 and 26.1). Read at the
     mean, the same y's give D_def 2.2 bits too small at 347 and C₂ 6 bits too small. So every at-mean budget of a
     slice without a real y is a lower bound on its cost, by about 2 bits on a benign slice and without limit on a
     rigid one.

3. **The obstruction is general, and a screen of every generated slice family applies it [P, M]** (§2,
   `theory/b599x4_cosetscreen.py`).
   - **Criterion.** Let N ≤ H be a subgroup of prime order and C_N the union of the slice orbits O with O + N ≠ O.
     The slice has no core when |C_N|(1 + 1/|N|) < v (|N| odd), or |C_N| < (1 − 1/|N|)·v (|N| even).
   - It bites whenever the slice group contains an element that acts on an axis Z_{p^e} of H (e ≥ 2) as some
     u ≡ 1 (mod p^{e−1}), u ≠ 1, and as the identity on every other axis: ×4 or ×7 on Z₉, −1 on Z₄, ×5 on Z₈, ×17
     on Z₃₂, and any group containing one of them. A multiplier coupled to a nontrivial action on another axis does
     not count. On the b = 1 sweeps this description matches the screen on all 116 slices.
   - Killed:

     | family | killed / screened |
     |---|---|
     | c = 4, 6, 8 batch configurations | 73 / 751 |
     | `ab_sweep` b ≥ 3 configurations | 152 / 1535 |
     | b = 0 configurations | 0 / 66 |
     | b = 1 multiplier slices (`logs/msweep_record/`) | 34 / 116 |
     | the 112 skipped configurations of `report/bge3_stage1_answer.md` | 24 / 112 |

   - The slices of known cores (347, 523 ×181, 647 ×353) pass.
   - Every killed slice of the 112 was already rated over budget or dead in stage 1. Only one had been ranked rich and
     closest to the budget: 599 ×4 at s = 2 (its s = 5 twin was at 2^159).
   - Also killed: 659 c = 4 with −1 on H′ = Z₄ (the −1 variant of the "best b ≥ 3 survivor" of
     `report/alpoge_design_answer.md` §11.2, margin 1.005), and 537's ×149 (b = 1, whole-slice estimate −0.6).

4. **The brief's fallbacks (its step 6) are no better: the Z₃² variants carry a stage-2 rigidity of at least 9.5 bits
   per y [P identity, M tables, E cost]** (§3, `theory/b599x4_z3a_rigidity.py`).
   - **The slices.** H′ = Z₃² with −1 on one axis ("Z₃²-a", 2^56.5–2^61.2 per core in the model), and −1 on Q's Z₃
     (2^61.3–2^62.0 per core).
   - **The rigid sub-problem.** The characters nontrivial on the −1 axis see only d = x(·, 0) − x(·, 1). The lag 0 of
     that sub-problem, minus its stage-1 compression, forces an exact count: D₁₁ = D*(y) differing free-sign pairs
     among the N₁₁ pairs of small letters (D* + 6 for Q's Z₃). Under uniform lifts D₁₁ ~ Bin(N₁₁, ½).
   - **The bound.** An exact DP over every letter-pair composition allowed by the two stage-1 energies (lag 0 = 2176
     and the d-energy 599, resp. 533) shows that the best possible stage-1 solution pays **9.82 bits (Z₃²-a) and
     9.50 bits (−1 on Q's Z₃) more than the at-mean model charges.** The median over compositions is 18.4 bits.
   - **Cost.** That is ≥ 2^66.3 and ≥ 2^70.8 per core at the model's solver cost. Correlations with the other
     equations cannot lower it: in the Gaussian model the joint Mahalanobis distance is at least the marginal one.
   - Reaching 2^47 would need an S-solver about 19 bits cheaper than the model's. At 523 and 347 no nested level
     saves anything (§5.1).
   - **A different U on Z₉ does not help:** ⟨2⟩ contains ×4 (no core [P]), −1 is the old lead (stage 1 2^−6 per class
     [E]).
   - **The non-split choice H′ = 3Z₂₇ admits only U ∈ {1, −1}.** Every other subgroup of Z₂₇* contains ×10 ≡ 1
     (mod 9).
   - **c = 10, Paley(59): stage 1 is empty by energy [P].** T1(0) = 272 is reached neither by ±1 letters (40) nor
     with any letter ±57, ±59 (each adds at least 3,248).

5. **(B) The multi-level stage-2 solver: no slice in reach needs one, and its calibration slices show why [M, E]**
   (§5).
   - **523.** The only sub-level of the Z₉ S-level is the Z₃ compression already used as the Wagner middle level. It
     is not constraint-dense: C_M − H_M/2 = 32.41 − 78.37/2 = −6.8 (exact H_M, true target). Its per-M-solution
     variant is the per-parent design of `report/paley_523_answer.md` §3.1, 2^18.7 per S-solution against the Wagner
     tree's 2^18.1.
   - **347.** None of the four order-3 sub-levels is dense (−8.8 to −9.9). With the old numpy D-check the S-join was
     about 3% of the run; with the compiled one it is about half (4.5 of 10 ms per S-solution).
   - **599 ×4.** The only dense level is the artefact of item 2.
   - So no 599 production solver was built. What §5 adds is a compiled D-check for the b ≥ 3 pipeline and a
     multi-core run of the 347 control that turns its one-event calibration into a measured S-solutions-per-core
     rate (in progress at the time of this draft).

6. **(C) The floor.** For ×4 there is none: the slice is empty at any compute. For the Z₃² variants it is about 2^66
   per core at the best y.
   - No GP slice of 599 at c = 6 (q = 11 is forced: h = 99 = 9·11) is within the budget.
   - What could still reach H(2396):
     - a slice family outside single-axis small-orbit multipliers (most of those are killed or rigid);
     - a c = 10 or c = 12 border with a flat corner (c = 10 then needs s = 29, a stage 2 of 30 parts per orbit);
     - b = 0 over Z₅₉₉ with the flat-pool k-tree, between the 491 and 859 rungs of `report/b0_flatpool_ladder.md`
       (beyond this machine there).

---

## 1. The obstruction [P, M]

### 1.1 Theorem

Let S act on G by automorphisms (here per-axis multipliers), N ≤ H a nontrivial subgroup, and
C = ⋃{O : O an S-orbit, O + N ≠ O}. Let x: G → {±1} be constant on S-orbits.
- (1) For every character χ with χ|_N ≠ 1: x̂(χ) = f̂(χ), f = x·1_C. *Proof.* An orbit O outside C is a union of
  N-cosets, x is constant on it, and Σ_{n ∈ N} χ(g₀ + n) = χ(g₀)·Σ_{n∈N} χ(n) = 0.
- (2) Σ_{χ|_N ≠ 1} |f̂(χ)|² = |G|·|C| − |G/N|·Σ_{c ∈ G/N} f̄(c)², where f̄(c) is the sum of f over the N-coset c
  (Parseval on G, and on G/N for the compression of f, whose characters are those trivial on N). C is a union of
  N-cosets; for |N| odd each coset sum inside C is odd, so the sum is at most |G|·|C|·(1 − 1/|N|²). For |N| even it
  is at most |G|·|C|.
- (3) Every χ with χ|_N ≠ 1 is nontrivial on H (N ≤ H), so a core needs Σ_k |x̂_k(χ)|² = 4v at each of the
  |G|(1 − 1/|N|) of them.

Hence no core exists when 4|C|(1 − 1/|N|²) < (1 − 1/|N|)·4v, i.e. |C|(1 + 1/|N|) < v (|N| odd), or |C| < (1 − 1/|N|)v
(|N| even). For b = 0 and b = 1 take H = G (every nontrivial character needs 4v). ∎

### 1.2 The 599 ×4 slice

- G = Z₂ × Z₃ × Z₉ × Z₁₁, S = ⟨×4 on Z₉⟩ × ⟨squares on Z₁₁⟩, 90 orbits (sizes 1, 3, 5, 15).
- N = 3Z₉: the orbits through h ∈ 1 + 3Z₉ or 2 + 3Z₉ are 3Z₉-stable; those through h ∈ 3Z₉ are not.
  C = Q × 3Z₉ × Z₁₁, |C| = 198, and 198·(4/3) = 264 < 599.
- In numbers: 396 characters of order 9 on Z₉, at most 104,544 per block among them, 948,816 needed in all.

`theory/b599x4_parseval.py` (numpy only, nothing imported from `theory/`; log `logs/b599x4/parseval_599x4.log`):
- builds the orbits from the generators;
- checks the coset structure;
- checks the character-sum identity exactly in Z[ζ₉] (coefficient vectors reduced modulo Φ₉);
- evaluates the bound;
- confirms it on 2,000 random invariant x (maximum 104,544, attained);
- runs the v = 19 brute-force control (exhaustive 2+2 join over all 2^10 blocks).

### 1.3 How the at-mean budget saw this slice

`theory/b599x4_levelprobe.py` (typical y's, at the mean, five seeds) gives:

| level | H_M | C_M | C_M − H_M/2 |
|---|---|---|---|
| S + 3Z₉ (onto Z₃ × Z₃) | 60.4–60.5 | 47.9–48.0 | +17.4 to +17.7 |
| S + order-3 subgroups of Z₃ × 3Z₉ | 43.4–43.5 | 25.6–26.1 | +3.9 to +4.3 |
| S + Z₉ | 44.0 | 21.0–21.4 | −0.6 to −1.0 |

The 3Z₉ level looks like the ideal nested level: the nested formula of `report/b599_answer.md` §5 gives 2^17.4 per
S-solution against 2^26.4 plain, about 9 bits cheaper. The formula is right about the counts at the mean and wrong
about the target. The fine level is cheap at the mean because it sees only the values at h ∈ {0, 3, 6}, and for the
same reason its target cannot be met.

---

## 2. The general screen [P, M]

`theory/b599x4_cosetscreen.py` builds every configuration of
- `theory/ab_cbatch.configs` (c = 4, 6, 8 at the 19 open orders);
- `theory/ab_sweep.configs_bge3` (maxU = 2) and `configs_b0`;
- the b = 1 multiplier sweeps `logs/msweep_record/m*.jsonl`;
- the 112 configurations of `logs/s1z/table.json`.

For every subgroup N ≤ H of prime order it computes C_N exactly and the margin v/(|C_N|(1 + 1/|N|)) (dead when > 1).
Output: `logs/b599x4/cosetscreen.jsonl`, summary `logs/b599x4/cosetscreen_summary.json`.

- Known-core slices, which must pass: 347 (margin 0.96), 523 ×181 (0.97), 647 ×353 (0.95).
- Killed among the b ≥ 3 families, by order: 599: 20, 643: 50, 659: 18, 739: 12, 839: 56, 883: 23, 917: 40, 947: 6.
- Killed b = 1 slices: 11 at 523, 9 at 537, 14 at 739. On all 116 b = 1 slices the kills are exactly the
  multiplier groups containing an element that is ≢ 1 but ≡ 1 mod 3 on Z₉ (mod 4 on Z₈), and ≡ 1 on the other CRT
  components. At 537, five groups contain ×5 on Z₈ only coupled with a nontrivial action on Z₆₇, and they pass. The
  523 production slice ×181 is ≡ 1 on Z₁₈ and passes.
- The 24 killed among the 112 skipped configurations:
  - 599 ×4 at s = 2 and 5;
  - 14 slices of 643 at c = 4 with multipliers on Z₃₂, Z₁₆, Z₈ or Z₄;
  - 8 slices of 839 at c = 8 with ×5 on Z₈ or −1 on Z₄.

  Only the 599 ×4 pair had been ranked as rich and near budget.

---

## 3. The Z₃² variants: stage-2 rigidity [P, M, E]

### 3.1 The identity [P]

In the slice Z₃²-a write L = Q × Z₃′ × Z₁₁ (198 points), r = x(·, 0, ·), s = x(·, 1, ·) = x(·, 2, ·) (the −1 axis is Z₃″),
and d = r − s ∈ {0, ±2}.
- For φ″ ≠ 1 on Z₃″, x̂(χ_L, φ″) = d̂(χ_L), since φ″(1) + φ″(2) = −1. These characters are nontrivial on H, so
  Σ_k PAF_L(d_k) = 4v·δ: at lag 0, W = Σ (d/2)² = v = 599.
- Compressed along Z₁₁ (stage 1): Σ (δ/2)² = v, with δ = y(·, 0) − y(·, 1) on L′ = Q × Z₃′.

Per fibre, with d/2 = (e₀ at z = 0, e₁ on the quadratic residues, e₂ on the non-residues):

    (δ/2)² − Σ_z (d/2)² = 20e₁² + 20e₂² + 10e₀(e₁ + e₂) + 50e₁e₂ =: cross.

By the letters (y = α + 5(p + n)):
- for a pair of small letters (±1, ±1), cross = −10 when the two free signs differ and 0 otherwise;
- for every other pair it is fixed by y: 10–30 for (small, large) and 70–110 for large letters of opposite signs
  (`logs/b599x4/z3a_rigidity.json`, `fixed_cross_values`).

So a core needs **D₁₁ = D*(y)**, D₁₁ the number of small-letter pairs whose free signs differ and D* one tenth of the
fixed part. With −1 on Q's Z₃ the same split gives W = (196v + 2b)/198 = 593 and Σ(δ/2)² = (16v + 2b)/18 = 533,
hence D₁₁ = D* + 6.

### 3.2 How far D* sits from the lifts [M, E]

`theory/b599x4_z3a_rigidity.py` does three things:
- tabulates cross per letter pair by brute force over the fibre patterns;
- checks the identity on 200 random lifts of random y;
- runs an exact DP over every multiset of letter pairs with the stage-1 lag-0 energy Σ(y_r² + 2y_s²) = 2176 and the
  stage-1 d-energy (599, resp. 533). These are both necessary for every stage-1 solution, so the DP is a relaxation and
  its minimum is a lower bound.

| slice | compositions | smallest excess over the mode | at | median (over compositions) |
|---|---|---|---|---|
| Z₃²-a | 1,725 | **9.82 bits** | N₁₁ = 63, D* = 46, P = 2^−13.15 | 18.4 |
| −1 on Q's Z₃ | 2,206 | **9.50 bits** | N₁₁ = 65, D₁₁ = 47, P = 2^−12.85 | 18.5 |

- The excess is −log₂ P(Bin(N₁₁, ½) = D*) minus the modal value, which is what an at-mean Gaussian charges for this
  equation.
- The best cases need almost every small pair to have opposite α, which is a special corner of stage 1. Typical
  stage-1 solutions pay much more.
- Reading the excess as a cost is [E]; the identity is [P].
  - It treats one equation of the stage-2 system on its own.
  - In the Gaussian model the squared Mahalanobis distance of the whole stage-2 vector is at least that of any single
    standardised coordinate, so the other equations can only add to it.

---

## 4. The at-mean model on real y's [E, M]

`theory/m2_levelcal.py` draws 6,000 uniform lifts of a real stage-1 solution and reads every lattice Gaussian twice,
at the true target (`q523_gauss.gauss_bits`, Mahalanobis included) and at the sample mean
(`ab_gauss.gauss_bits_at_mean`). Logs: `logs/m2/levelcal_523.json`, `logs/m2/levelcal_347.json`.

- **523** (`logs/q523/y523_best.json`, CRT to Z₂ × Z₉):
  - the two readings agree to 0.1 bit;
  - the record's figures are reproduced: C₂ 146.15 (record 146.2), D_def 26.13 (record 26.14);
  - the Z₃ sub-level costs C_M = 32.41.
- **347** (the published core's y):
  - true C₂ 117.92 against 111.95 at the mean, D_def 16.14 against 13.90, Mahalanobis 2.88;
  - the four order-3 sub-levels of G″ = Z₃² cost 36.4–37.2 bits (true) against 33.0–33.7 (at the mean).

So the model's two calibration points are true-target readings. The at-mean reading used for every slice without a
real y is optimistic by about 2 bits of D_def on a benign slice (347), and by ∞ (×4) or ≥ 9.5 bits (the Z₃² variants)
on a rigid one.

---

## 5. The stage-2 solver (B)

### 5.1 Where a nested level could help, measured on the real y's [M entropies, E bits]

The nested-level cost of `report/b599_answer.md` §5 per S-solution is
2^(C_S − H_S + H_M/2) + 2^((C_S − C_M) − (H_S − H_M)/2). It saves work only when the level is constraint-dense,
C_M > H_M/2.
- `theory/m2_levelcal.py` now also computes H_M exactly: the Shannon entropy of the compressed S-values under uniform
  lifts, by a convolution over the orbits of each coset orbit.
- At the S level it returns H₂ − h_D exactly (131.88 at 347, 143.25 at 523), an independent check of `sd_freedom`.
- Everything is at the true target.

| slice | level | H_M | C_M | C_M − H_M/2 | nested 2^· per S-solution | realised now |
|---|---|---|---|---|---|---|
| 523 | S + 3Z₉ (the Z₃ middle level) | 78.37 | 32.41 | −6.8 | 2^18.7 | 2^18.1 pair visits (Wagner, K = 32) |
| 347 | S + each order-3 subgroup of Z₃² | 91.9–92.6 | 36.4–37.2 | −8.8 to −9.9 | 2^12.1–2^13.2 | 2^13.1 pair visits (Wagner, this run) |
| 599 ×4 | S + 3Z₉ | 60.4–60.5 | 47.9–48.0 (at the mean) | +17.4 to +17.7 | 2^17.4 | (no core, §1) |

- **523.** A nested level adds nothing: the only sub-level is the Z₃ middle level the record already uses, and its
  per-M-solution variant (the per-parent design of `report/paley_523_answer.md` §3.1) costs 0.6 bits more than the
  Wagner tree.
- **347.** None of the four sub-levels is dense. Their nested cost equals the realised Wagner cost, and the S-join was
  not the bottleneck anyway (§5.3).
- **599 ×4.** The only dense level in this work is the artefact of the obstruction.

So no production multi-level S-solver was built: it has no slice to run on.

### 5.2 A compiled D-check for the b ≥ 3 pipeline [M]

`theory/m2_dcheck.py` replaces the numpy D-check of `theory/ab_stage2.py` (`_dcheck_fast`). It uses the same
algorithm and returns the same output format.
- Per S-solution, each block's D-configurations are enumerated with per-orbit-pair bilinear tables (numba).
- The 2+2 join on a 64-bit linear hash of the D-level signature uses the block pairing with the least work and an
  open-addressing table on the smaller pair side; every key match is re-checked exactly.

Validation (`logs/m2/dcheck_test347.log`):
- the published 347 core and the a1 core of `logs/ab/prod347/` are returned from their own S-solutions;
- 20 of 20 planted lifts (D-target set to the planted D-signature) are found;
- identical solution sets to the reference on all 600 real S-solutions of a middle value.

Speed: 5.6 ms per S-solution against 54 ms (one core; 1.6 ms of it builds the signatures).

`theory/m2_prod.py` is the restart-safe driver:
- S-join by `ab_stage2.FastSolver2`, then the compiled D-check;
- the S-solutions of each middle value are saved before any D-check, and D-check progress is logged, so a restarted
  worker resumes inside a middle value;
- every D-solution is assembled, checked on G, completed with the 347 template border and its forced corner, and
  checked H Hᵀ = 1388·I;
- `theory/m2_tally347.py` deduplicates S-solutions across workers and cores by hash.

### 5.3 The 347 control with many cores [M counts, E model]

*Run m2a (4 workers) in progress; the tally is filled in when it ends.*

### 5.4 523 [M]

- Model at the true target on the production y (`theory/m2_levelcal.py`): C₂ 146.15, C_S 83.59, D_def 26.13 (record:
  146.2, 26.14).
- D-success on the record's 1,000 saved production S-solutions (`theory/q523_dsuccess_real.py`, re-run, 51 s,
  `logs/m2/dsuccess523_repro.log`): **2^−26.02 per S-solution** (record 2^−26.02). Also reproduced:
  - 2^7.84 D₃-solutions and 2^11.21 D-candidates per S-solution;
  - 37.23 fine bits, Mahalanobis 0.05;
  - this against the model's D_def 26.13 above.

---

## 6. What remains for 599

| route | status |
|---|---|
| c = 6, ×4 or ⟨2⟩ on Z₉ | no core [P] (§1) |
| c = 6, Z₃²-a; −1 on Q's Z₃ (H′ = Z₉ or Z₃²) | ≥ 2^66.3 / ≥ 2^70.8 per core [P identity, E cost] (§3) |
| c = 6, −1 on Z₉ (the old lead); −1 on both axes of Z₃² (neg23) | stage 1 2^−6.0 per class, resp. ≪ 1 distinct solution [E] (`report/bge3_stage1_answer.md` §4.4, `report/b599_answer.md` §0) |
| c = 6, −1 on Q's Z₃ and on H′ (28 orbits) | stage 1 2^−39.5 / 2^−40.1 per class [E] (`report/bge3_stage1_answer.md` §4.2) |
| c = 6, U = 1 | 2^84 per core at the mean [E] (so ≥ 2^84) |
| c = 6, H′ = 3Z₂₇ (non-split) | only U ∈ {1, −1} escape §2; not costed |
| c = 8 | dead [E] (`report/b599_answer.md` §6) |
| c = 10, Paley(59) | stage 1 empty by energy [P] |
| c = 10, s = 29; c = 12 | need a border with a flat corner first (none known) |
| b = 1 (m = 598) | dead [M] (`report/b599_answer.md` §7) |
| b = 0 over Z₅₉₉ | the flat-pool k-tree route of `report/b0_flatpool_ladder.md`; untried, between its 491 and 859 rungs |

---

## 7. Reproduce

```
python3 theory/b599x4_parseval.py --out logs/b599x4/parseval_599x4.json         # §1: proof checks, v = 19 control (1 s)
python3 theory/b599x4_cosetscreen.py logs/b599x4/cosetscreen.jsonl --summary S.json   # §2: every family (about 15 min)
python3 theory/b599x4_z3a_rigidity.py --case z3a  --out FILE                     # §3 (1 s)
python3 theory/b599x4_z3a_rigidity.py --case negQ3 --out FILE
python3 theory/m2_levelcal.py 523 --N 6000 --out FILE                             # §4, §5.1 (a few minutes each)
python3 theory/m2_levelcal.py 347 --N 6000 --out FILE
python3 theory/b599x4_levelprobe.py SEED 4000                                     # §1.3 (typical y, at the mean)
python3 theory/m2_dcheck.py test347 --n 600                                       # §5.2 validation (about 3 min)
for w in 0 1 2 3; do python3 theory/m2_prod.py --tag m2a --worker $w --middles 400 --nsample 8000000 & done   # §5.3
python3 theory/m2_tally347.py m2a                                                 # §5.3 tally
python3 scripts/verify_packed_h.py logs/m2/prod347/H1388_<core>.npz 1388          # independent check of each H(1388)
python3 theory/q523_dsuccess_real.py 18 29 4 logs/q523/y523_best.json logs/q523/ssol_sample_p523f_w2_c7.npy --n 1000 --seed 1   # §5.4
```
