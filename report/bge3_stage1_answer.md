# b ≥ 3 stage 1 for the 117 skipped configurations: calibrated z-conditioned counts and exact stage-2 screens

Date: 2026-10-17. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_p2_bge3_stage1.md`.
Scripts: `theory/s1z_*.py`. Logs: `logs/s1z/`.

**Tags.** [P] proved here (proof in the text); [M] exact computation (script named); [E] estimate (estimator and its
calibration named).

---

## 0. Bottom line

1. **No H(4v) was built.** No configuration among the 117 skipped entries has a stage 1 that is rich per symmetry class
   together with a stage 2 under the 2^47 budget, and the one configuration that came close is proved to have no core.
   So the prompt's sampling and stage-2 steps had no input (§7).

2. **The 117 skipped entries are 112 distinct configurations, 107 slices.** Five 599 c = 6 entries of
   `logs/ab/c6_bordered.jsonl` appear twice (once under the older key without c); five pairs of the 112 are one slice
   each (one U given by two generator lists, or a swap of equal axes of H'). Verdicts of the 112 (§3): stage 1 empty
   [M] 3 (all three also [P]); stage 1 empty [E] 1 (a mod-3 lattice obstruction); stage 1 at most
   2^−4 per symmetry class [E] 19; no core [M] 1; stage 2 over 2^47 per core [E] 88, among them all 85 rich ones; open 0.

3. **Calibration of the z-conditioned count (§2) [M vs E].**
   - The count is now taken over **all** compression quadruples, grouped into multisets of symmetry classes with exact
     multiplicities (`theory/s1z_lib.py`), not over first-found or random quadruples. The exact part (Q-level and
     lag 0) is cross-checked against two independent exact counters and the old canonical counts of the v = 35, 47
     controls.
   - **Totals must be read per symmetry class.** Every ordered stage-1 solution comes with an orbit of block orders,
     signs, U-fixed translations, reversals and automorphisms of order 2^13.6 (347's slice) to 2^46 (643 with U = 1).
     "Rich" in ordered terms can mean less than one solution up to symmetry: the 599 lead's 2^18.5 ordered solutions
     are 2^−6.0 per class (this is `report/b599_answer.md`'s "≪ 1 distinct solution"), the 719 lead's 2^27.3 are
     2^−2.1 per class.
   - Against exact totals on 347's slice (all 27 multisets) and 27 same-shape rungs: where solutions exist the
     estimate is within ±3.6 bits (347: 1.9 low; v = 203: 3.6 low; v = 127: 2.7 high). **Every rung whose exact
     count is 0 was estimated at ≤ 2^−4 solutions per symmetry class** (5 of 5; seven more rungs have no
     compression quadruple at all). Per multiset the estimate is unreliable (347: 12 of 27 multisets exactly empty
     although estimated at 2^2–2^8 each).
   - Estimating from K multisets drawn in proportion to their exact lag-0 weight: at K = 40 the 5–95% range is within
     ±0.4 bit of the all-multiset estimate on 15 of 20 rungs and within ±0.6 on 17; on three exactly empty rungs a
     few multisets carry the weight and the subsample can miss them (down to −7.5 bits; up at most +1.4).
   - Symmetric solutions are a small share where solutions are plentiful: 4.9% of 347's (all blocks fixed by a
     nontrivial affine map), 0% at v = 203 and v = 127.
   - **Margin:** a configuration counts as rich when its estimate is at least 2^+4 per symmetry class, as dead [E]
     when it is at most 2^−4; in between it is undecided.

4. **Stage-2 budgets of all 112 (§4) [E].** Only six configurations cost at most 2^48.1 per core by the model
   calibrated at 347 and 523 (2^{D_def + C_S − H_S/2}):
   - 739 c = 4 (q = 23, s = 2, H' = Z8): 2^29.5, but **stage 1 is empty [M, P]** (screen A; §4.1 gives the
     parity reason: an odd number of large letters cannot be balanced across the cosets);
   - 599 c = 6 with −1 on Q's Z3 and on H' (Z9 or Z3²): 2^34.3 each, but stage 1 is 2^−40 per symmetry class [E];
   - 919 c = 8 (q = 19, s = 2, −1 on Z3; the order is already built): 2^38.7, stage 1 outside the lattice in every
     multiset [E] (a mod-3 congruence, §4.3);
   - the 599 lead (−1 on Z9): 2^43.5, stage 1 2^−6.0 per symmetry class [E] (and the symmetric families searched
     in `report/b599_answer.md` §4 are empty [M]);
   - **719 c = 6, q = 17, s = 2, −1 on Q's Z3 (H' = Z7): 2^48.1, the only one with a non-negligible stage 1.**

5. **The 719 lead has no core [M] (§5).** Its stage 1 is 2^27.3 ordered solutions [E] (2^24.0 in its 18 bordered
   multisets, 13 in 347's template class and 5 monomial/Singer). But stage 2 is decided without finding a single y:
   the exact compressions of stage 2 along H' (Z17 level per character, Z2 × Z17, Z3 × Z17 and Q × Z17) depend on a
   stage-1 solution only through its block **types** (the multiset of letters per orbit class), and **all
   242,977 type quadruples of all 39 compression multisets fail one of them** (5,991 at Z17, 231,936 at Z2 × Z17,
   5,050 at Q × Z17). So no stage-1 solution of this slice lifts to a core
   (`theory/s1z_typescreen.py`, validated on the published 347 core, on planted lifts against an independent exact
   count, and on the v = 47 control, whose Q × Z7 count at its core's y is exactly the 5,184 cores found there earlier).

6. **The rich configurations are all over budget [E] (§6).** The cheapest rich one is 599, c = 6, q = 11, s = 2 with
   U = ×4 on Z9: stage 1 2^+30.8 per symmetry class, stage 2 2^49.9–2^52.6 per core over six typical y's, 2^+15 to
   2^+19 expected cores per class. Within the 599 c = 6, s = 2 family stage-1 richness and the stage-2 deficit rise
   together and step over the window "rich and under 2^47" (from −1 on Z9, 2^−6.0 per class at 2^43.5 per core, to ×4
   on Z9, 2^+30.8 at 2^51.5). The other rich 599 s = 2 slices with U ≠ 1 cost 2^56.5–2^62, 719 q = 7 with ×4 on Z17
   2^71.5, everything else that is rich 2^84 and up (c = 8 and s ≥ 3 2^116 and up, the 643 c = 4 slices 2^167 and up).

7. **Exact stage-2 screens elsewhere [M].** 719, c = 6, q = 17, s = 2 with U = 1 (stage 1 2^−9.8 per class [E], 2^60.7
   per core): 85 of its 89 compression multisets have no lift at Z17 or Z2 × Z17; the other 4 pass Z17, Z2 × Z17 and
   Z3 × Z17. The cheapest rich slice (×4 on Z9) passes Z11 and Z2 × Z11 in all 85 multisets; the joint levels do not
   fit in memory for its letters (§8.2).

8. **Frontier (§8).** Beyond these 112: (i) rich-but-expensive 599 c = 6 slices (2^49.9–2^62 per core) need a
   multi-level stage-2 solver and the joint compression screens, which do not yet fit in memory for their letters;
   (ii) configurations with more than 60 G'-orbits within the open orders' other subgroups; (iii) c ≥ 10, where no
   border with a flat corner is known (669, 789); (iv) non-split H and noncyclic Q as in 359.

---

## 1. Setting and new tools

Setting as in `report/alpoge_design_answer.md` §1: c = b + 1, G = G' × Z_q, H = H' × Z_q, Q = G'/H' (split here:
G' = Q × H'), slice group S = U × M, stage 1 the Z_q compression y on G' with
Σ_k PAF_{G'}(y_k) = T1, stage 0 its Q-compression z with Σ_k PAF_Q(z_k) = T_Q.

| tool | what it does |
|---|---|
| `theory/s1z_lib.py` | `ZOrb`: rows (every achievable row sum in some quadruple with Σ r² = 4v; no truncation), compressions with the Parseval and per-character bounds and their exact energy distributions, classes under the per-block symmetries, class multisets with Σ PAF_Q = T_Q by a 2+2 join and their exact lag-0 counts, Gaussian of the remaining lags (finite-sample log-det correction), direct exact counts; `SDCount`: exact counts by the S/D split; `log2_symmetry`: order of the symmetry group of the stage-1 solution set |
| `theory/s1z_screen.py` | the production z-screen (replaces the never-run `theory/ab_s1_zscreen.py`, whose rows came from screen A's list truncated at 40 quadruples and empty when screen A gives up); restart-safe, `--steal` shares one queue between workers |
| `theory/s1z_calib.py`, `s1z_calib_batch.py`, `s1z_calib_report.py` | calibration against exact counts, rungs, summary |
| `theory/s1z_budget_batch.py` | typical-y stage-2 budget of the 112 configurations, no stage-2 orbit cap |
| `theory/s1z_border.py` | which multisets carry the 347 template border or a monomial/Singer border |
| `theory/s1z_qzq.py`, `s1z_typescreen.py` | exact compression screens of stage 2 per block type, over all type quadruples |
| `theory/s1z_zqscreen.py`, `s1z_char_batch.py` | the same aggregated over all types of a compression quadruple (cheap levels; per Q-character) |
| `theory/s1z_cpsat.py` | CP-SAT stage-1 sampler for a prescribed compression quadruple |
| `theory/s1z_table.py` | the table of §3 |
| `theory/s1z_validate.py` | the validations of §1, §2.1 and §5.3 (logs in `logs/s1z/validate/`) |

**Symmetries [P].** The stage-1 solution set is closed under, per block, y ↦ −y, y ↦ y(· − t) for t ∈ G' fixed by
U, y ↦ y(−·); under permutations of the blocks; and under the automorphisms of G' that preserve H' and commute with
U (per-axis unit scalars; for Q = Z2^r also GL(r, 2) on Q). All of them keep T1, the letters and U-invariance, and act
on compressions by translation, sign, reversal and Aut(Q), keeping PAF_Q and every count below. So a multiset of
compression classes stands for (4!/Π multiplicities!) × Π class sizes ordered compression quadruples, all with the same
exact and estimated counts (checked on 347: 24 random members of four multisets, identical exact counts).

---

## 2. Calibration [M vs E]

(`theory/s1z_calib.py --exact`, `logs/s1z/calib_347.jsonl`, `logs/s1z/calib/calib_*.jsonl`,
summary `logs/s1z/calib/calib_summary.json` by `theory/s1z_calib_report.py`.)

### 2.1 The exact counters

- Direct: per block every completion with the compression (energy-capped) enumerated, then the energy-split hashed 2+2
  join on the full PAF vectors (`theory/ab_s1sd.ejoin4`), every solution re-checked.
- S/D (`SDCount`): σ = (y(0,·) + y(1,·))/2, δ = (y(0,·) − y(1,·))/2 on G'' with their own compressions; S-join, then a
  D-join per S-solution, every solution re-checked.
- They agree on all 27 multisets of 347's slice (0 to 4,164 solutions each).
- On the H' = 1 controls (v = 35, q = 5; v = 47, q = 7, s = 2 and 3), where z = y, the exact totals reproduce the
  earlier canonical counts (`logs/ab/paley_c6.jsonl`): 622,080 ordered = **1,728** canonical, 6,912 = **36**,
  622,080 = **1,728**.

### 2.2 Totals

`log2 sym` is the order of the symmetry group (§1); "per class" = log2(total) − log2 sym.

| rung | v | L1 | rest lags | multisets | log2 sym | estimate | exact | per class: est / exact | error (bits) |
|---|---|---|---|---|---|---|---|---|---|
| 347's slice | 347 | 10 | 6 | 27 | 13.6 | 2^20.97 | 2^22.85 | +7.4 / +9.3 | −1.9 |
| c = 6, −1 on Q's Z3, H' = Z3, q = 11 | 203 | 12 | 4 | 4 | 22.9 | 2^20.54 | 2^24.15 | −2.4 / +1.2 | −3.6 |
| c = 8, H' = Z3 (−1), q = 5 | 127 | 16 | 8 | 8 | 28.0 | 2^35.56 | 2^32.90 | +7.6 / +4.9 | +2.7 |
| c = 8, H' = Z3 (−1), q = 7 | 175 | 16 | 8 | 24 | 28.0 | 2^1.23 | **0** | −26.8 / – | – |
| c = 6, −1 on H' = Z3, q = 7 | 131 | 12 | 4 | 27 | 22.9 | 2^18.93 | **0** | −4.0 / – | – |
| c = 6, −1 on H' = Z5, q = 11 | 335 | 18 | 8 | 62 | 23.9 | 2^9.59 | **0** | −14.3 / – | – |
| c = 6, −1 on H' = Z5, q = 13 | 395 | 18 | 8 | 17 | 23.9 | 2^14.38 | **0** | −9.6 / – | – |
| c = 6, −1 on Q's Z3 and H' = Z5, q = 13 | 395 | 16 | 12 | 8 | 14.6 | 2^−21.5 | **0** | −36.1 / – | – |
| c = 4, H' = Z4, q = 13 | 211 | 16 | 8 | 1 | 31.2 | 0 (outside the lattice) | **0** | – | – |

Seven more rungs have no compression quadruple at all (stage 1 empty [M]: c = 6 at v = 203, 239 with −1 on H' = Z3;
v = 239, 311 with −1 on Q's Z3; c = 8 at v = 271, 319, 463). Eleven rungs were too rich for the exact counters
(estimates 2^13–2^75 ordered, 2^−13.5 to 2^+45 per class), and the twelfth, v = 335 with −1 on Q's Z3 and H' = Z5,
was counted exactly in 11 of its 25 multisets, all empty (estimate 2^−26 per class); they enter only the
K-subsampling test below.

- **Where solutions exist the estimate is within ±3.6 bits, in both directions.** The sign is not stable (347 and
  v = 203 low, v = 127 high), so no correction is applied.
- **The ordered totals of the empty rungs look large (2^9.6–2^18.9) but are ≤ 2^−4 per symmetry class.** Read per
  class, the estimator never called an empty rung rich.
- These agree with `report/b599_answer.md` §3, whose numbers are per class (it summed one estimate per symmetry class
  of quadruples): its v = 335 rung (2^−2.3 per its classes) is exactly empty here too, and its v = 215 rich quadruples
  were within 0.6 bit.

### 2.3 Per multiset, number of multisets, symmetric solutions

- **Per multiset the estimate is not reliable.** On 347's slice, 12 of the 27 multisets are exactly empty although
  each was estimated at 2^2–2^8 solutions; where both are positive the error runs from −3.7 to +3.4 bits (median
  −1.7). At v = 203, 2 of 4 multisets are empty; the others are 2.9–8.2 bits low. The totals average this out.
- **Degenerate multisets** (blocks with one or a few completions, as at v = 131: block lists of 222, 1, 1, 1) are where
  the Gaussian is meaningless; the production screen counts those exactly instead (any block list below 30,000).
- **Number of multisets.** Drawing K multisets with probability ∝ (ordered quadruples × exact lag-0 count) and scaling
  the mean Gaussian probability by the exact lag-0 total: the 5–95% range of log2(K-estimate / all-multiset estimate)
  over 2,000 draws is within ±0.4 bit at K = 40 on 15 of the 20 rungs with multisets and within ±0.6 on 17. The three
  exceptions are exactly empty rungs whose weight sits in a few multisets: v = 335 (−1 on H' = Z5; −1.8 … +1.3),
  v = 175 (c = 8; −3.1 … +1.4) and v = 395 (−1 on H' = Z5; −7.5 … +1.3). The subsample misses such multisets more
  often than it overweights them, so the error is mostly downward; at K = 20 the upper end stays below +1.9 bits. The
  production screen uses all multisets when there are at most 60, else K = 40 (the first 49 configurations) or K = 24
  (`logs/s1z/calib/calib_summary.json`, `K_bootstrap_vs_all`). The one importance-sampled verdict near the margin,
  the 599 lead, was rerun over all 85 multisets (§4.4).
- **Symmetric solutions** (every block fixed by an affine map g ↦ φ(g) + t of G' that moves some U-orbit): 4.9% of
  347's 7,550,976 ordered solutions; 0% at v = 203 and v = 127. The scarce-regime excess of
  `report/b599_answer.md` §3 (93% symmetric, Gaussian 1.3–5.5 bits low) did not recur in these rungs.

### 2.4 Margin

A configuration is called **rich [E]** when the z-conditioned estimate is at least 2^+4 solutions per symmetry class,
**dead [E]** when at most 2^−4 (every exactly empty rung lies there; the largest observed error where solutions exist is
3.6 bits), and undecided in between. Exact zero (no compression quadruple, or no energy split) is [M].

---

## 3. The 112 configurations

`theory/s1z_table.py` (`logs/s1z/table.json`): one row per distinct configuration, sorted by v, c, q, s and L1.

**Verdicts of the 112** (first that applies):

| verdict | configurations |
|---|---|
| stage 1 empty [M] (no compression quadruple; all three also [P]: 739 by the parity argument of §4.1, the two 919 s = 3 slices by the energy, §8.1) | 3 |
| no core [M] (exact stage-2 compression screens, §5) | 1 |
| stage 1 at most 2^−4 per symmetry class [E] | 19 |
| stage 1 empty [E] (target outside the lattice in every multiset, §4.3) | 1 |
| stage 2 over 2^47 per core [E] | 88 |
| open | 0 |

Among the 88 over budget, 85 are rich (at least 2^+4 per symmetry class) and 3 are in the
undecided band (the 719 q = 7 pair with −1 on Q's Z3 and on Z17 at 2^+3.9 per class, 2^70 per core, and 839 c = 8,
s = 2, ×5 on Z8 at 2^+3.1, 2^142); the band's fourth member, the 719 lead, has no core [M]. A symmetric-block search
for these small totals would decide stage 1 only, and stage 2 is out of reach for all three, so none was run.

**Columns.**
- v, c; q, s (f): the order, c = b + 1, the staging prime, the index s of M in F_q^* and f = (q − 1)/s.
- H' and U: H' as a product of cyclic groups; U by its generators ("Q's Z3" is the Z3 of Q = Z2 × Z3; an axis of H'
  is named by its order, with its number when H' has equal factors).
- L1: number of G'-orbits.
- stage 1: log2 of the estimated number of ordered stage-1 solutions, then per symmetry class [E] (§1, §2.4): over all
  compression multisets when there are at most 60 (84 configurations), else an importance sample of 24 or 40 multisets
  (25 configurations; in-sample relative standard error 0.8%–3.3%). "empty [M]": no compression quadruple
  meets the Q-level and lag-0 equations; "−∞ [E]": the target is outside the lattice in every multiset.
- stage 2 per y, D_def, cost per core [E]: the typical-y budget of §4.
- cores per class [E]: stage 1 per class + stage 2 per y − log2 s.
- exact stage-2 screens [M]: compression multisets with no lift at some level of §5.3 (aggregated screens, or the per-type
  screen for the 719 lead), else the per-character Z_q screens.

**Duplicates.** Five pairs are one slice each: two generator lists of the same U (719, q = 7, s = 2 and 3: −1 on Q's
Z3 with −1 on Z17, and −1 on Q's Z3 × Z17 with −1 on Z17), and a swap of the two equal axes of H' (599 with
H' = Z3², s = 2 and 5; 643 with H' = Z4 × Z4 × Z2). The 719 pairs give identical rows (same group, same seeds); the
axis-swap pairs differ by estimator noise only: stage 1 by at most 0.6 bit, the typical-y cost (one random
typical y each) by 3.9–9.9 bits, more for the larger budgets (§6). So the 112 configurations are 107
slices.

| v, c | q, s (f) | H' | U | L1 | stage 1: log2 ordered; per sym. class | stage 2 per y, D_def, cost/core [E] | cores per class [E] | exact stage-2 screens | verdict |
|---|---|---|---|---|---|---|---|---|---|
| 599, 6 | 11, 2 (5) | 9 | −1 on Q's Z3, Z9 | 28 | -24.9; -40.1 | 4.7, 17.2, 2^34.3 | -36 | per character: 0/37 dead | stage 1 < 2^-4 per symmetry class [E] |
| 599, 6 | 11, 2 (5) | 3×3 | −1 on Q's Z3, Z3(axis 1), Z3(axis 2) | 28 | -24.9; -39.5 | 4.7, 17.2, 2^34.3 | -36 | per character: 0/37 dead | stage 1 < 2^-4 per symmetry class [E] |
| 599, 6 | 11, 2 (5) | 9 | −1 on Z9 | 30 | 18.5; -6.0 | -2.2, 21.5, 2^43.5 | -9 | per character: 0/85 dead | stage 1 < 2^-4 per symmetry class [E] |
| 599, 6 | 11, 2 (5) | 9 | ×4 on Z9 | 30 | 61.1; +30.8 | -12.3, 25.8, 2^51.5 | +18 | 0/85 multisets dead; 85 pass Z11, Z2×Z11 (Z3×Z11 too large) | stage 2 over budget [E] (2^52 per core) |
| 599, 6 | 11, 2 (5) | 9 | −1 on Q's Z3 | 36 | 85.7; +54.8 | -13.8, 31.9, 2^62.0 | +40 | – | stage 2 over budget [E] (2^62 per core) |
| 599, 6 | 11, 2 (5) | 3×3 | −1 on Z3(axis 1) | 36 | 59.8; +29.5 | -11.9, 30.8, 2^60.4 | +17 | – | stage 2 over budget [E] (2^60 per core) |
| 599, 6 | 11, 2 (5) | 3×3 | −1 on Z3(axis 2) | 36 | 59.8; +29.5 | -6.5, 30.2, 2^56.5 | +22 | – | stage 2 over budget [E] (2^56 per core) |
| 599, 6 | 11, 2 (5) | 3×3 | −1 on Q's Z3 | 36 | 85.6; +55.4 | -13.1, 31.2, 2^61.3 | +41 | – | stage 2 over budget [E] (2^61 per core) |
| 599, 6 | 11, 2 (5) | 9 | 1 | 54 | 116.1; +77.9 | -9.4, 44.8, 2^84.4 | +68 | – | stage 2 over budget [E] (2^84 per core) |
| 599, 6 | 11, 2 (5) | 3×3 | 1 | 54 | 116.1; +78.5 | -9.0, 44.5, 2^84.0 | +68 | – | stage 2 over budget [E] (2^84 per core) |
| 599, 6 | 11, 5 (2) | 9 | −1 on Q's Z3, Z9 | 28 | 103.9; +88.8 | -279.4, 181.5, 2^419.8 | -193 | – | stage 2 over budget [E] (2^420 per core) |
| 599, 6 | 11, 5 (2) | 3×3 | −1 on Q's Z3, Z3(axis 1), Z3(axis 2) | 28 | 104.0; +89.4 | -279.2, 181.3, 2^419.5 | -192 | – | stage 2 over budget [E] (2^420 per core) |
| 599, 6 | 11, 5 (2) | 9 | −1 on Z9 | 30 | 151.5; +127.0 | -83.3, 123.7, 2^234.4 | +41 | – | stage 2 over budget [E] (2^234 per core) |
| 599, 6 | 11, 5 (2) | 9 | ×4 on Z9 | 30 | 175.2; +144.9 | -10.9, 97.6, 2^159.2 | +132 | – | stage 2 over budget [E] (2^159 per core) |
| 599, 6 | 11, 5 (2) | 9 | −1 on Q's Z3 | 36 | 244.9; +214.1 | -22.2, 119.7, 2^200.6 | +190 | – | stage 2 over budget [E] (2^201 per core) |
| 599, 6 | 11, 5 (2) | 3×3 | −1 on Z3(axis 1) | 36 | 217.8; +187.5 | -24.0, 121.3, 2^202.1 | +161 | – | stage 2 over budget [E] (2^202 per core) |
| 599, 6 | 11, 5 (2) | 3×3 | −1 on Z3(axis 2) | 36 | 217.8; +187.6 | -12.8, 118.1, 2^195.0 | +172 | – | stage 2 over budget [E] (2^195 per core) |
| 599, 6 | 11, 5 (2) | 3×3 | −1 on Q's Z3 | 36 | 244.9; +214.6 | -22.3, 119.6, 2^200.7 | +190 | – | stage 2 over budget [E] (2^201 per core) |
| 599, 6 | 11, 5 (2) | 9 | 1 | 54 | 365.4; +327.2 | 31.3, 156.2, 2^240.8 | +356 | – | stage 2 over budget [E] (2^241 per core) |
| 599, 6 | 11, 5 (2) | 3×3 | 1 | 54 | 365.4; +327.8 | 31.4, 156.2, 2^240.7 | +357 | – | stage 2 over budget [E] (2^241 per core) |
| 643, 4 | 5, 2 (2) | 32 | ×3 on Z32 | 36 | 88.0; +60.8 | -134.1, 79.1, 2^166.9 | -74 | – | stage 2 over budget [E] (2^167 per core) |
| 643, 4 | 5, 2 (2) | 32 | ×5 on Z32 | 40 | 24.7; -6.4 | -128.8, 82.1, 2^164.6 | -136 | – | stage 1 < 2^-4 per symmetry class [E] |
| 643, 4 | 5, 2 (2) | 32 | ×7 on Z32 | 52 | 200.5; +172.4 | -175.9, 109.0, 2^229.7 | -4 | – | stage 2 over budget [E] (2^230 per core) |
| 643, 4 | 5, 2 (2) | 16×2 | ×3 on Z16 | 56 | 202.6; +171.5 | -207.8, 126.4, 2^263.5 | -37 | – | stage 2 over budget [E] (2^264 per core) |
| 643, 4 | 5, 2 (2) | 32 | ×9 on Z32 | 64 | 263.1; +226.9 | -159.3, 108.6, 2^221.1 | +67 | – | stage 2 over budget [E] (2^221 per core) |
| 643, 4 | 5, 2 (2) | 16×2 | ×5 on Z16 | 64 | 236.0; +200.8 | -203.4, 132.5, 2^265.1 | -4 | – | stage 2 over budget [E] (2^265 per core) |
| 643, 4 | 5, 2 (2) | 32 | ×15 on Z32 | 68 | 263.2; +234.0 | -234.3, 146.7, 2^306.6 | -1 | – | stage 2 over budget [E] (2^307 per core) |
| 643, 4 | 5, 2 (2) | 32 | −1 on Z32 | 68 | 181.6; +155.5 | -330.4, 179.5, 2^403.4 | -176 | – | stage 2 over budget [E] (2^403 per core) |
| 643, 4 | 5, 2 (2) | 16×2 | ×7 on Z16 | 72 | 275.7; +243.5 | -261.1, 162.1, 2^336.6 | -19 | – | stage 2 over budget [E] (2^337 per core) |
| 643, 4 | 5, 2 (2) | 16×2 | −1 on Z16 | 72 | 195.5; +166.4 | -357.1, 202.3, 2^433.6 | -192 | – | stage 2 over budget [E] (2^434 per core) |
| 643, 4 | 5, 2 (2) | 8×4 | ×3 on Z8 | 80 | 361.0; +324.9 | -223.1, 150.4, 2^306.9 | +101 | – | stage 2 over budget [E] (2^307 per core) |
| 643, 4 | 5, 2 (2) | 8×4 | −1 on Z8 | 80 | 322.4; +286.3 | -269.4, 171.1, 2^353.9 | +16 | – | stage 2 over budget [E] (2^354 per core) |
| 643, 4 | 5, 2 (2) | 8×2×2 | ×3 on Z8 | 80 | 294.3; +259.2 | -307.2, 186.8, 2^391.0 | -49 | – | stage 2 over budget [E] (2^391 per core) |
| 643, 4 | 5, 2 (2) | 8×2×2 | −1 on Z8 | 80 | 217.0; +184.9 | -395.0, 224.7, 2^479.5 | -211 | – | stage 2 over budget [E] (2^480 per core) |
| 643, 4 | 5, 2 (2) | 32 | ×17 on Z32 | 96 | 489.6; +448.4 | -195.8, 151.8, 2^297.5 | +252 | – | stage 2 over budget [E] (2^298 per core) |
| 643, 4 | 5, 2 (2) | 16×2 | ×9 on Z16 | 96 | 466.7; +426.5 | -224.4, 166.2, 2^327.6 | +201 | – | stage 2 over budget [E] (2^328 per core) |
| 643, 4 | 5, 2 (2) | 8×4 | ×5 on Z8 | 96 | 466.5; +426.3 | -220.2, 162.5, 2^323.4 | +205 | – | stage 2 over budget [E] (2^323 per core) |
| 643, 4 | 5, 2 (2) | 8×4 | −1 on Z4 | 96 | 448.4; +408.2 | -241.2, 169.1, 2^341.7 | +166 | – | stage 2 over budget [E] (2^342 per core) |
| 643, 4 | 5, 2 (2) | 8×2×2 | ×5 on Z8 | 96 | 421.4; +382.2 | -277.1, 186.3, 2^380.3 | +104 | – | stage 2 over budget [E] (2^380 per core) |
| 643, 4 | 5, 2 (2) | 4×4×2 | −1 on Z4(axis 1) | 96 | 385.6; +346.5 | -307.5, 200.2, 2^410.2 | +38 | – | stage 2 over budget [E] (2^410 per core) |
| 643, 4 | 5, 2 (2) | 4×4×2 | −1 on Z4(axis 2) | 96 | 385.1; +345.9 | -315.3, 206.5, 2^420.1 | +30 | – | stage 2 over budget [E] (2^420 per core) |
| 643, 4 | 5, 2 (2) | 4×2×2×2 | −1 on Z4 | 96 | 259.3; +224.1 | -456.8, 262.4, 2^559.5 | -234 | – | stage 2 over budget [E] (2^560 per core) |
| 643, 4 | 5, 2 (2) | 32 | 1 | 128 | 689.3; +643.1 | -227.6, 189.2, 2^367.3 | +415 | – | stage 2 over budget [E] (2^367 per core) |
| 643, 4 | 5, 2 (2) | 16×2 | 1 | 128 | 668.8; +623.6 | -252.4, 201.0, 2^392.2 | +370 | – | stage 2 over budget [E] (2^392 per core) |
| 643, 4 | 5, 2 (2) | 8×4 | 1 | 128 | 668.8; +623.6 | -253.1, 200.9, 2^392.9 | +370 | – | stage 2 over budget [E] (2^393 per core) |
| 643, 4 | 5, 2 (2) | 8×2×2 | 1 | 128 | 627.6; +583.4 | -302.4, 220.8, 2^442.2 | +280 | – | stage 2 over budget [E] (2^442 per core) |
| 643, 4 | 5, 2 (2) | 4×4×2 | 1 | 128 | 627.6; +583.4 | -305.7, 224.0, 2^445.5 | +277 | – | stage 2 over budget [E] (2^445 per core) |
| 643, 4 | 5, 2 (2) | 4×2×2×2 | 1 | 128 | 545.1; +502.0 | -400.9, 260.0, 2^540.6 | +100 | – | stage 2 over budget [E] (2^541 per core) |
| 643, 4 | 5, 2 (2) | 2×2×2×2×2 | 1 | 128 | 380.3; +341.1 | -595.1, 338.6, 2^734.9 | -255 | – | stage 2 over budget [E] (2^735 per core) |
| 719, 6 | 7, 2 (3) | 17 | ×4 on Z17 | 30 | 50.3; +25.4 | -37.0, 35.1, 2^71.5 | -13 | – | stage 2 over budget [E] (2^72 per core) |
| 719, 6 | 7, 2 (3) | 17 | −1 on Q's Z3; −1 on Z17 | 36 | 19.4; +3.9 | -26.9, 36.5, 2^70.3 | -24 | – | stage 2 over budget [E] (2^70 per core) |
| 719, 6 | 7, 2 (3) | 17 | −1 on Q's Z3, Z17; −1 on Z17 | 36 | 19.4; +3.9 | -26.9, 36.5, 2^70.3 | -24 | – | stage 2 over budget [E] (2^70 per core) |
| 719, 6 | 7, 2 (3) | 17 | −1 on Q's Z3, Z17 | 52 | 37.7; +21.2 | -29.1, 48.2, 2^94.9 | -9 | – | stage 2 over budget [E] (2^95 per core) |
| 719, 6 | 7, 2 (3) | 17 | −1 on Z17 | 54 | 127.4; +101.5 | -50.9, 57.5, 2^116.4 | +50 | – | stage 2 over budget [E] (2^116 per core) |
| 719, 6 | 7, 2 (3) | 17 | −1 on Q's Z3 | 68 | 264.9; +228.9 | -78.8, 83.5, 2^162.8 | +149 | – | stage 2 over budget [E] (2^163 per core) |
| 719, 6 | 7, 2 (3) | 17 | 1 | 102 | 388.3; +345.0 | -100.9, 118.7, 2^231.2 | +243 | – | stage 2 over budget [E] (2^231 per core) |
| 719, 6 | 7, 3 (2) | 17 | ×4 on Z17 | 30 | 109.9; +85.0 | -108.8, 88.6, 2^177.3 | -25 | – | stage 2 over budget [E] (2^177 per core) |
| 719, 6 | 7, 3 (2) | 17 | −1 on Q's Z3; −1 on Z17 | 36 | 94.0; +78.4 | -276.5, 157.2, 2^359.8 | -200 | – | stage 2 over budget [E] (2^360 per core) |
| 719, 6 | 7, 3 (2) | 17 | −1 on Q's Z3, Z17; −1 on Z17 | 36 | 94.0; +78.4 | -276.5, 157.2, 2^359.8 | -200 | – | stage 2 over budget [E] (2^360 per core) |
| 719, 6 | 7, 3 (2) | 17 | −1 on Q's Z3, Z17 | 52 | 151.8; +135.2 | -386.7, 222.3, 2^507.2 | -253 | – | stage 2 over budget [E] (2^507 per core) |
| 719, 6 | 7, 3 (2) | 17 | −1 on Z17 | 54 | 245.6; +219.6 | -176.7, 153.6, 2^302.6 | +41 | – | stage 2 over budget [E] (2^303 per core) |
| 719, 6 | 7, 3 (2) | 17 | −1 on Q's Z3 | 68 | 407.5; +371.6 | -121.4, 156.5, 2^276.7 | +249 | – | stage 2 over budget [E] (2^277 per core) |
| 719, 6 | 7, 3 (2) | 17 | 1 | 102 | 613.7; +570.4 | -134.2, 212.8, 2^370.6 | +435 | – | stage 2 over budget [E] (2^371 per core) |
| 719, 6 | 17, 2 (8) | 7 | −1 on Q's Z3 | 28 | 27.3; -2.1 | -8.6, 20.6, 2^48.1 | -12 | 39/39 multisets dead (all 242,977 type quadruples; §5) | no core [M] (exact stage-2 compression screens) |
| 719, 6 | 17, 2 (8) | 7 | 1 | 42 | 26.9; -9.8 | -0.5, 28.6, 2^60.7 | -11 | 85/89 multisets dead; 4 pass Z17, Z2×Z17, Z3×Z17 (Q×Z17 too large) | stage 1 < 2^-4 per symmetry class [E] |
| 719, 6 | 17, 4 (4) | 7 | −1 on Q's Z3 | 28 | 103.2; +73.8 | -16.6, 71.0, 2^126.0 | +55 | – | stage 2 over budget [E] (2^126 per core) |
| 719, 6 | 17, 4 (4) | 7 | 1 | 42 | 149.0; +112.3 | 16.1, 92.1, 2^148.9 | +126 | – | stage 2 over budget [E] (2^149 per core) |
| 719, 6 | 17, 8 (2) | 7 | 1 | 42 | 313.0; +276.2 | 135.9, 174.3, 2^253.9 | +409 | – | stage 2 over budget [E] (2^254 per core) |
| 739, 4 | 23, 2 (11) | 8 | 1 | 32 | empty [M] | 17.3, 15.2, 2^29.5 | – | – | stage 1 empty [M] |
| 739, 4 | 23, 11 (2) | 8 | 1 | 32 | 256.0; +219.8 | 66.7, 223.7, 2^367.6 | +283 | – | stage 2 over budget [E] (2^368 per core) |
| 823, 8 | 17, 2 (8) | 3×2 | −1 on Z3 | 32 | -57.2; -89.2 | -85.7, 61.0, 2^131.4 | -176 | – | stage 1 < 2^-4 per symmetry class [E] |
| 823, 8 | 17, 2 (8) | 3×2 | 1 | 48 | 3.0; -39.3 | -38.6, 48.9, 2^108.1 | -79 | – | stage 1 < 2^-4 per symmetry class [E] |
| 823, 8 | 17, 4 (4) | 3×2 | −1 on Z3 | 32 | 7.7; -24.3 | -229.0, 170.3, 2^353.5 | -255 | – | stage 1 < 2^-4 per symmetry class [E] |
| 823, 8 | 17, 4 (4) | 3×2 | 1 | 48 | 134.0; +91.7 | -91.1, 153.2, 2^279.7 | -1 | – | stage 2 over budget [E] (2^280 per core) |
| 823, 8 | 17, 8 (2) | 3×2 | −1 on Z3 | 32 | 129.2; +97.2 | -410.6, 350.4, 2^706.0 | -316 | – | stage 2 over budget [E] (2^706 per core) |
| 823, 8 | 17, 8 (2) | 3×2 | 1 | 48 | 321.4; +279.1 | -78.0, 303.8, 2^522.9 | +198 | – | stage 2 over budget [E] (2^523 per core) |
| 839, 8 | 13, 2 (6) | 8 | ×3 on Z8 | 40 | -6.6; -42.6 | -99.0, 69.0, 2^153.6 | -143 | – | stage 1 < 2^-4 per symmetry class [E] |
| 839, 8 | 13, 2 (6) | 8 | −1 on Z8 | 40 | -53.9; -86.8 | -132.5, 79.7, 2^186.8 | -220 | – | stage 1 < 2^-4 per symmetry class [E] |
| 839, 8 | 13, 2 (6) | 8 | ×5 on Z8 | 48 | 43.0; +3.1 | -75.5, 66.6, 2^142.0 | -74 | – | stage 2 over budget [E] (2^142 per core) |
| 839, 8 | 13, 2 (6) | 4×2 | −1 on Z4 | 48 | -51.6; -87.6 | -157.6, 96.8, 2^223.9 | -246 | – | stage 1 < 2^-4 per symmetry class [E] |
| 839, 8 | 13, 2 (6) | 8 | 1 | 64 | 72.1; +27.1 | -56.3, 70.1, 2^147.3 | -30 | – | stage 2 over budget [E] (2^147 per core) |
| 839, 8 | 13, 2 (6) | 4×2 | 1 | 64 | 23.6; -20.3 | -98.3, 85.3, 2^189.3 | -120 | – | stage 1 < 2^-4 per symmetry class [E] |
| 839, 8 | 13, 2 (6) | 2×2×2 | 1 | 64 | -73.1; -113.1 | -188.2, 121.9, 2^279.2 | -302 | – | stage 1 < 2^-4 per symmetry class [E] |
| 839, 8 | 13, 3 (4) | 8 | ×3 on Z8 | 40 | 6.5; -29.5 | -163.0, 132.1, 2^265.9 | -194 | – | stage 1 < 2^-4 per symmetry class [E] |
| 839, 8 | 13, 3 (4) | 8 | −1 on Z8 | 40 | -39.0; -72.0 | -235.6, 161.3, 2^338.6 | -309 | – | stage 1 < 2^-4 per symmetry class [E] |
| 839, 8 | 13, 3 (4) | 8 | ×5 on Z8 | 48 | 31.7; -8.3 | -121.0, 127.3, 2^244.1 | -131 | – | stage 1 < 2^-4 per symmetry class [E] |
| 839, 8 | 13, 3 (4) | 4×2 | −1 on Z4 | 48 | -64.7; -100.6 | -278.9, 194.7, 2^402.5 | -381 | – | stage 1 < 2^-4 per symmetry class [E] |
| 839, 8 | 13, 3 (4) | 8 | 1 | 64 | 121.2; +76.2 | -111.2, 147.4, 2^275.9 | -37 | – | stage 2 over budget [E] (2^276 per core) |
| 839, 8 | 13, 3 (4) | 4×2 | 1 | 64 | 72.9; +28.9 | -194.9, 182.7, 2^359.5 | -168 | – | stage 2 over budget [E] (2^360 per core) |
| 839, 8 | 13, 3 (4) | 2×2×2 | 1 | 64 | -23.8; -63.8 | -363.4, 255.0, 2^528.0 | -429 | – | stage 1 < 2^-4 per symmetry class [E] |
| 839, 8 | 13, 4 (3) | 8 | ×3 on Z8 | 40 | 106.7; +70.8 | 20.3, 85.9, 2^132.6 | +89 | – | stage 2 over budget [E] (2^133 per core) |
| 839, 8 | 13, 4 (3) | 8 | −1 on Z8 | 40 | 60.0; +27.0 | 38.1, 77.4, 2^116.1 | +63 | – | stage 2 over budget [E] (2^116 per core) |
| 839, 8 | 13, 4 (3) | 8 | ×5 on Z8 | 48 | 179.4; +139.4 | 16.2, 106.1, 2^168.0 | +154 | – | stage 2 over budget [E] (2^168 per core) |
| 839, 8 | 13, 4 (3) | 4×2 | −1 on Z4 | 48 | 82.9; +46.9 | 45.3, 94.0, 2^139.2 | +90 | – | stage 2 over budget [E] (2^139 per core) |
| 839, 8 | 13, 4 (3) | 8 | 1 | 64 | 271.7; +226.7 | 26.3, 142.5, 2^221.0 | +251 | – | stage 2 over budget [E] (2^221 per core) |
| 839, 8 | 13, 4 (3) | 4×2 | 1 | 64 | 223.4; +179.4 | 40.4, 135.3, 2^206.8 | +218 | – | stage 2 over budget [E] (2^207 per core) |
| 839, 8 | 13, 4 (3) | 2×2×2 | 1 | 64 | 126.9; +86.9 | 68.7, 121.2, 2^178.6 | +154 | – | stage 2 over budget [E] (2^179 per core) |
| 839, 8 | 13, 6 (2) | 8 | ×3 on Z8 | 40 | 197.2; +161.2 | -271.4, 277.1, 2^528.3 | -113 | – | stage 2 over budget [E] (2^528 per core) |
| 839, 8 | 13, 6 (2) | 8 | −1 on Z8 | 40 | 150.6; +117.6 | -432.7, 344.1, 2^688.7 | -318 | – | stage 2 over budget [E] (2^689 per core) |
| 839, 8 | 13, 6 (2) | 8 | ×5 on Z8 | 48 | 285.6; +245.7 | -164.1, 265.3, 2^471.7 | +79 | – | stage 2 over budget [E] (2^472 per core) |
| 839, 8 | 13, 6 (2) | 4×2 | −1 on Z4 | 48 | 189.4; +153.4 | -498.8, 410.2, 2^809.1 | -348 | – | stage 2 over budget [E] (2^809 per core) |
| 839, 8 | 13, 6 (2) | 8 | 1 | 64 | 419.3; +374.4 | -97.5, 295.9, 2^509.1 | +274 | – | stage 2 over budget [E] (2^509 per core) |
| 839, 8 | 13, 6 (2) | 4×2 | 1 | 64 | 371.1; +327.1 | -281.0, 375.2, 2^692.6 | +44 | – | stage 2 over budget [E] (2^693 per core) |
| 839, 8 | 13, 6 (2) | 2×2×2 | 1 | 64 | 274.6; +234.6 | -647.0, 532.0, 2^1058.7 | -415 | – | stage 2 over budget [E] (2^1059 per core) |
| 919, 8 | 19, 2 (9) | 3×2 | −1 on Z3 | 32 | −∞ [E] | 6.8, 19.2, 2^38.7 | – | per character: 0/15 dead | stage 1 empty [E] (target outside the lattice of every multiset) |
| 919, 8 | 19, 2 (9) | 3×2 | 1 | 48 | -13.4; -55.8 | 9.6, 30.3, 2^59.9 | -47 | – | stage 1 < 2^-4 per symmetry class [E] |
| 919, 8 | 19, 3 (6) | 3×2 | −1 on Z3 | 32 | empty [M] | – | – | – | stage 1 empty [M] |
| 919, 8 | 19, 3 (6) | 3×2 | 1 | 48 | empty [M] | – | – | – | stage 1 empty [M] |
| 919, 8 | 19, 6 (3) | 3×2 | −1 on Z3 | 32 | 60.5; +28.5 | 86.3, 89.7, 2^124.1 | +112 | – | stage 2 over budget [E] (2^124 per core) |
| 919, 8 | 19, 6 (3) | 3×2 | 1 | 48 | 219.3; +176.9 | 110.0, 144.8, 2^207.0 | +284 | – | stage 2 over budget [E] (2^207 per core) |
| 919, 8 | 19, 9 (2) | 3×2 | −1 on Z3 | 32 | 132.5; +100.6 | -465.6, 405.3, 2^806.9 | -368 | – | stage 2 over budget [E] (2^807 per core) |
| 919, 8 | 19, 9 (2) | 3×2 | 1 | 48 | 329.8; +287.5 | -75.9, 341.9, 2^590.3 | +208 | – | stage 2 over budget [E] (2^590 per core) |

---

## 4. The six configurations with a stage-2 budget of at most 2^48.1

Costs by `theory/s1z_budget_batch.py` (`logs/s1z/budget_skipped.jsonl`): the typical-y level budget of
`theory/ab_typical_levels.py` with no stage-2 orbit cap, and the model 2^{D_def + C_S − H_S/2} calibrated at 347 and 523
(it reproduces `report/b599_answer.md`'s 2^43.5 for the 599 lead). Stage 1 by `theory/s1z_screen.py`
(`logs/s1z/zscreen_skipped.jsonl`).

| configuration | G'-orbits | stage 2 per y, D_def, cost per core [E] | stage 1, ordered; per symmetry class | status |
|---|---|---|---|---|
| 739, c = 4, q = 23, s = 2, H' = Z8, U = 1 | 32 | +17.3, 15.2, 2^29.5 | no compression quadruple | **stage 1 empty [M, P]** (§4.1) |
| 599, c = 6, q = 11, s = 2, H' = Z3², −1 on Q's Z3 and H' | 28 | +4.7, 17.2, 2^34.3 | 2^−24.9; 2^−39.5 | dead [E] |
| 599, c = 6, q = 11, s = 2, H' = Z9, −1 on Q's Z3 and H' | 28 | +4.7, 17.2, 2^34.3 | 2^−24.9; 2^−40.1 | dead [E] |
| 919, c = 8, q = 19, s = 2, H' = Z3 × Z2, −1 on Z3 (built order) | 32 | +6.8, 19.2, 2^38.7 | target outside the lattice in all 15 multisets | dead [E] (mod-3 congruence, §4.3) |
| 599, c = 6, q = 11, s = 2, H' = Z9, −1 on Z9 (the lead) | 30 | −2.2, 21.5, 2^43.5 | 2^18.5; 2^−6.0 | dead [E]; symmetric families empty [M] (b599 §4) |
| **719, c = 6, q = 17, s = 2, H' = Z7, −1 on Q's Z3** | 28 | −8.6, 20.6, 2^48.1 | 2^27.3; 2^−2.1 | **no core [M]** (§5) |

### 4.1 739, c = 4, q = 23, s = 2: stage 1 is empty [P]

Letters ±1 + 11·{−2, 0, 2} = {±1, ±21, ±23}; G' = Z2² × Z8, all 32 G'-orbits are points; T1(0) = 2680; 4v = 2956.
- Energy: Σ_k Σ_g y_k(g)² = 2680 over 128 entries, so 440·n21 + 528·n23 = 2552, i.e. 5·n21 + 6·n23 = 29: exactly
  **one** entry of size 21 and **four** of size 23 — five large entries in all.
- Cosets: each block has four Q-cosets of eight entries. If the signed sum of a coset's large entries is in
  [−4, 4], the coset sum has |sum| ≤ 4 + (number of small entries) ≤ 8, and the number of large entries is even (an
  odd number of values from ±21, ±23 sums to an odd number of absolute value at least 19). Otherwise that signed sum
  has absolute value at least 19 (odd count) or 40 (even count), and |sum| ≥ 14.
- Q-level: Σ_k |ẑ_k(ψ)|² = 12 for ψ ≠ 1, so each block has Σ_q (z(q) − z̄)² ≤ 36/4 = 9 (Parseval), hence all four
  coset sums lie within 6 of each other. A block with both kinds of coset needs sums 14 and 8 (same sign) and then has
  Σ (z − z̄)² ≥ 27 > 9; a block with four cosets of the second kind has |row sum| ≥ 56 > √2956 = 54.4. So every coset
  of every block is of the first kind, with an even number of large entries, and the total is even. Five is odd. ∎
The argument uses only the coset size, the letters, T1(0) and 4v, so the other H' of order 8 with U trivial (not in the
skipped list) are empty too. `theory/ab_sweep.py` screen A and `theory/s1z_screen.py` (no compression quadruple)
agree [M].

### 4.2 The two 599 slices with −1 on Q's Z3 and on H'

Both have 37 compression multisets (296 compressions); the Gaussian of the 24 remaining lag classes leaves
2^−24.9 ordered solutions, 2^−40 per symmetry class [E] — 36 bits under the margin of §2.4. (U acts on Q, so
`theory/b599_symgen.py` does not apply; the per-character Z11 screens leave all 37 multisets alive [M],
`logs/s1z/agg/char_screens.jsonl`.) The joint compression screens were not run here: 36 bits below the calibrated band, stage 1 is decided [E].

### 4.3 919, c = 8, q = 19, s = 2, −1 on Z3 (the order is built)

All 15 compression multisets have the stage-1 target outside the integer lattice spanned by 800 sampled
y-quadruples. The obstruction is a congruence mod 3 of λ·Σ_k PAF_k, λ = the sum over the 14 lags with Z3-coordinate 0
and nonzero Q-part (`logs/s1z/lattice_919_probe.log`): every sample gives 1 (mod 3), the target 2. Since all letters
are ±1 mod 3, λ·PAF(y) is a function of y mod 3 and of the Z3-compression of y; a proof that the residue is fixed by
the compressions was not completed, so this stays [E]. The per-character Z19 screens leave all 15 multisets alive [M].

### 4.4 The 599 lead

Stage 1 is 2^18.5 ordered solutions, 2^−6.0 per symmetry class [E], over all 85 multisets
(`logs/s1z/zscreen_599lead_all.jsonl`; the importance sample of 40 in `logs/s1z/zscreen_skipped.jsonl` gives the same
2^18.51, because the Gaussian of the 16 remaining lag classes costs 94.4–96.0 bits in every multiset). This is
`report/b599_answer.md`'s 2^4.9, which was a sum of one estimate per symmetry class of compression quadruples. Its
symmetric families (order-3 elements over all quadruples; order 3 and the 15-dimensional involutions over the bordered
classes) are empty [M] (b599 §4). The per-character Z11 screens leave all 85 multisets alive [M] (`logs/s1z/agg/char_screens.jsonl`).

---

## 5. The 719 lead: no core [M]

### 5.1 The slice

v = 719, b = 5, c = 6, h = 119 = 7·17, G = Z2 × Z3 × Z7 × Z17, H = Z7 × Z17, Q = Z2 × Z3. Staging prime q = 17,
M = the squares of F17^* (s = 2, f = 8), U = −1 on Q's Z3 (U acts on Q, as in 347's slice). G' = Z2 × Z3 × Z7 has
28 G'-orbits (14 points (a, 0, h) and 14 pairs {(a, ±1, h)}), 84 S-orbits per block; letters ±1 + 8·{−2, 0, 2} =
{±1, ±15, ±17}; T1 = 2536 at 0, −340 on H' ∖ 0, +68 off H' (`logs/s1z/c719/cfg_719_q17_s2_negQ3.json`). This is the
only skipped configuration whose stage 1 is not negligible and whose stage 2 costs about 2^48 per core
(H2 103, C2 111.6, −8.6 per y, D_def 20.6, C_S 67.0, h_D 24.0) [E].

### 5.2 Stage 1 [M for the counts, E for the total]

- 162 compressions, 47 classes, 50 class multisets with Σ PAF_Q = T_Q, 39 of them with a non-zero lag-0 count;
  56,064 ordered compression quadruples; exact lag-0 count 2^101.0.
- Gaussian of the 12 remaining lag classes over all 39 multisets: **2^27.3 ordered stage-1 solutions**
  (+3.4 to +17.0 bits per multiset; `logs/s1z/c719/zest_719_q17_s2_negQ3.jsonl`). The symmetry group has order
  2^29.4 (14 U-fixed translations, sign and reversal per block, block order, three multipliers of Z7 modulo ±1), so
  this is **2^−2.1 per symmetry class**: inside the undecided band of §2.4.
- Borders (`theory/s1z_border.py`, `logs/s1z/c719/border_719_q17_s2_negQ3.jsonl`): 18 of the 39 multisets have a
  member with a ±1 forced corner, 13 against 347's template frame and 5 monomial ones against the Singer templates
  (profile (4, 2, 0, 0)); they carry 2^24.0 of the 2^27.3.
- Sampling failed: the α-list sampler `theory/ab_s1_sample.py` needs 1.5M α-vectors per block and runs out of memory
  (962 MB for one array); CP-SAT with the compression quadruple of the best bordered multiset fixed
  (`theory/s1z_cpsat.py`) found nothing in 4 × 200 s (`logs/s1z/c719/cpsat_ms5_4x200s.log`). By the energy, every
  stage-1 solution has exactly ten large entries (two of size 17, eight of size 15, counted on G'), all the rest ±1.

### 5.3 Stage 2 by type-level compression screens [P, M]

**The conditions [P].** For a quotient Qbar of Q (Qbar = 1, Z2, Z3, Q), compress a core along the kernel of
Q → Qbar and along H': x̄_k(p, t) = Σ of x_k over the points of G above (p, t) ∈ Qbar × Z17. Summing the core identity
over the fibres of G → Qbar × Z17 gives Σ_k PAF_{Qbar × Z17}(x̄_k) = T', T' the fibre sums of T
(`report/paley6_answer.md` §3 for c = 4). Each is a necessary condition.

**Types [P].** A lift of a stage-1 solution chooses a fibre pattern independently for every G'-orbit j and block, and
every point of orbit j carries the same pattern (S acts trivially on the staging coordinate of G'), so orbit j adds
n_{j,p} times its pattern to x̄_k(p, ·). Hence the distribution of x̄_k over the lifts depends only on the multiset of
(orbit class, letter) pairs of block k, its **type** (orbit class = coset profile and fibre parts), and the number of
lift quadruples meeting a condition is a function of the **type quadruple**. A type quadruple with count 0 has no
lift, for every stage-1 solution that has it.

**Search [M]** (`theory/s1z_qzq.py`, `theory/s1z_typescreen.py`). For each of the 39 multisets (one compression
quadruple each, which suffices by §1), every block type with that compression (per coset component, letter multisets
per orbit class with the coset sums and energy at most T1(0) − 3·E_min), every type quadruple with Σ_k E_k = T1(0),
and for each the exact count at Qbar = 1, Z2, Z3, Q in turn (DP per type with exact integer counts, PAF signatures,
energy-split 2+2 join on (energy, 64-bit hash); a collision can only add, so zero is exact).

**Result [M]:** 242,977 type quadruples; **every one has count 0 at some level** — 5,991 at the Z17 level, 231,936
at Z2 × Z17, 0 at Z3 × Z17, 5,050 at Q × Z17 (`logs/s1z/c719/typescreen_719_summary.jsonl`). **The slice has no core.**

**Cross-check [M]** (`theory/s1z_zqscreen.py`): per compression quadruple, a DP over the orbits of each block gives
the set of (energy, x̄) over all letter assignments and all lifts at once; one energy-bucketed join per level decides
whether any type quadruple survives. It kills all 39 multisets (4 at Z17, 33 at Z2 × Z17, 2 at Q × Z17;
`logs/s1z/c719/aggscreen_full_719.jsonl`), and at every level it agrees with the per-type screen
(`logs/s1z/c719/aggscreen_validation.txt`).

**Validation [M].**
- The compression identity at all four levels holds exactly on the published 347 core (q = 19).
- Planted lifts of random S-invariant x on three 719-shaped analogs (v = 131, 95, 83), every level: the hashed
  count equals an independent count (every lift enumerated, PAF computed directly on the full vector, exact integer
  join without hashing) in 36 of 36 cases, and every planted y has count ≥ 1.
- Positive control: at v = 47 (c = 6, q = 7, s = 2, H' = 1, where Q × Z7 is all of stage 2) the count at the y of
  the known core is 5,184, the number of cores the earlier exact stage-2 run found for that y
  (`logs/ab/mono6/ctrl_v47_q7_s2.jsonl`); the type screen keeps its type quadruple.

### 5.4 What this says about the estimators

The typical-y Gaussian put stage 2 at −8.6 bits per y and the z-conditioned estimate put stage 1 at 2^27.3 ordered.
Read per symmetry class, the expected number of cores was 2^(27.3 − 8.6 − 29.4 − 1) = 2^−11.7 (the last bit: the
s = 2 multipliers of Z17 outside M, which fix y). So the estimates never predicted a core here; the exact screens
show there is none.

---

## 6. Rich but over budget [E]

"Rich" is an estimated stage 1 of at least 2^+4 solutions per symmetry class (§2.4). Costs are the model
2^{D_def + C_S − H_S/2} at a typical y (`logs/s1z/budget_skipped.jsonl`); for the two cheapest rich slices, six
typical y's each (`logs/s1z/budget_seeds_599rich.log`). Expected cores per symmetry class = stage 1 per class + stage 2
per y − log2 s (§5.4).

| configuration | G'-orbits | stage 1: ordered; per class | stage 2 per y | D_def | cost per core | cores per class | exact screens [M] |
|---|---|---|---|---|---|---|---|
| 599, c = 6, q = 11, s = 2, H' = Z9, ×4 on Z9 | 30 | 2^61.1; 2^+30.8 | −10.4 to −14.4 | 24.3–27.5 | **2^49.9–2^52.6** (six y's) | +15 to +19 | all 85 multisets pass Z11 and Z2 × Z11; Z3 × Z11 does not fit in memory |
| 599, c = 6, q = 11, s = 2, H' = Z3², −1 on one Z3 (two keys, isomorphic) | 36 | 2^59.8; 2^+29.5 | −6.5 to −12.5 | 29.7–30.9 | 2^56.5–2^61.2 (six y's) | +16 to +22 | – |
| 599, c = 6, q = 11, s = 2, H' = Z3², −1 on Q's Z3 | 36 | 2^85.6; 2^+55.4 | −13.1 | 31.2 | 2^61.3 | +41 | – |
| 599, c = 6, q = 11, s = 2, H' = Z9, −1 on Q's Z3 | 36 | 2^85.7; 2^+54.8 | −13.8 | 31.9 | 2^62.0 | +40 | – |
| 719, c = 6, q = 7, s = 2, H' = Z17, ×4 on Z17 | 30 | 2^50.3; 2^+25.4 | −37.0 | 35.1 | 2^71.5 | −13 | – |
| 599, c = 6, q = 11, s = 2, H' = Z3², U = 1 | 54 | 2^116.1; 2^+78.5 | −9.0 | 44.5 | 2^84.0 | +68 | – |

- **Every rich configuration costs more than 2^47 per core.** The cheapest, 599 with ×4 on Z9, is 2.9 to 5.6 bits over
  at the six typical y's. The spread of the cost over typical y's is 2.7 bits there and 4.7 bits for the Z3² slice; the
  two isomorphic Z3² keys of the skipped list (one typical y each: 2^56.5 and 2^60.4) differ by 3.9 bits, inside that
  spread.
- **Within the 599 c = 6, s = 2 family, stage-1 richness and the stage-2 deficit rise together**, and no member falls
  in the window "rich and under 2^47":

  | U (and H') | G'-orbits | stage 1 per class | D_def | cost per core |
  |---|---|---|---|---|
  | −1 on Q's Z3 and on H' (Z9 or Z3²) | 28 | 2^−40 | 17.2 | 2^34.3 |
  | −1 on Z9 (the lead) | 30 | 2^−6.0 | 21.5 | 2^43.5 |
  | ×4 on Z9 | 30 | 2^+30.8 | 24.3–27.5 | 2^49.9–2^52.6 |
  | −1 on one Z3 of Z3² | 36 | 2^+29.5 | 29.7–30.9 | 2^56.5–2^61.2 |
  | −1 on Q's Z3 (Z9 or Z3²) | 36 | 2^+55 | 31.2–31.9 | 2^61.3–2^62.0 |
  | 1 | 54 | 2^+78.5 | 44.5 | 2^84 |

  The step from −1 on Z9 to ×4 on Z9 (same 30 orbits) gains 37 bits of stage 1 per class and costs 8 bits of stage 2.
- **Cores are expected in the 599 rich slices** (+15 to +41 per symmetry class [E]); what is missing is a stage-2 solver
  3 to 6 bits cheaper than the model (§8.3). The 719 slice with ×4 on Z17 is the opposite case: a rich stage 1 but a stage
  2 of −37 bits per y, 2^−13 cores per class.
- **Everything else that is rich costs 2^84 or more**: 599 with U = 1 (2^84), the other 719 slices (2^95 and up),
  c = 8 (839, 823, 919: 2^116 and up), the s ≥ 3 slices (2^116 and up) and the 643 c = 4 slices (2^167 and up).

---

## 7. Sampling

The prompt's step 3 (sample bordered stage-1 solutions of every rich configuration under the budget, screen each y
exactly, rank by stage-2 yield) has no input: no configuration is both rich (≥ 2^+4 per symmetry class) and under
2^47 per core.
- The one configuration under the budget with a non-negligible stage 1, the 719 lead, was decided without sampling:
  the type-level screens of §5 cover every stage-1 solution of every compression quadruple at once.
- On that slice both samplers failed anyway (§5.2): `theory/ab_s1_sample.py` runs out of memory on its α-lists, and
  CP-SAT with the compression fixed finds nothing in 800 s. A sampler for these sparse letter sets would
  enumerate the ten large entries first (their positions are pinned down by the coset sums) and then solve for the
  ±1 part.
- For the cheapest rich-but-expensive slice of §6 (×4 on Z9) the exact screens at the cheap levels leave all 85
  multisets alive, so a stage-2 campaign there would first need the joint levels (§8) and a stage-1 sampler.

---

## 8. Obstructions and the remaining b ≥ 3 frontier

### 8.1 Exact obstructions found here

- **719, c = 6, q = 17, s = 2, H' = Z7, −1 on Q's Z3: no core [M]** (§5), by two exact methods (per type
  quadruple and aggregated per compression quadruple). This was the only skipped configuration with a stage-2 budget
  near 2^48 and a stage 1 that is not negligible.
- **719, c = 6, q = 17, s = 2, H' = Z7, U = 1: 85 of the 89 compression multisets have no lift [M]** (aggregated
  screen, 8 at Z17 and 77 at Z2 × Z17; `logs/s1z/agg/agg_719_c6_q17_s2_H7_U[].jsonl`). The other four pass Z17,
  Z2 × Z17 and Z3 × Z17, and the joint Q × Z17 level (a product of 2.4·10^10 states) or the per-type screen (about
  600,000 type quadruples) was not run: the slice is already dead in stage 1 (2^−9.8 per symmetry class [E]) and costs
  2^60.7 per core.
- **739, c = 4, q = 23, s = 2, H' of order 8, U = 1: no stage-1 solution [P]** (§4.1): the energy forces five large
  letters, the Q-level forces an even number in every coset.
- **919, c = 8, s = 3 (U = 1 and −1 on Z3): no stage-1 solution [P]**: every letter of the f = 6 fibre is at least 5
  in absolute value, and 4·Σ_j w_j·min y² exceeds T1(0) (`logs/s1z/budget_skipped.jsonl`, `energy_ok: false`).
- Configurations with no compression quadruple meeting the Q-level and lag-0 equations: stage 1 empty [M] (§3).
- Calibration rungs, exact (§2): stage 1 empty for c = 6, H' = Z3 at v = 203, 239 (−1 on H') and v = 239, 311 (−1 on
  Q's Z3); c = 8, H' = Z3 (−1) at v = 175, 271, 319, 463; c = 6, H' = Z5 at v = 335, 395 (−1 on H') and v = 395 (−1
  on Q's Z3 and H'); c = 4, H' = Z4 at v = 211; c = 6 at v = 131 (−1 on H' = Z3).

### 8.2 What the exact screens can and cannot do

- The type-level and aggregated screens decide a configuration without any stage-1 solution when the blocks'
  compressed distributions are small: sparse letter sets (s = 2 with large f, as at 719) and few G'-orbits per coset.
- For the 599 c = 6 slices (letters ±1, ±9, ±11) every multiset survives the Z11 and Z2 × Z11 levels; the
  Z3 × Z11 level needs a product of three coset-group sets (1.4·10^9 states per block for U = ×4 on Z9) and the full
  Q × Z11 level six, and the per-type screen faces about 2·10^10 type quadruples (counted for the Z3² variant with
  −1 on Q's Z3 and H'). These are undecided here.

### 8.3 The remaining b ≥ 3 frontier

1. **Rich but expensive, inside the 112 (§6).** The best is 599, c = 6, q = 11, s = 2, H' = Z9, U = ×4 (order 3) on
   Z9: stage 1 2^+30.8 per symmetry class [E], at least 17 bordered multisets, stage 2 2^49.9–2^52.6 per core
   over six typical y's [E], and about 2^17 cores per symmetry class expected [E]. It needs (a) a stage-1 sampler (the
   α-list sampler and CP-SAT both failed on the similar 719 slice), (b) the joint Z3 × Z11 and Q × Z11 screens on the
   sampled y's (cheap per y, unlike per compression quadruple), and (c) a stage-2 solver with a nested middle level
   (`report/b599_answer.md` §5) to bring 2^50 toward 2^45. The Z3² variants are 6–10 bits dearer.
2. **Configurations outside the batch generators.** `theory/ab_cbatch.py` enumerated split G with one cyclic
   multiplier group on one axis of H' (and, for c = 6, −1 on Q's Z3, alone or with H'). Not covered: products of
   multiplier groups on several axes, multipliers acting diagonally across axes (479's ×13), Q-automorphism slices
   beyond `theory/ab_qaut.py`, and the other abelian types of H' at 739 (only H' = Z8 was in the skipped list; the
   proof of §4.1 covers every H' of order 8 with U = 1).
3. **More than 60 G'-orbits.** The largest skipped configurations (643 with 128 orbits, 839 with 64, 719 with 102)
   are already decided by their stage-2 budgets (2^147 and up [E]); configurations with more orbits at the other open
   orders were never listed.
4. **c ≥ 10.** For 669 and 789 (c = 10 only) and the c = 12 options of 599, 659, 719, 839, 863, 947, 971, frames exist
   (`report/alpoge_design_answer.md` §3.1) but no frame pair with a flat forced corner is known; the Singer ansatz
   with an order-9 symmetry is exhausted at c = 10, and the c = 12 probe covered 250 of 75,573 first blocks.
5. **Non-split H and noncyclic Q** (359's Q = Z2³ with H = 2Z8 × Z11, and its real-character row twist): the staged
   systems and the screens here assume G' = Q × H'; the non-split case needs the coset map of `theory/b599_cond.py`
   (which supports it for the Q-level) carried through the type screens.

---

## 9. Reproduce

```
python3 theory/s1z_budget_batch.py logs/s1z/budget_skipped.jsonl --part 0/2          # stage-2 budgets (two workers)
python3 theory/s1z_screen.py logs/s1z/zscreen_skipped.jsonl --steal --N 1500 --K 24    # z-screen (one per core; the first
                                                                                      # 49 rows were run with --N 3000 --K 40)
python3 theory/s1z_table.py --md TABLE.md                                             # §3
python3 theory/s1z_calib.py logs/ab/cfg_347_s6.json OUT.jsonl --N 2000 --exact        # 347 calibration (1 min)
python3 theory/s1z_calib_batch.py                                                     # the 27 rungs (logs/s1z/calib/rungs.json)
python3 theory/s1z_calib_report.py                                                    # §2 summary, K-bootstrap
python3 theory/s1z_border.py logs/s1z/c719/cfg_719_q17_s2_negQ3.json OUT.jsonl        # bordered multisets
python3 theory/s1z_typescreen.py logs/s1z/c719/cfg_719_q17_s2_negQ3.json OUT.jsonl --part i/4   # §5.3 (about 1 h on 4 cores)
python3 theory/s1z_zqscreen.py logs/s1z/c719/cfg_719_q17_s2_negQ3.json OUT.jsonl --levels '[[],[0],[1],[0,1]]'   # cross-check (40 min)
python3 theory/s1z_agg_batch.py --part 0/1                                            # aggregated screens, second tier
python3 theory/s1z_char_batch.py OUT.jsonl                                            # per-character Z_q screens
python3 theory/s1z_cpsat.py CFG.json '<z>' OUT.jsonl --want 3 --secs 200             # CP-SAT stage-1 sampler
python3 logs/s1z/lattice_919_probe.py                                                 # §4.3
```
The validations of §1, §2.1 and §5.3 are `python3 theory/s1z_validate.py canon|sd|invar|screens`
(logs in `logs/s1z/validate/`): canonical counts of the v = 35, 47 controls; direct against S/D counter on all 27
multisets of 347; invariance over random members of four multisets; the 347-core identity at all levels, planted
lifts at v = 131, 95, 83, and the v = 47 positive control.
