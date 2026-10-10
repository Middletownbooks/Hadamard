# Tractability after H(2292): what changed for the open orders below 3000 (including 2092)

Date: 2026-10-10. Inputs:
- `report/RESULT_2292.md`;
- the Paley design session (`report/paley_design_answer.md`);
- the 571 session (`report/paley_571_answer.md`);
- the all-subgroup slice sweeps (`logs/msweep/`, `theory/p571_msweep.py`) run here for every open order.

Tags as elsewhere: [M] exact, [P] proved, [E] Gaussian estimate.

## 1. What H(2292) and the two design sessions changed

1. **The working recipe is a multiplier slice plus an exact split.** Every success has the same shape:
   - the multiplier group fixes a cofactor of G and acts on the rest through cyclotomic classes;
   - the core then splits exactly: stage 1 = the coset sums y (a small integer GS problem on the fixed part), stage 2 =
     the class values given y;
   - each stage is solved by an exact 2+2 join with class splits from pointwise-fixed factors.

   Instances: 359 and 419 (the authors), 573 (ours), and the controls H(164), H(172) and H(316).
2. **Tractability is decided by the split, not by total constraint bits.**
   - 573: stage 2 at about 2^−12 per y was enough; 1 core in 141 y's.
   - 571: stage-2 yields of 2^−39 or worse per y made it unreachable at p = 3 and 5 with the old split.
   - The 571 session's corrected p = 3 split (ternary Eisenstein stage 2) moves 571 ×27 to 2^−3 … 2^−13 per y. Its
     calibrated cost is still 2^52–2^56 per core, because no pointwise-fixed factor is left for a class split.
3. **Gaussian estimates are unreliable in both directions, so exact screens come first.**
   - Exact kills of slices with positive Gaussians:
     - Theorem 5 / the w-screen kills every Paley(q) slice with large q (599 at q = 23, 659 at q = 47, 739 at q = 41,
       523 at q = 29, 537 at q = 67);
     - the θ₅ criterion kills 211 ×27 (+70);
     - the Paley(19) slice of 571 (+39.9) is 94 % refuted by CP-SAT.
   - In the other direction, 573's stage 2 completed far more often than predicted.
   - New here: `scripts/theta_screen.py`, a slice-free exact screen. For each odd prime q | v − 1 it decides whether
     any slice whose multiplier group maps onto (Z/q)* can exist.
4. **Joint local search is not a solver for these systems [M].** Annealing and tabu on the 419 orbit-label problem fail
   (`scripts/slice_sa.py`), so the solving must come from exact joins and splits.

## 2. All-subgroup sweep of the remaining open orders below 3000 (b = 1) [E] + exact screens

- One run per order (`theory/p571_msweep.py`), covering subgroups with 16–90 orbits per block.
- Each top slice is classified by which factor of Z_{v−1} it fixes pointwise and the cyclotomic index s it induces on
  the remaining prime factor.

| v (order) | best slices (gen; orbits; log₂E) | structure | status |
|---|---|---|---|
| 523 (2092) | ×181; 90; **+92.1** | fixes Z₁₈; **quartic (s = 4) classes on Z₂₉** | **alive; new candidate (§3)** |
| | ×91; 54; +11.5 | fixes Z₁₈; Paley(29) (s = 2) | dead (w-screen, q = 29) [M] |
| | ×53; 50; +8.3 | fixes Z₂; index 12 on Z₂₆₁ | unexplored |
| 537 (2148) | ×9; 56; +59.1 | fixes Z₈; **sextic (s = 6) on Z₆₇** | **stage 1 exactly empty** [M] (§3) |
| | ×59 / ×15; 35; +18.5 / +11.5 | fixes Z₂; index 6 on Z₂₆₈ | unexplored (35 orbits: joinable size) |
| | ≤ 30-orbit slices | | exhaustively empty (theory-arm session) |
| 599 (2396) | ×27; 78; +27.3 | fixes Z₂₆; Paley(23) | dead (Theorem 5) |
| 643 (2572) | ×13; 18; −1.8 | | nothing positive |
| 647 (2588) | ×353 / ×183; 70 / 64; +36.2 / +24.6 | fixes Z₂; index 24 on Z₃₂₃ = Z₁₇ × Z₁₉ | unexplored (no fixed factor beyond Z₂) |
| 659 (2636) | ×71; 42; +2.9 | fixes Z₁₄; Paley(47) | dead (w-screen) |
| 669 (2676), 719 (2876) | none in the 16–90 orbit range | | — |
| 739 (2956) | ×215; 90; +61.9 | fixes Z₂; index 24 on Z₃₆₉ | unexplored |
| | ×127; 90; **+57.3** | fixes Z₁₈; **quartic on Z₄₁** | **alive; candidate** |
| | ×289; 54; +1.5 | fixes Z₁₈; Paley(41) | dead (w-screen) |

**θ_q screen** (`scripts/theta_screen.py`): every large-prime q fails, which kills all slices that map onto (Z/q)* for
q = 29 (523), 67 (537), 13 and 23 (599), 107 (643), 17 and 19 (647), 47 (659), 167 (669), 359 (719) and 41 (739).
The quartic and sextic slices above do not map onto (Z/q)*, so the screen does not apply to them. As a check, it
passes 573 and 419: their Paley slices use only the squares mod 11.

## 3. Generalised-Paley (index s ≥ 3) slices: the new candidates

**Ansatz.**
- G = Z_F × Z_q, with Z_F fixed pointwise and Z_q split into s cyclotomic classes of size f = (q − 1)/s.
- Per coset and block there are 1 + s signs: α at z = 0 and c_j on class C_j.
- Stage 1 is y = α + f·Σc_j on Z_F with Σ_k PAF(y_k) = T_J. This is the same exact identity as in the Paley case,
  because it only uses the Z_q-trivial character.

**Stage split** (`scripts/cyclo_split.py`, Gaussian; total from the sweep):

| slice | stage-1 solutions (log₂) | stage-2 completions per y (log₂) |
|---|---|---|
| 523 ×181 (Z₁₈ × quartic Z₂₉) | **+44.8** | **+47.3** |
| 739 ×127 (Z₁₈ × quartic Z₄₁) | +28.7 | +28.6 |
| 537 ×9 (Z₈ × sextic Z₆₇) | +7.2 | +51.9 |
| *control* 523 ×91 (Paley(29), exactly dead) | −41.5 | — |
| *control* 571 ⟨61⟩ (Paley(19), 94 % refuted) | −18.9 | +58.8 |

- The two controls come out negative at stage 1, consistent with their exact refutations.
- **The index-4 slices are in the opposite regime to 571.**
  - Stage 2 is rich: every y is expected to complete many times. So random sublists of the stage-2 lists suffice, and
    the join cost per core is roughly the square root of the list product divided by the yield.
  - The burden moves to finding stage-1 solutions on a small fixed group: Z₁₈, with 10 letters per position.

**Exact checks [M]:**
- **537 ×9 has no stage-1 solution** (`scripts/cyclo_stage1_exact.py 8 67 6`, 20 s).
  - Lag 0 forces either one |y| = 43 and 31 × |y| = 1, or three 21s, one 23 and 28 × |y| = 1. No placement meets the
    off-zero lags.
  - Positive control: restricted to the lag-0 multiset of the verified H(164) core (v = 41, Paley(5) on Z₈), the same
    enumerator returns 1.5·10⁶ verified stage-1 solutions.
  - So the Gaussian +7.2 was again an overcount; 537 ×9 is dead.
- **523 ×181:**
  - lag 0 gives 3n₁₃ + 4n₁₅ + 13n₂₇ + 15n₂₉ = 34, with many admissible compositions (not rigid);
  - the stage-1 annealer (`scripts/cyclo_stage1_sa.py 18 29 4`) did not find a solution in 6 × 5·10⁷ steps (best cost
    784). This is not informative, since annealing also fails the 419 control.
  - Stage 1 here is 4 × 18 letters, which exact methods should handle (§4).

## 4. Implications for 2092 (v = 523)

- **Current route:** the b = 0 d = 9 k-tree on Kaggle, at n = 142 bits (flat-restricted about 91). There have been
  many runs and no hit.
- **New route [E]:** the b = 1 quartic-Paley(29) slice on Z₁₈ × Z₂₉ (generator ×181).
  - Gaussian whole-slice +92, stage 1 +44.8, stage 2 +47.3 per y.
  - The b = 1 border always exists.
  - Its stage 1 is a small exact problem: four length-18 sequences over 10 letters with 10 lag equations. It has the
    same structure as 537's stage 1, which was decided in 20 s, but is larger.
  - **First computation:** an exact stage-1 solver for Z₁₈.
    1. Enumerate the lag-0 compositions.
    2. Build per-block lists by composition (about 2^30 per block for the larger ones).
    3. Join 2+2 with class splits from the rational characters of Z₁₈: compressions to Z₂, Z₃, Z₆ and Z₉, all pointwise
       fixed, as in the 573 Z₄ split.
    4. Alternatively, run CP-SAT per composition. It decides rigid cases fast but may not find generic solutions.
  - **If stage 1 is non-empty, stage 2 is the easy side** (+47 per y). It needs a join over class values given y, with
    per-coset choices C(4, ·) under the fixed sum, class-split by Z₁₈ characters.
  - **Main risk:** Gaussian overcounts at stage 1. They happened at 537 ×9 (+7.2, empty) and in all the Paley(q) cases.
    523's stage-1 margin (+44.8) is far larger than those (+7.2, −18.9, −41.5), but only an exact solver settles it.
- **739 ×127** is the same shape (Z₁₈ × quartic Z₄₁; stage 1 +28.7, stage 2 +28.6) and would reuse the same solver.

## 5. Still unexplored in the sweep

- The index-6/12/24 slices on composite cofactors: 537 ×59 (35 orbits), 647 ×353 (+36), 739 ×215 (+62), 523 ×53.
- These need a stage split over a composite Z_q (several Gauss-period fields). The 571 session's level machinery
  covers part of it.

## 6. 523 quartic slice: exact stage 1 found; stage 2 needs a structured solver (2026-10-10) [M/E]

- **Stage 1 is non-empty [M].** `scripts/cyclo_stage1_cpsat.py 18 29 4` finds a verified y on Z₁₈ in 323 s
  (`logs/cyclo/stage1_cpsat_523.json`): Σ PAF = 1976 at lag 0 and −116 elsewhere, rows (2, 0, 0, 0). The letters
  used are 62 × |1|, 6 × |13| and 4 × |15|.
  - Controls: the same model finds a v = 41 solution in 3.8 s, and proves the 537 sextic slice INFEASIBLE in 3.5 s.
    The latter matches the exhaustive enumerator.
- **Stage 2 model** (`scripts/cyclo_stage2_cpsat.py`):
  - The exact reduction is R(h, t) = Σ_k Σ_g [α_g c_{g+h,t} + c_{g,t′} α_{g+h} + Σ_ij N_t(i, j) c_{g,i} c_{g+h,j}] = −4,
    with N the order-4 cyclotomic numbers of Z₂₉ and t′ the class of −s. That is 36 independent equations; the lags
    (h, 0) follow from stage 1.
  - Control: given the y of the verified H(164) core, it recovers a full core (FFT-checked) in 0.1 s.
- **On 523's y: UNKNOWN after 1800 s** (`logs/cyclo/stage2_cpsat_523_y0.json`).
  - A planted instance of the same size (random α, c; targets computed from it) is also UNKNOWN after 600 s, while a
    planted v = 41 instance is solved instantly (`logs/cyclo/stage2_cpsat_controls.jsonl`).
  - So CP-SAT cannot solve stage 2 at 523 size, and the UNKNOWN on the real y carries no information.
- **Stage-2 size for this y [E]:**
  - 6 choices per |y| = 1 coset and 4 per |y| = 13 or 15 coset, so H₂ ≈ 180 bits;
  - mean yield +47 per y (Gaussian), so C₂ ≈ 133 constraint bits;
  - plain meet in the middle costs about 2^(C₂/2) ≈ 2^66 per y. Not tractable as is.
- **What is needed:** a hierarchical exact stage-2 solver.
  - Z₁₈ = Z₂ × Z₉ is fixed pointwise, so the condition splits exactly into an S level (x₀ + x₁) and a D level
    (x₀ − x₁) on Z₉ × Z₂₉. Entries lie in {0, ±2}, with exactly one of S and D nonzero at each point.
  - Z₉'s rational characters (Z₃-compression) give further class splits.
  - The 571 session's S/D machinery (`theory/p571_hsolve.py`) is the template. Its cost floor is about
    (C₂ − h_D)/2, where h_D is the D-level freedom. With 2^47 completions per y, the slack is very large, and sublists
    can be chosen freely, which the 571 case lacked.
  - Costing this split exactly is the next design step.
