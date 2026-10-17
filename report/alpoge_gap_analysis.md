# What Alpöge's examples say our methods still lack

Date: 2026-10-17. Branch `claude/focused-planck-j1645f`. Answers: "what our methods still lack, which further analysis
of Alpöge's examples might find additional methods for, now that we have successful methods of our own for
comparison"; in particular b = 0 (where his method beats ours), b ≥ 3 (no new example of ours), and a revisit of the
167/179 half-shift template now that we have working b = 1 methods.
Scripts: `forensics/near_structure.py`, `forensics/b0_fingerprints.py`, `forensics/template_opt.py`,
`forensics/template_calib.py`. Logs: `logs/gap/`.

**Tags.** [M] exact computation (script named); [P] proved here; [E] Gaussian estimate (calibration named);
[I] interpretation.

Our methods that built a matrix: Paley(11) + multiplier staged join (H(2292); also 419 from scratch), quartic-Paley
two-level S/D solver (H(2092)), doubly-staged cross-level solver (H(2588)), flat-pool k-tree (H(3532); also two new
487 cores), 3-cycle-twisted b = 0 slice (H(3676)), b = 0 GP dual level (H(3732)), cyclotomic T-matrices with
Cooper–Wallis (H(3668)); plus the staged b ≥ 3 pipeline, which reproduced 347's slice (new H(1388)).

---

## 0. Bottom line

1. **No hidden structure is left to find in the published cores [M].** A census of six families of near-structure
   ("template + a few deviations") over all 12 cores, with 10 of our own cores as controls, recovers every known
   structure: the half-shift arcs of 167 and 179, the Legendre(11) twists of 359 and 419, 479's GL(2,2) twist, and,
   in the controls, our Paley twists and 919's block 3-cycle. It finds **nothing new** in the five b = 0 cores (223,
   283, 311, 487, 491), in 347, or in 443 (§1). Whatever made those orders tractable for Alpöge is in the solver,
   not in a slice we have not seen.
2. **b = 0: his solver is the same family as our flat-pool k-tree, run at a larger scale [M, E] (§2).**
   - Every published b = 0 block is a local minimum of its own PAF energy under row-sum-preserving swaps (0 improving
     swaps in 19 of 20 blocks; 1 in the last). So are the blocks of our k-tree cores (pool blocks come from per-block
     descent). Our cores built without pools (919, 933) have 6–84 improving swaps per block.
   - At 487 our k-tree cores and Alpöge's are indistinguishable (mean PAF²/v 0.11–0.13 against 0.12–0.14).
   - At 223, 283 and 311 his blocks are much deeper (mean PAF²/v 0.082–0.105; 45–68% of lag classes at the least
     possible |PAF|), and the slices are larger: 157, 207 and 139 constraint bits, against 131 for 487 and 112 for
     our 883. The flat-restricted cost at 283 is about 119 bits, which a 4-list k-tree could only meet with lists of
     about 2^40 (ours: 2^24–2^28).
   - **What we lack is scale and depth of the pools, not a new idea.** The cheapest calibration of how far the
     method goes is 311 (flat cost ≈ 84 bits [E], 487's level), then 223 (≈ 90, where our 523 runs failed).
3. **b ≥ 3: what we lack is a stage 1 for larger configurations [M] (§4).** Alpöge's three b ≥ 3 cores are small
   configurations (10–20 G'-orbits). At the open orders every configuration up to 26 G'-orbits that our exact
   batches finished is empty, and **117 of the 270 configurations (more than 26 orbits) were never examined**,
   because exhaustive stage-1 enumeration does not scale there. (The one larger configuration costed in detail, 599
   at c = 6, has essentially no stage 1 once it is counted conditioned on the Q-compressions: `report/b599_answer.md`.)
4. **The half-shift template, revisited [P, M, E] (§3).**
   - It is our S/D split along the Z₂ axis with the S/D support pattern σ fixed in advance, so the two halves share
     no variable [P]. Our b = 1 solvers (523, 647) use the same split but solve S first and pay the D-deficit.
   - One profile controls both halves: W_t = V_t + 2 at every lag [P]. Parity forces every V_t even, which exact arcs
     violate; the published σ is arcs plus a few parity-fixing exceptions and is a local optimum of the Gaussian.
   - **Calibrated against exact counts** (h = 11–23), the Gaussian overstates the S-half by 2.4–3.8 bits and the
     D-half by 0.4–1.2 [M]. With that correction the halves cross zero between v = 359 and 419 and are at −22.6 and
     −15.3 bits at 523, worse beyond [E]. **The earlier verdict stands: not a route to the open orders.**
   - It does not combine with our slices: arcs need a long pointwise-fixed cyclic factor, and our successful b = 1
     slices fix only Z₁₈ (523), Z₄ (573) or Z₂ (419, 647).
   - The halves were not solved with flat pools (their blocks are not local minima) [M]. The natural tool is the
     classical aperiodic complementary-sequence search (lengths 41–46 here), which we do not have; but by the scaling
     above it would only matter below v ≈ 420.

---

## 1. Near-structure census [M] (`forensics/near_structure.py`, `logs/gap/near_structure.{log,json}`)

Each core in its own group, about the common fixed point of its multiplier stabiliser U. Families:
(a) near-(anti-)invariance under every multiplier class a mod U (any order; −1 is the near-symmetric/skew test for
b = 0); (b) half-shift: boundary of σ(g) = x(g)x(g + τ) for every translation τ of order 2; (c) b = 0 only: boundary
of the orbit-label function and of σ_s on the cyclotomic classes (log coordinate), every shift s; (d) cross-block
relations x_k(aO) vs x_j(O), every ordered pair and class; (e) b = 0 only: T-product and Legendre agreement.
Null: orbit labels permuted within orbit-size classes (U-invariance and label counts kept), 300 permutations;
z = (observed − mean)/sd; "fam" = tests in the family (|z| = 4 is about 6·10⁻⁵ per test).

| core | (a) z (fam) | (b) z | (c) z (fam) | (d) z (fam) | reading |
|---|---|---|---|---|---|
| 167 b=1 | −1.9 (81) | **−16.1** | – | 2.9 (984) | half-shift arcs (known) |
| 179 b=1 | 2.7 (87) | **−18.5** | – | 2.9 (1056) | half-shift arcs (known) |
| 223 b=0 | 2.9 (73) | – | −2.4 (73) | 3.5 (888) | nothing |
| 283 b=0 | −3.5 (93) | – | −2.7 (93) | 4.0 (1128) | nothing (Bonferroni p ≈ 0.04, 0.07) |
| 311 b=0 | −2.9 (61) | – | −2.9 (61) | −3.4 (744) | nothing |
| 347 b=5 | −2.8 (5) | −0.7 | – | 2.6 (72) | nothing |
| 359 b=7 | **−8.9** (7) | −5.2 | – | 3.9 (96) | Legendre(11) twist (known); (b) is its consequence |
| 419 b=1 | **−4.9** (11) | 0.5 | – | 3.4 (144) | Legendre(11) twist (known) |
| 443 b=1 | −2.6 (15) | 1.6 | – | −3.1 (192) | nothing |
| 479 b=3 | 4.1 (23) | 0.1 | – | −4.3 (288) | (a): GL(2,2) twist (A and D agree 32/32 on the trivial Z₂²-coset); (d) A←D, a = 29: weak, unexplained (Bonferroni p ≈ 0.005) |
| 487 b=0 | 2.8 (53) | – | −2.3 (53) | −3.1 (648) | nothing |
| 491 b=0 | 1.8 (69) | – | −2.2 (69) | 3.3 (840) | nothing |

(e) at b = 0: |z| ≤ 2.4 everywhere. Controls: our 573 and 419 Paley cores show their twist in (a) (−6.5, −5.2), our
919 its block 3-cycle in (d) (+7.9) and in (a) (+4.9: A alone is invariant under the order-51 group); every other
control family stays within |z| ≤ 4.5, the noise level of families this size (twist structure inflates (d) slightly,
because whole cosets move together). Djoković's 463 core gives (a) +11.2 as an artefact: its blocks are invariant
about different points, so the common-fixed-point stabiliser is trivial and its d = 21 symmetry reappears as a class.

---

## 2. b = 0: same family as ours, at larger scale [M, E]

(`forensics/b0_fingerprints.py`, `logs/gap/b0_fingerprints.txt`)

| core | mean PAF²/v per block | lag classes at least \|PAF\| | improving swaps per block | constraint bits n (own slice) |
|---|---|---|---|---|
| Alpöge 223 (d = 3) | 0.091–0.103 | 0.54–0.68 | 0, 0, 0, 0 | 156.6 |
| Alpöge 283 (d = 3) | 0.082–0.090 | 0.51–0.64 | 0, 0, 0, 0 | 207.0 |
| Alpöge 311 (d = 5) | 0.094–0.105 | 0.45–0.65 | 0, 0, 0, 0 | 139.4 |
| Alpöge 487 (d = 9) | 0.124–0.138 | 0.33–0.52 | 1, 0, 0, 0 | 130.9 |
| Alpöge 491 (d = 7) | 0.097–0.127 | 0.43–0.46 | 0, 0, 0, 0 | 169.1 |
| ours 487 k-tree | 0.112–0.133 | 0.19–0.44 | 0, 0, 0, 0 | 130.9 |
| ours 883 k-tree (d = 21) | 0.146–0.180 | 0.19–0.29 | 1, 1, 0, 1 | 111.9 |
| ours 919 (twisted, no pools) | 0.34–0.40 | 0.11 | 70–84 | 49.6 |
| ours 933 (GP, no pools) | 0.33–0.64 | 0.19–0.38 | 6–43 | – |

- **Local minima [M].** A swap flips one +1 and one −1 U-orbit (row sum kept). Blocks drawn from pools built by
  per-block descent on the PAF energy are local minima; blocks of solutions found otherwise are not (919, 933).
  Alpöge's b = 0 blocks are local minima, like ours. Together with flatness (the W0 fingerprint of `REPORT.md`) and
  the identical 487 statistics, this places his b = 0 method in the flat-pool-and-join family [I].
- **His pools were deeper where the slices were larger.** Mean PAF²/v falls to 0.082–0.105 at 223, 283 and 311
  (ours never below 0.11). Every PAF is ≡ v (mod 4); for v ≡ 3 (mod 4) the least |PAF| values are −1 and 3, and half
  to two thirds of his lag classes sit there [M].
- **No cheaper b = 0 slice exists at these orders [M]** (`scripts/b0_open_screen.py 223,283,311,487,491`): with 15–130
  orbits, each order's own d is its only feasible untwisted slice (491 also has d = 5 at 235.9 bits); the 3-cycle
  twists cost the same (223: 156.4, 283: 206.7, 487: 130.4); 487 d = 27 and 491 d = 35 fail on row sums.
- **Scale [E].** With the measured flat gain of 1.7–1.9 bits per lag class (283, 487), the flat-restricted costs are
  about 84 (311), 90 (223), 106 (491) and 119 (283, measured) bits, against 83 (487) and 78 (883) for our successes
  and about 91 for our failed 523 runs. A 4-list k-tree needs lists of about 2^(n′/3): 2^28 at 311, 2^30 at 223,
  2^35 at 491, 2^40 at 283. `REPORT.md` §0 (item 5) put the whole ρ ≤ 0.30 flat set at 283 near 2^33 per row sum (±5 bits),
  short of 2^40; so either the flat sets are larger than that estimate, or the join was better than a plain k-tree.
- **No join pairing is visible** (pair sums |PAF_A + PAF_B| = 2 in all three pairings alike), unlike our 487 core,
  whose AB|CD pairing shows (15 of 27 classes against 8 and 8) from our best-first τ-buckets [M].

**What to build [I]:** pools of near-floor local minima at 2^30–2^33 per row sum (disk- or GPU-backed), a
measurement of n′(ρ) and flat-set size at 283 that converges, and a calibration ladder 311 → 223 → 491 → 283. Each
step tells how far the method reaches before it is aimed at an open order.

---

## 3. The half-shift template revisited [P, M, E]

### 3.1 It is our S/D split with σ fixed first [P]

For m = 2h (h odd), y ↦ (y mod 2, y mod h) identifies the pair {y, y + h} with the two Z₂-layers at one point of
Z_h. So σ(y) = x(y)x(y + h) is the indicator of the S-cells (x₀ = x₁) of the S/D split along Z₂, and the two
halves of `report/alpoge_methods_forensics.md` §4.2 are the S-level and the D-level of that split. Our 523 and 647
solvers use the same split, but enumerate S-solutions and check D afterwards; the cost floor is the D-deficit
(2^26 S-solutions per core at 523). The template fixes σ before any search, so the halves are independent.

### 3.2 One profile controls both halves [P]

With V_t = Σ_k |S_k ∩ (S_k − t)| and W_t = Σ_k |D_k ∩ (D_k − t)| (cyclic, t ≠ 0): **W_t = V_t + 2.**
Per block, |D ∩ (D − t)| = |D| − |D ∩ (S − t)| = |D| − |S| + |S ∩ (S − t)|; summing, Σ|D_k| − Σ|S_k| = (2h + 1) −
(2h − 1) = 2. Every lag sum of a half has the parity of its term count, so σ is admissible only if every V_t is even
(then every W_t is too) and V_t ≥ 2; summing over t, Σ_k C(|S_k|, 2) must be even. Four exact arcs violate it at
t = 1 (the parity lemma of the forensics report).

### 3.3 What arcs buy [E] (`forensics/template_opt.py`, `logs/gap/template_scaling.txt`)

Gaussian per half: Var(lag sum) = term count (terms are pairwise uncorrelated), lattice spacing 2.

| v | arcs, S / D | random σ (parity ignored), S / D | arcs, calibrated (§3.4) | published σ |
|---|---|---|---|---|
| 167 | +26.0 / +30.4 | +20.1 / +22.8 | +22.6 / +29.5 | +25.8 / +29.2 (parity exact) |
| 179 | +25.9 / +30.3 | +19.3 / +22.1 | +22.5 / +29.4 | +25.7 / +29.2 (parity exact) |
| 251 | +22.3 / +26.8 | +12.0 / +14.7 | +18.9 / +25.9 | |
| 307 | +16.9 / +21.4 | +3.6 / +6.3 | +13.5 / +20.5 | |
| 359 | +10.1 / +14.8 | −5.9 / −3.2 | +6.7 / +13.9 | |
| 419 | +0.7 / +5.4 | −18.5 / −15.8 | −2.7 / +4.5 | |
| 443 | −3.5 / +1.2 | −24.0 / −21.3 | −6.9 / +0.3 | |
| 523 | −19.2 / −14.4 | −44.0 / −41.3 | −22.6 / −15.3 | |
| 571 | −29.7 / −24.9 | −57.1 / −54.4 | −33.1 / −25.8 | |
| 647 | −47.8 / −42.9 | −79.3 / −76.6 | −51.2 / −43.8 | |
| 739 | −71.8 / −66.9 | −108.3 / −105.6 | −75.2 / −67.8 | |

- The aperiodic gain of arcs over random supports grows from about 6 bits per half at 167 to about 26 at 523, but
  the halves still go negative near 400.
- Simulated annealing over σ (sizes free, parity penalised) never improves on the published σ, and from arc or
  random starts it cannot satisfy every parity condition at once (4–9 lags stay odd). The published arcs plus 5
  exceptions are a local optimum, consistent with the parity analysis in the forensics report.

### 3.4 Calibration against exact counts [M] (`forensics/template_calib.py`, `logs/gap/template_calib.json`)

Parity-feasible templates (four balanced arcs, 1–3 flips); exact count of each half by a hashed 2+2 join.

| h (v) | S-half: exact / Gaussian (log₂) | D-half: exact / Gaussian |
|---|---|---|
| 11 (23) | 7.32 / 9.71 (one template: 0 solutions against 9.75) | 12.21 / 12.80 |
| 15 (31) | 8.58 / 12.34 | 14.44 / 15.57 |
| 19 (39) | 11.38 / 14.41 | 16.98 / 17.40 |
| 21 (43) | 12.23 / 15.29 | 17.41 / 18.36 |
| 23 (47) | 12.64 / 16.30 | 18.19 / 19.44 |

The S-half excess (2.4–3.8 bits, not growing) is the per-block row sums: Σ_k r_k² = 1 is implied by the lag
equations, but the independent-lag Gaussian spreads it over one aggregate dimension instead of four discrete row
sums. The scaling table uses −3.4 (S) and −0.9 (D).

### 3.5 Combination with our b = 1 slices [I]

- Arcs exist only along a cyclic coordinate fixed pointwise by the slice group. Our successful slices fix Z₁₈ = Z₂ ×
  Z₉ (523 ×181), Z₄ (573 ×29), Z₂ (419 ×7, 647 ×353): after the Z₂ axis at most a Z₉ remains, too short for an
  aperiodic gain.
- Without arcs, fixing σ first buys parity consistency only, which the S-first search already has (an S-solution
  exists only for a parity-consistent σ), and E_S and E_D move together (§3.2), so S-first already favours σ with
  good D-halves. No gain expected.
- **Where the template would help, it needs a solver we lack:** each half is an aperiodic quadruple problem of length
  41–46 at 167/179 (136–150 constraint bits). Its blocks are not local minima (4–151 improving swaps per block; energy
  40–60% of random signs on the same support) [M], so they did not come from flat pools: a different solver from
  b = 0. The classical tool is end-to-middle search with spectral filters, as for base and Turyn-type sequences.

---

## 4. b ≥ 3: where our coverage stops [M]

- Alpöge's b ≥ 3 configurations are small: 347 has 10 G'-orbits (q = 19) and 479 has 16–20; 359 is a Paley(11) GP
  slice with a Singer-type c = 8 border (`report/alpoge_design_answer.md` §1.4).
- At the open orders our exact batches (`theory/ab_cbatch.py`, ≤ 26 G'-orbits, 1,200 s and 4.5 GB caps) found no
  survivor at c = 6 outside 659 or at c = 8; c = 4 left two f = 2 slices. Of 270 configurations, **117 have more
  than 26 G'-orbits and were never examined**, and 76 hit the caps (§11.2–11.3 of the design answer).
- Small configurations at the open orders are empty, and the larger ones are beyond exhaustive stage 1. The one
  larger configuration costed in detail, 599 at c = 6 (30 G'-orbits), turned out to have essentially no stage-1
  solutions: its "+48 bits" was an artefact of sampling y with the lag-0 energy only; conditioned on the 420
  admissible Q-compression quadruples it holds about 2^4.9 ordered solutions, and its rich neighbour (H' = Z₃², −1 on
  one axis, stage 1 about 2^45.7) costs 2^58–2^61 per core in stage 2 (`report/b599_answer.md`, 2026-10-15).
- **What to build:** first a calibrated stage-1 count per configuration (the z-conditioned estimator of
  `theory/b599_flevels.py`, which corrected 599's +48 to 2^4.9), to find configurations whose stage 1 is rich and
  whose stage 2 is cheap; then sampling there (`theory/ab_s1_cpsat.py --enumerate`, `theory/ab_s1_sample.py` per
  compression quadruple), ranked by stage-2 Gaussian yield, feeding the existing stage-2 solver.

---

## 5. Ranked: what to add

1. **Scale the flat-pool k-tree (b = 0) [§2].** Deeper pools of local minima (mean PAF²/v about 0.09) at 2^30–2^33
   per row sum; calibrate on 311 (≈ 84 flat bits) and 223 (≈ 90) before aiming at an open b = 0 slice. This is the
   one place where Alpöge's examples show a concrete, reproducible method that is simply stronger than ours.
2. **Stage-1 counting and sampling for b ≥ 3 configurations with 30–60 G'-orbits [§4].** First a calibrated
   (z-conditioned) stage-1 count over the unexamined configurations; 599 at c = 6 is already known to be dead.
3. **Template-first decoupling: not recommended [§3].** It is calibrated dead at v ≥ 450, and inside our slices it
   reduces to the S-first search we already run.
4. **Low priority:** the weak A←D relation in 479 (§1).

---

## 6. Reproduce

```
python3 forensics/near_structure.py --P 300 --json logs/gap/near_structure.json      # §1, about 10 min
python3 forensics/b0_fingerprints.py                                                 # §2
python3 scripts/b0_open_screen.py 223,283,311,487,491 20000                          # §2, slice costs
python3 forensics/template_calib.py 11,15,19,21,23 --per-h 4                         # §3.4, seconds
python3 forensics/template_opt.py 167,179 --steps 150000 --seeds 1                   # §3.3 (annealing)
```
The §3.3 scaling table is `logs/gap/template_scaling.txt` (`template_opt.evaluate` on arcs and on 20 random σ per v).
