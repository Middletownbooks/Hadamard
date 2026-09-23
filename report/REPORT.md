# v=283 forensic workup → solver reconstruction: status report

Date: 2026-09-23. All scripts are in `scripts/`, raw logs in `logs/`, and generated tables in `report/tables.md` (regenerate with `python3 scripts/make_tables.py`).
**[M]** = measured in this session. **[E]** = estimate or model output. Each claim carries one of the two tags.

## 0. Headline

1. **[M] Flatness is a family signature, not a universal one.** All five Alpöge-group b=0 cores (v = 223, 283, 311, 487, 491) have individually flat blocks. Their share variance is 0.0053–0.0080 against a null of ≈0.036. The four b=1 cores (167, 179, 419, 443) and the three b≥2 cores (347, 359, 479) look like uniform solutions: share variance 0.030–0.040 and Fisher p ≈ 0.05–0.6. Djoković's 1992 v=463 core (b=0, d=21) sits in between: two blocks are flat and two are not. Per the W0 decision rule, the b=0 family and the b≥1 families were produced by different solvers, or at least different objectives.
2. **[M] Within b=0, flatness weakens as the orbit count r falls.** In r order, share variance goes 0.0053 (r=95), 0.0057 (75), 0.0072 (71), 0.0062 (63), 0.0080 (55). Relative to the single-block flatness floor that tabu reaches, the real blocks sit at 1.15× the floor (283), 1.22× (223), 1.43× (491), 1.45× (311) and 1.56× (487).
3. **[M] Stop rule triggered.** No candidate reached v=151 faster than MITM with the fingerprint reproduced:
   - **S3** (joint tabu with an annealed flatness penalty) is the only candidate that works at all. It reproduces partial flatness at little extra cost. Its cost grows as roughly 2^(0.55·n) over n = 33–50. MITM already beats it by 20–30× at v=131. It found no solution at v=151 in 13 runs of 40M flips.
   - **S2** (residual sharing, two implementations) fails at v=79.
   - **S4** (thin-set RRR, 10 variants) fails at v=61.
   - Breakout weighting, fixed-λ and hinge variants were no better.
   - **W3 was not started.**
4. **[M] Flat-subspace constraint bits (S1).** Partial joins on real flat pools confirm the Gaussian model to within 1–3% per class. The result is n′ ≈ 119 ± 2 bits at v=283 (profile 19,19,19,7, ρ ≤ 0.30) and n′ ≈ 83 bits at v=487 (profile 37,17,17,1, ρ ≤ 0.36). The unrestricted values are n = 207 and 130.
5. **[M/E] Solver-mechanism test: flat pools plus a k-tree (Wagner) 4-list join. Consistent at v=487, unresolved at v=283.**
   - At v=487 the Wang–Landau (WL) flat-block counts converged, and n′(ρ) was measured from pools.
     - Flat solutions in the real profile first exist at ρ ≈ 0.30.
     - A k-tree join over the whole flat set first becomes feasible at ρ ≈ 0.39.
     - The real blocks sit at 0.345–0.369: far above the existence edge, a little below join balance.
   - At v=283 the WL counts below ρ ≈ 0.34 did not converge (4 runs), so the same test is inconclusive there.
   - The rough numbers that exist for v=283 (~2^33 flat blocks per row sum at ρ ≤ 0.30, hence only ~2^10 flat solutions in the profile) would put the real v=283 blocks near the existence edge, which would argue against a join. Treat that as unverified, uncertain by ±5 bits or more.
   - **Practical consequence for W3.** Whatever the group used, a flat-subspace k-tree join at d=9 is the only approach in this workup with an estimated cost within budget: lists of ~2^28–2^29 flat blocks per slot at v=487 and ~2^30 at v=523 [E]. It is not a reconstruction of their method.

## 1. Data provenance and core extraction [M]

- **Source.** The 13 full matrices from `bbeartheancient/hoa64/matrices`: 12 Alpöge orders plus Djoković's order 1852. All are normalised: first row and first column all +1.
- **Method** (`scripts/extract_cyclic.py`, `scripts/extract_group.py`).
  - Each GS block is D_r·C·D_c, where C is circulant or back-circulant, or group-developed f(h−g) / f(g+h).
  - Each block's sequence is recovered by propagating sign ratios, then checked with an exact rank-one test on the whole block.
  - The candidate (−1)^i twists, and the order-2 character twists for groups, are resolved by requiring the complementarity equations to hold.
  - Exactly one valid combination exists for each core. For v=283 the result is identical to `seq_v283_g.txt` (the Appendix A sequence).
- **b≥2 group structure.** The group is recovered by testing mixed-radix translation invariance.
  - v=347: G = Z2×Z3×Z57 (b=5).
  - v=359: G = Z2×Z2×Z8×Z11 (b=7).
  - v=479: G = Z2×Z2×Z119 (b=3).
  - In each case the PAF sum takes the values −4b on H∖{e} and +4 off H, with [G:H] = b+1. The core occupies rows and columns 4b onward.
- **Outputs.** Cores in `data/core_v*_b*.txt`. Every one passes `Core.verify()`.

## 2. W0 — fingerprint census [M]

The full table is in `report/tables.md`, and the raw data in `logs/census.json` (script `scripts/census.py`).

| v | b | d | r | row sums | special | mean PAF std/√m | share var | log10 Fisher p | nbrs (raw ≤2 / orbit ≤2) |
|---|---|---|---|---|---|---|---|---|---|
| 167 | 1 | 1 | 166 | (2,0,0,0) | – | 0.696 | 0.0301 | −1.0 | 0 / – |
| 179 | 1 | 1 | 178 | (2,0,0,0) | – | 0.731 | 0.0331 | −0.6 | 0 / – |
| **223** | 0 | 3 | 75 | (11,11,11,23) | yes | **0.302** | **0.0057** | **−25.9** | 0 / 0 |
| **283** | 0 | 3 | 95 | (19,19,19,7) | yes | **0.290** | **0.0053** | **−31.0** | 0 / 0 |
| **311** | 0 | 5 | 63 | (21,19,21,1) | no | **0.316** | **0.0062** | **−17.9** | 0 / 0 |
| 347 | 5 | 6 | 70 | (6,10,26,24) | no | 0.753 | 0.0354 | −1.0 | 0 / 0 |
| 359 | 7 | 5 | 96 | (14,18,30,4) | no | 0.820 | 0.0401 | −0.3 | 0 / 0 |
| 419 | 1 | 15 | 42 | (2,0,0,0) | – | 0.697 | 0.0319 | −0.7 | 0 / 0 |
| 443 | 1 | 12 | 50 | (2,0,0,0) | – | 0.752 | 0.0355 | −0.2 | 0 / 0 |
| 463 (Djoković) | 0 | 21 | 23 | (43,1,1,1) | yes | 0.577 (A, C 0.44; B, D 0.70) | 0.0219 | −1.0 | 0 / 0 |
| 479 | 3 | 4 | 128 | (6,6,38,20) | no | 0.703 | 0.0309 | −1.3 | 0 / 0 |
| **487** | 0 | 9 | 55 | (17,1,17,37) | no | **0.358** | **0.0080** | **−15.4** | 0 / 0 |
| **491** | 0 | 7 | 71 | (29,27,15,13) | no | **0.339** | **0.0072** | **−17.1** | 0 / 0 |

Notes. These are measurements unless a note says otherwise.

- **d and r.** d is the order of the largest multiplier group U common to all four blocks, after a per-block translation. The r values reproduce the prompt's (75, 95, 63, 55, 71, 42, 50). For the b≥2 cores, U contains only scalar multipliers; GL-twisted automorphisms were not searched.
- **Row sums** are given up to sign and translation convention. For b≥2 cores, "special" is not meaningful.
- **Parseval identity.** For b=0 and b=1, the mean PAF std/√m equals 4·√(share variance). The two fingerprint numbers are therefore one statistic, not two independent confirmations.
- **Pair-share variances** always come in equal pairs: AB≡CD, AC≡BD, AD≡BC, because the four shares sum to 1. The prompt's "AD and BC both 0.010" is this identity, not an independent pairing hint. After that, the v=283 pair pattern (AD=BC=0.010 vs 0.005–0.006) does not repeat across the flat family. At 223 the largest pair is AB/CD; at 487 it is AC/BD. There is no consistent pairing.
- **Rarity metric correction.** The prompt's "product of per-class p-values" (reproduced here: −76.8 at v=283) has null expectation −0.434·(number of classes), which is −20.4 at v=283. It is not a p-value. I report the Fisher combination instead.
  - W1 shows Fisher is somewhat anti-conservative in the tail: P(log10 p ≤ −2) = 0.017–0.023 against a nominal 0.01, and P(log10 p ≤ −3) = 0.001–0.004 against 0.001. Discount the b=0 Fisher values by roughly 1 order of magnitude. They remain decisive.
- **Neighbour counts.** No core has another solution within Hamming radius ≤ 2, whether counting arbitrary entry flips or whole-U-orbit flips. This includes swaps within one block and flips across two blocks. The count is exact, by hashing the ΔPAF vectors.
- **Half-shift structure in 167/179: not reproduced.** I tested:
  - agreement of x_k(i+m/2) with ±x_l(u·i+t)·χ(i), for all block pairs, all units u, all t, and χ ∈ {1, (−1)^i};
  - PAF_k(s) against PAF_l(s+m/2);
  - P_k(ω) against P_l(ω+m/2).

  The best raw result was p ≈ 10^−6.9, before a Bonferroni correction over ~10^6 tests, i.e. nothing. The strongest half-shift effect in the whole set is block A of v=419: PAF_A(209) = +130, a 274/418 self-agreement at half shift, p ≈ 10^−9.7 uncorrected. **I need the definition behind the 10^−55 and 10^−59.5 figures** (which sequence representation, which statistic) to relate it to flatness. The b=1 cores are not flat, so whatever that structure is, it is not the flatness fingerprint.
- **The largest block share is not a filter signature.** In the flat cores, max P_k(ω)/v ≤ 1.97 (share ≤ 0.494). A hard "share ≤ ½" filter alone would only bring the share variance down to 0.0187 (Dirichlet conditioned on it), so it cannot explain 0.005. The observed maxima are what general flatness predicts.

## 3. W1 — null calibration [M]

Uniform solutions were drawn with an unbiased pool-MITM (`scripts/w1_null.py`):
- random uniform pools of 2^11–2^13 blocks per slot, per row-sum profile;
- a linear 64-bit hash of the lag-class PAF vector, so the AB/CD join is an outer sum;
- a uniform subsample over all matches;
- every sampled solution re-verified.

At v=43 the enumeration is exhaustive: 3.3·10^11 ordered solutions, so the exact n = 21.72 bits against the Gaussian model's 21.3 [M]. Summary in `logs/w1_summary.json`.

| v | d | n [E] | #sampled | share var mean (sd) | min | 0.0375(1−1.45d/m) | P(log10 Fisher ≤ −2) | KS D (per-class T vs Dirichlet) |
|---|---|---|---|---|---|---|---|---|
| 43 | 3 | 21.3 (exact 21.7) | 1200 | 0.0337 (0.0100) | 0.0109 | 0.0337 | 0.023 | 0.061 |
| 61 | 3 | 32.9 | 7200 | 0.0351 (0.0082) | 0.0109 | 0.0348 | 0.017 | 0.036 |
| 73 | 3 | 41.3 | 6786 | 0.0353 (0.0076) | 0.0136 | 0.0353 | 0.019 | 0.032 |
| 79 | 3 | 45.4 | 1063 | 0.0356 (0.0071) | 0.0176 | 0.0354 | 0.020 | 0.032 |
| 131 | 5 | 50.5 | 1019 | 0.0353 (0.0073) | 0.0169 | 0.0354 | 0.018 | 0.035 |
| 151 | 5 | 59.6 | 6 | 0.030 (n too small) | – | 0.0355 | – | – |

v = 97, 103 and 109 at d=3 were also attempted. With 2^12 pools and 240 s each they gave no solutions, so the ladder stops at v=79 for d=3.

- **Dirichlet(1,1,1,1) is a good approximation, with a finite-size correction.** The mean share variance is ≈ 0.0375·(1 − 1.45·d/m) [fit, M]. The per-class T statistic is detectably non-Dirichlet (KS D ≈ 0.03–0.06), but the effect is small.
- **Implications for W0.** At the real cores' parameters the null mean is 0.0359–0.0367. Scaling the null sd by √(classes) gives z ≈ 8 at v=283 and z ≈ 6 at v=487 [E]. **The b=0 rarities are trustworthy in direction and are conservative after the Fisher discount.**
- **Unbiased solutions are not flat.** Out of ~17,000 unbiased solutions, the flattest has share variance 0.0109. That is still 2× the value of any real b=0 core. Solutions at v=283 flatness do not occur by chance at small m.
- **n estimator.** `scripts/nbits.py` is validated at v=43 (0.4-bit underestimate). Its values: 283 → 206.8, 223 → 156.6, 311 → 139.6, 487 → 130.1, 491 → 169.2, **523 (d=9) → 141.9**. Two of these differ from the prompt's figures (≈200 at v=283, 53 at v=151); mine is validated only at v=43.

## 4. W2 — solver reconstruction on the ladder (b=0) [M]

Engine: `scripts/engine.py`, numba. It works on orbit-level states with incremental O(L·r) delta and flip updates, verified against FFT. Throughput is ≈3·10^5 flips/s at v=79–151 and ≈10^4–5·10^4 at v=283 per core.

### 4.1 Candidates and outcomes

| candidate | what it is | result |
|---|---|---|
| S3 tabu, λ=0 (control) | tabu on ΣR² | Solves v ≤ 131. Median flips: 9.6k (n=33), 87k (41), 405k (45), 7.1M (50.5). v=151: 0/4 at 40M. Solutions are null-like (share var 0.031–0.036). |
| **S3 annealed λ with a floor** | tabu on ΣR² + λ(t)·Σ_k ΣPAF_k², λ decaying 3 → λ_min | λ_min=0.3: 6/6 solved at 73, 79, 131, with flat-ish solutions (share var 0.018, 0.017, 0.021) at ≤1.2× the λ=0 cost. λ_min=0.6–1.0: flatter (share var 0.012–0.016), but cost rises 2–5× and success falls at 131. v=151: 0/9 at 40M (λ_min 0.3, 0.6). |
| S3 fixed λ | constant λ = 1 or 3 | Mostly fails (R=0 is not a minimum of E). Solved at v=61; 2/3 at v=79 with share var 0.012. |
| S3 hinge | λ·max(0, F_k − F0) | Far fewer solutions (flat solutions are rarer); v=79 solved 1/9. |
| S3 breakout | adaptive class weights | No gain over the control (v=131: 5.5M, fail, 3.2M flips). |
| S2 (two versions) | block coordinate descent on (PAF_k − t_k)², t_k = PAF_k − R/4 (or −R), Jacobi or Gauss–Seidel | Fails at v=79 (0/10 in 2M flips; residual stalls at ΣR² ~ 10³). Restricting moves to one block per turn is much worse than joint moves. |
| S4 RRR | per-block thin sets \|X̂_k\|² = T_k, T_k updated by residual sharing (EMA η ∈ {1, 0.1, 0.01}, β ∈ {0.3, 0.5, 0.7}) | Fails at v=61 (0/12, 20k iterations). The target feedback makes the B-set chase the iterate. |
| S1 | flat pools + join | Measurement only; see 4.3. A generic 2-list join is hopeless at v=283. The 4-list k-tree join is the proposed next experiment (§5). |
| MITM (baseline) | unbiased pool-MITM, same code as W1 | v=131: 1,220 solutions / 1,800 s ≈ **1.5 s per solution**. v=151: 6 solutions / ~2,600 s ≈ **7 min per solution** (shared CPU). |

Scaling plot: `report/w2_scaling.png`. It shows median flips-to-solution against n. Hollow markers mean fewer than half the runs solved at the cap.

### 4.2 Why the stop rule fires

- **S3 is slower than MITM already at v=131.** It takes 25–50 s per solution against ~1.5 s. At v=151 it failed 13/13 runs at 40M flips (~4–5 min each), against ~7 min per solution for MITM.
- **Steps grow with n and with constraint density.** Over the ladder, log2(steps) rises about 0.55 bits per bit of n. The step from v=79 to v=131 has the same orbit count r=27 but higher density δ = n/(4r) (0.42 → 0.47), and costs 17×. The solved cores sit at δ = 0.52–0.60, so extrapolation is strictly worse [E].
- **The fingerprint is only partly reproduced.** Annealed-λ S3 gives share variance 0.012–0.021 at v = 79–131, i.e. PAF ratio ≈ 1.5–2× the single-block floor at that m. That matches the looser real cores (487: 1.56×, 311: 1.45×), not v=283 (1.15×). Nothing reached the "≈0.29·√m" level, and at small m that level may not be reachable relative to the floor.

### 4.3 S1 measurements

- **Flat pools** (`scripts/s1_pools.py`). These are distinct blocks with PAF ratio ≤ ρ, sampled from tabu trajectories. They are not uniform samples of the flat set.
  - v=283, ρ ≤ 0.30: ~5,500 blocks in 15 min on one core, spread over row sums.
  - v=487, ρ ≤ 0.36: ~4·10^6 blocks in 10 min.
- **n′ from partial joins** (`scripts/s1_nprime.py`). Measured bits over k lag classes against the Gaussian prediction from the pools' own covariance:
  - v=283 (19,19,19,7), ρ=0.30, pools 283–765 per row sum, 2^37.7 four-tuples. Measured 2.69 / 5.28 / … / 28.93 bits for k = 1 / 2 / … / 11, against Gaussian 2.63 / 5.23 / … / 28.71. The measured excess is +0.02–0.08 bits per class. **n′ = 118.1 (Gaussian) → ≈119–120 (corrected)** [M/E]. That confirms the prompt's rough 120.
  - v=487 (37,17,17,1), ρ=0.36, pools of 4,096. k = 1–3: 3.25 / 6.51 / 9.77 bits measured, against 3.20 / 6.41 / 9.61 Gaussian. **n′ = 81.6 (Gaussian) → ≈83 (corrected).**
- **Flat-block counts** (Wang–Landau over F = Σ_c PAF(c)², `scripts/s1_flat_measure.py`). The method is validated against exact enumeration at v=61: agreement within 0.2 bits at every threshold.
  - **v=487, d=9.** Converged to ln f = 10⁻⁵ in 62M flips; the coarse 10⁻³ run agrees within ~0.3 bits for ρ ≥ 0.34. Counts are log2 of the number of all 2^55 blocks with ratio ≤ ρ:

    | ρ | 0.30 | 0.32 | 0.34 | 0.36 | 0.38 | 0.40 | 0.44 |
    |---|---|---|---|---|---|---|---|
    | log2 count | 21.5 | 23.6 | 25.8 | 28.1 | 29.9 | 31.7 | 34.5 |

    The row-sum split (≈ −2.3 bits for each of 1, 17, 37) comes from the pools and is approximate.
  - **v=283, d=3: not resolved.**
    - Two full-range runs at 10⁻³ tolerance agree only to ~2 bits at ρ = 0.34–0.36: log2 count ≈ 43–45 at 0.34 and 48–50 at 0.36.
    - Below 0.32 the walker barely visits, and two range-restricted runs failed to converge (`logs/wl_283_combined.log`).
    - The one number available, ≈ 2^35 total at ρ ≤ 0.30, is **unreliable**.
- **Join-balance / existence test** (`scripts/rho_star.py`, `logs/rho_star_*.log`). Definitions:
  - ρ_exist: 4 × log2 N_flat per slot ≥ n′.
  - ρ_join: log2 N_flat per slot ≥ n′/3.

  | v | ρ_exist | ρ_join | real blocks | real / tabu floor |
  |---|---|---|---|---|
  | 487 | ≈ 0.30 | ≈ 0.39 | 0.345–0.369 | 1.56× |
  | 283 | ≲ 0.30 (unreliable) | ≈ 0.36 (±0.01) | 0.286–0.299 | 1.15× |

  - At v=487 the real solution is ~2^19 deep in the flat-solution set; at ρ = 0.36 about 2^19 flat solutions exist in the profile. Its flatness is 1–2 WL bins below join balance.
  - At v=283 the real blocks are well below ρ_join. If the rough count is right, a join at ρ_join would have produced blocks near 0.36, not 0.29.

### 4.4 Long-cap S3 at v=151 [M]

`logs/grid3.jsonl`, 400M-flip cap, on a CPU shared with 3–5 other jobs:

| λ_min | result | flips | wall time | share var |
|---|---|---|---|---|
| 0.0 | solved | 218M | 21 min | 0.040 (null-like) |
| 0.0 | **unsolved** | 400M | 33 min | – |
| 0.3 | solved | 67.5M | 9 min | 0.020 |

Together with 13/13 failures at 40M flips, the S3 time-to-solution at v=151 is ≳15–30 min, against ≈7 min per solution for crude MITM. **The stop rule holds.**

## 5. What I recommend next (not run: the stop rule fired)

1. **Decide whether a flat-subspace join is worth building as a *solver* at d=9.** It does not have to match the group's method.
   - **Balance at v=487.** Measured n′(ρ) is 84–85 at ρ ≈ 0.39–0.40, and the WL count gives ~2^27.4–2^29.1 per row sum there. The k-tree join therefore needs lists of ~2^28.3 per slot, with time and memory of the same order: 2^28 × 8 B ≈ 2 GB of compact keys [E].
   - **Stages.** Stage 1 matches the pair sums on k ≈ 10 lag classes; stage 2 joins the remaining ~17.
   - **Bottleneck: pool generation.** The current tabu sampler yields ~4·10^5 distinct flat blocks per core-minute, so 2^28.3 per row sum would take ~13 core-hours per row sum. That needs a ≥10× faster generator (stride-1 sampling, enumerating the flat neighbours of each visited state) before it becomes a Kaggle job.
   - **Validation.** Ground truth first: at v=487 at least one solution is known to exist at ρ ≤ 0.37, and ~2^19–2^27 flat solutions are expected.
   - **Transfer to v=523 (d=9, L=29, n=141.9 [E]).** Measure n′(ρ) and WL there first; the expectation is lists of ~2^30 [E].
2. **Fix WL at v=283 in the low-flatness tail.** Use replica-exchange WL or multicanonical sampling seeded from the tabu flat pools. That resolves ρ_exist and whether the v=283 solution is "one of the flattest ~2^10".
3. **For the b=1 and b≥2 families,** which carry no flatness fingerprint, the census gives no handle on their solver. Any reconstruction there has to start from different evidence, such as the half-shift statistic, once it is defined.

## 6. What I could not determine / explicit uncertainties

- The method the Alpöge group actually used. Everything in §5 is a hypothesis that fits the census; it is not evidence of their method.
- The half-shift statistic in 167/179 (definition needed).
- The flat-block count at v=283 below ρ ≈ 0.34. Plain WL did not converge there, so the prompt's ~2^38 per slot (and ~2^32 flat solutions) is neither confirmed nor refuted. My rough number (~2^33 per row sum) is unreliable.
- **Pools are trajectory samples, not uniform over flat blocks.** The partial-join/Gaussian agreement suggests the covariance is representative, but a bias of a few percent in n′ cannot be excluded.
- **S3 cost at v=151** rests on 3 long runs (67M, 218M, >400M flips), so the median is only roughly known.
- **The MITM baseline is my own crude numpy implementation.** A tuned MITM would be faster, which only strengthens the stop-rule conclusion.
- **GL(3,2)-twisted automorphisms for the b≥2 groups were not searched.** The d reported for 347, 359 and 479 covers scalar multipliers only.
