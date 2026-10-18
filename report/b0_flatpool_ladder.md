# b = 0 flat-pool k-tree: the calibration ladder (311, 223, 491, 283) and the 859/13 decision

Session of 2026-10-17 (prompt `report/prompt_p1_flatpools.md`; running log at the end of
`report/prompt_flatpool_stage1.md`). Tags: [M] exact computation or measurement, [P] proof, [E] estimate.

## Bottom line

- **859/13: no-go on this machine, and on Kaggle [E].**
  - Its ρ-flat sets hold about 2^−12 cores at ρ ≤ 0.34 and about 2^−4.6 at ρ ≤ 0.36.
  - Cores first exist near ρ ≈ 0.37–0.38, where the flat sets are about 2^33 blocks per slot.
  - A plain 4-list k-tree there needs pools of that size: roughly 1 TB of prep data and 2^45 or more operations per
    core (§6).
- **Rung 3 (491) is stopped by the rule of the prompt [E].**
  - At ρ = 0.34, pools of 2^27 (the disk ceiling here) hold 2^−1.7 solution tuples.
  - The whole flat set (2^30.6 per slot) would need about 2^46.4 operations and roughly 1,500 hours per core.
- **Rungs 1 and 2 reached: new, inequivalent cores at 311 and 223 [M]** (`data/core_v311_b0_d5_0.txt`,
  `data/core_v223_b0_d3_0.txt`). Both were found at the rate the model predicts: buckets 1,118 and 773, against about
  1,200 and 600 per core.
- **Rung 1 (311): the method works at scale; the rate is slow [M/E].**
  - The pools are the whole reachable ρ ≤ 0.32 flat set.
  - A production-scale planted test recovers Alpöge's core.
  - The measured and modelled rate is one core per about 1,200 k = 9 buckets, about 2 hours of join per core on 4 cores.
  - Production result: **a new, inequivalent 311 core** after about 1,100 k = 9 buckets (§3), at the predicted rate.
- **283 (rung 4, measurement) [M/E]:** the ρ ≤ 0.30 flat set is about 2^35 per row sum (converging
  capture–recapture). It holds about 20 cores, and a plain k-tree needs about 2^53 operations per core, so Alpöge's
  283 is reachable by a plain k-tree only at cluster scale (§7).
- **New quantitative understanding [E, checked against measurements]:**
  1. **Energy-shell drift.** Flat blocks sit on a PAF-energy sphere. This makes the true n′ about
     **n′_Gauss + 0.1 · n_ind**: 0.098 bits per independent class, from an exact sphere model with no free parameter
     that reproduces the measured partial-join drift. Every earlier ladder estimate was about 3 bits optimistic.
  2. **Bucket order.** Inside a bucket the rest classes carry extra variance, about 0.08–0.15 bits per rest class.
     The separable best-first score is already the expected-count order; the energy tilt is negligible.
  3. **Clusters.** A core occupies a bucket as a cluster of 32 tuples (single-block reversals × slot swaps).
     Reversal-canonical pools halve the work at no loss.
  4. **Existence edge.** Alpöge's cores sit just above the existence edge of their slices: the ρ at which the whole
     flat set first contains a core.



## 1. What changed in the method

| change | effect | where |
|---|---|---|
| Full rotation closures in all four slots (no symmetry breaking) | each core present in all e rotations (31–37 distinct buckets) | `build_slot_pools*.py --closure` |
| Reversal-canonical slots | one block per reversal pair: pools ÷ 2, pairs ÷ 4, cores per bucket unchanged | `--revcanon` |
| 6-bit (k ≤ 10) and 5-bit (k ≤ 12) packed keys | keyed classes beyond 8 (k = 8 was OOM at 311: 2^31 pairs per bucket) | `kt_join_slots*.py`, `kt_join_disk.py` |
| Partner-key order independent of τ [P] | A and C sorted once, not per bucket: 18 s → 6 s per bucket at 311 | same |
| Disk-capable join | presorted memory-mapped pools, keys computed in the merge, shared files for equal pools | `kt_join_disk.py` |
| Skip histogram merge when one pass fits | one third less pair work | `kt_join_slots.py` (fast mode) |
| Radix-partitioned exact key match instead of the 1 GB byte bloom | 311 k = 9 bucket 18.6 s → 6.5 s (the bloom was ~61% of a bucket, by py-spy) | `scripts/pmatch.py` |

Validation: planted cores recovered by every join variant at 311 (k = 10, 12) and 223 (two-word masks, shared slots),
and at production scale at 311. `kt_join_disk.py` reproduces `kt_join_slots.py`'s bucket counts exactly.

## 2. Models of the rate

**(a) Sphere (energy-shell) drift [E]** (`scripts/sphere_drift.py`).
- Model: four independent uniform points on spheres S^{n−1}(R); the exact density at 0 of their sum (projected on j
  coordinates), against the Gaussian with the same covariance.

| j (n = 30) | 8 | 10 | 12 | 14 | 20 | 30 |
|---|---|---|---|---|---|---|
| sphere model, bits | 0.23 | 0.35 | 0.49 | 0.66 | 1.31 | 2.92 |
| 311 measured, global partial joins | 0.23 | 0.35 | 0.54 (in bucket) | — | — | — |

- The excess over all n coordinates is 0.098 n bits (1.94 at n = 20, 2.92 at 30, 3.51 at 36, 4.50 at 46).
- It is scale-free, so true n′ ≈ n′_Gauss + 0.1 · n_ind for every slice.

**(b) In-bucket rate [M/E]** (`scripts/bucket_pjoin.py`, 311, k = 10).
- Method: the AB and CD pair lists of one bucket, with pair-pairs counted on j more classes.

| bucket \|τ\|² | pairs AB × CD | in-bucket Gaussian, 20 rest classes | E (tuples) |
|---|---|---|---|
| 40 (top) | 2^25.8 × 2^24.6 | 55.88 bits (2.79 per class) | 2^−5.5 Gaussian, ~2^−7 with drift |
| 136 | 2^24.8 × 2^23.6 | 55.54 | 2^−7.1 Gaussian |
| 232 | 2^23.8 × 2^22.6 | 55.19 | 2^−8.7 Gaussian |
| 392 | 2^22.1 × 2^20.8 | 54.54 | 2^−11.6 Gaussian |

- Inside a bucket the chunk classes take less energy than typical, |a_ch|² ≈ |τ|²/4 + kσ²/2 against kσ², so the rest
  carries more. In the sphere model that costs 0.15 bits per rest class; 0.11 was measured.
- The tilt in favour of larger |τ|² is about 0.004 bits per unit of |τ|². Pair counts fall much faster, so the
  separable best-first score is the expected-count order.
- The Gaussian correction of `scripts/bucket_rank.py` is about 0 at 311 and at 223.
- The ranks of the two known cores' tuples (311, 223) are deeper than the model's shares suggest (§4), so the
  top-bucket rate may be up to about 3× optimistic.

**(c) Cost model** (`scripts/ladder_cost.py`, calibrated on (b)):
- log₂E_top = Σ l − 2pk − b_r (n_ind − k) − 0.08 (n_ind − k);
- buckets per core = m / E_top, where m is the number of tuples per core in its bucket (32 for closure pools, 2 or 1
  for reversal-canonical pools);
- bucket time ≈ 6·10^−8 s × (|A| + |C| + pairs).

## 3. Rung 1: v = 311, d = 5, profile (21, 21, 19, 1)

| quantity | value |
|---|---|
| ρ | 0.32 (Alpöge's blocks: 0.3245, 0.306, 0.3245, 0.3154; his two row-sum-21 blocks lie above the cut) |
| flat set [M] | 207,555 / 361,882 / 436,396 orbits (row sums 1 / 19 / 21) = 2^23.62 / 2^24.42 / 2^24.69 blocks; the union of 4 generation chunks exceeds the capture–recapture estimate, and Alpöge's two blocks with ρ ≤ 0.32 were already in the pools |
| n′ | 81.26 Gaussian; partial joins +0.03 per class at k = 8–10; sphere 84.2 |
| solution tuples in the closure pools | 2^13.2, i.e. about 10 distinct cores [E] |
| run A: closure pools, k = 10 | 842 buckets, 0 hits (cluster of 32: one core per about 4,700 buckets expected; consistent) |
| run B: reversal-canonical pools (2^23.69 / 2^23.69 / 2^23.42 / 2^22.62), k = 9 | 74M AB and 33M CD pairs per bucket, 6.5 s with pmatch (18.6 s before); model one core per about 1,200 buckets (about 2 h) |

**Result of run B [M]: a new 311 core**, `data/core_v311_b0_d5_0.txt`.
- Found at bucket 1,118 (τ = (6, −2, 2, 2, 2, 2, −2, 2, −2)) after about 2.3 h of k = 9 join. The model predicts one
  core per about 1,200 buckets. In total 1,118 k = 9 buckets plus 842 k = 10 closure buckets (worth about 0.2 cores)
  were searched.
- The bucket held 2 tuples, the core and its A ↔ B swap: the cluster multiplicity m = 2 of reversal-canonical pools.
- Verified by `scripts/verify_b0_core.py`, which shares no code with the search (`logs/k311/verify_core_v311_0.json`):
  - integer PAF sum 1244·δ;
  - invariant under u = 6 (order 5);
  - Goethals–Seidel H Hᵀ = 1244·I exactly.
- **Inequivalent to Alpöge's 311 core** under shifts, negation and reversal of single blocks, block permutations and a
  common multiplier: different PAF-multiset class key [P, M].
- Block ρ = 0.298, 0.314, 0.317, 0.257 (row sums 21, 21, 19, 1). It is deeper than Alpöge's (0.306–0.325) and well
  inside the flat set, as the existence count (about 10 cores at ρ ≤ 0.32) allows.
- Measured rate: one core in 1,118 + about 170 (k = 10 equivalent) buckets, against the model's about 1,200. The
  per-bucket model of §2 (sphere drift, in-bucket inflation, m = 2) is consistent with the first hit; with one event
  the calibration is good to about a factor of 3.

## 4. Rung 2: v = 223, d = 3, profile (23, 11, 11, 11)

| ρ | flat set (row sums 11 / 23) [M, lower estimate] | n′ Gaussian | n′ sphere | distinct cores in the whole flat set [E] |
|---|---|---|---|---|
| 0.30 | 2^25.9 / 2^23.7 | 84.47 | 88.0 | ~3 |
| 0.31 | 2^27.9 / 2^26.1 | 86.73 | 90.2 | — |
| 0.32 | 2^29.0 / 2^27.4 | 87.77 | 91.3 | ~2,500 |
| 0.34 | 2^32.4 / 2^30.7 | 91.56 | 95.1 | — |

- **Bucket order** (`logs/k223/rank_core_k*.jsonl`): Alpöge's 444 tuples (74 rotations × 6 slot assignments) rank
  142 / 9k / 86k (best / 10% / median) at k = 8. So "not among the first 352" (the earlier planted test) is a typical
  draw. 4 of the 444 tuples fall in the first 1,000 buckets, where the model predicts 6.3%.
- **Plan** (`logs/k223/run_prod223.sh`): ρ = 0.32, generation to 2^28 closures (rs 11), reversal-canonical slots
  (2^24.8 / 2^27 / 2^27 / 2^27, with B = C = D one shared pool) and the disk join at k = 12. Model: about 2.7 h per
  core after about 4 h of generation.

**Production (2026-10-17, from 22:35 UTC):**
- Generation stopped after 26 chunks: 2,814,949 / 663,818 orbits (row sums 11 / 23); three container reclaims cost
  about 1.5 h.
- Reversal-canonical slots: A = 2^24.55, B = C = D = 2^26.63 (one shared pool).
- Disk join at k = 12: about 46 s per bucket (77M AB × 273M CD pairs, 2 passes).
- Model: about 600 buckets (about 7.7 h) per core.
- **Result [M]: a new 223 core**, `data/core_v223_b0_d3_0.txt`.
  - Found at bucket 773 (τ = (2, 2, −2, −2, 2, 2, −2, 2, −2, 2, 2, 2)); the model predicts one core per about 600
    buckets. The bucket held 2 tuples: the core and its C ↔ D swap.
  - Verified by `scripts/verify_b0_core.py` (`logs/k223/verify_core_v223_0.json`): integer PAF sum 892·δ, invariant
    under u = 39 (order 3), Goethals–Seidel H Hᵀ = 892·I.
  - **Inequivalent to Alpöge's 223 core** (different PAF-multiset class key).
  - Block ρ = 0.309, 0.295, 0.301, 0.308 (row sums 23, 11, 11, 11); Alpöge's are 0.3215, 0.301, 0.301, 0.301.
  - Join time: about 6.5 h at 47 s per bucket. The 6 GB pass buffers had been squeezing the 9 GB of memory-mapped
    prep out of RAM. With `--mem 4e9` a bucket takes 15–18 s, so about 3 h per core.

## 5. Rung 3: v = 491, d = 7, profile (29, 27, 15, 13) — stopped [E]

| ρ | flat set per row sum [M] | n′ Gaussian | n′ sphere | distinct cores in the whole flat set |
|---|---|---|---|---|
| 0.32 | 2^27.7–2^27.9 | 103.48 | 106.8 | 2^−4.8 |
| 0.34 | 2^30.5–2^30.7 | 106.35 | 109.7 | ~11 |
| 0.36 | 2^33.0–2^33.3 | 109.01 | 112.3 | ~2,400 |

- Pools of 2^27 (about 4 GB per slot of prep on disk; the session has about 30 GB) hold 2^−1.7 solution tuples at
  ρ = 0.34.
- The whole flat set needs about 1,500 h and 2^46.4 operations per core.
- Alpöge's 491 blocks have ρ = 0.355, 0.338, 0.356, 0.311, just above the existence edge (about 0.33–0.34; the whole
  ρ ≤ 0.32 set holds about 2^−4.8 cores).

## 6. 859/13 (open): estimate and decision — no-go

Measurements: `logs/k859/` (capture–recapture at ρ = 0.34 and 0.36 with 65k–185k recaptures; at ρ = 0.38–0.40 the sets
are too large for recaptures, so their flat-set sizes are extrapolated [E]).

| ρ | flat set (row sums 1 / 25 / 53) | n′ Gaussian [M] | n′ sphere [E] | distinct cores in the whole flat set [E] |
|---|---|---|---|---|
| 0.34 | 2^28.7 / 2^28.75 / 2^27.55 [M] | 112.53 | 115.7 | 2^−12 |
| 0.36 | 2^31.1 / 2^31.3 / 2^30.2 [M] | 115.13 | 118.3 | 2^−4.6 |
| 0.38 | ≈ 2^33.4 / 2^33.6 / 2^32.7 [E] | 117.39 | 120.5 | ≈ 6 |
| 0.40 | ≈ 2^35.5 / 2^35.7 / 2^34.9 [E] | 119.38 | 122.5 | ≈ 550 |

- Multiplicity of a core among the tuples: 66 rotations × 2³ independent block reversals × 2 (C ↔ D).
- **Cores exist only from ρ ≈ 0.37.** Pools must hold most of a flat set of 2^33–2^35 blocks per slot; random subsets
  lose a factor (L/N)^4.
- **Cost of a plain 4-list k-tree** (`logs/ladder/cost_inputs_859.json`, `scripts/ladder_cost.py`), using the whole
  ρ ≤ 0.40 flat set as reversal-canonical pools and a key wider than 64 bits (k = 10, pairs about L per bucket):
  - E_top ≈ 2^−15.5 tuples, i.e. about 2^16.5 buckets per core;
  - about 2^53 element operations per core;
  - roughly 3 TB of prep data (4 slots × 2^34.5 blocks × about 30 B);
  - about 170,000 h on this machine.
  - At ρ = 0.38 the same costs about 2^55.
  - With pools of 2^28–2^30 (what this machine or one GPU could hold), the pools contain 2^−8.5 to 2^−2.5 solution
    tuples, so no core is expected at all.
- **Kaggle:** 30 GB RAM, 2 × T4 (16 GB), about 20 GB of output per notebook, 12 h sessions, about 30 GPU-h per week.
  That is three orders of magnitude short in storage and about 10^4 short in compute. No Kaggle plan is viable.
- **What it would need:** a cluster-scale k-tree over 2^34–2^35 flat blocks per slot (about 3 TB, distributed by key
  prefix), about 10^16 element operations per core (about 10 A100-days at 10^10 merged elements per second), and
  about 300 core-hours of flat-block generation at ρ ≈ 0.40 (multiplier closure included).
  - Or a better join: the cost exponent here is (1 + b_r/p) ≈ 2.1 in the pool size and b_r · n_ind ≈ 122 bits in
    the slice.
- 859/11 (n = 204.4), 739/9 (209.7) and 659/7 (236.6) are larger still. They are not pursued.

## 7. Rung 4: v = 283, d = 3, profile (19, 19, 19, 7) — measurement

- **n′ [M]:** 118.33 Gaussian at ρ ≤ 0.30, with +0.031 bits per class at k = 12 (`logs/k283/nprime_r0.30_m10.log`;
  the old value was 118.1). About 122.8 with the sphere excess (n_ind = 46).
- **Flat-set size, universal curve [E].** Across every slice measured here, log₂(flat set per row sum)/r is nearly a
  function of ρ alone:

  | ρ | 0.30 | 0.31 | 0.32 | 0.34 | 0.36 |
  |---|---|---|---|---|---|
  | log₂ N / r (mean, range) | 0.337 (0.316–0.349) | 0.360 | 0.382 (0.365–0.392) | 0.426 (0.410–0.433) | 0.464 (0.451–0.469) |

  - Data: 311, 223, 491, 859 by capture–recapture, 487 by Wang–Landau.
  - The low ends belong to unbalanced row sums (223's 23, 859's 53).
  - At 283 (r = 95), ρ ≤ 0.30 gives 2^32.0–2^33.2 blocks per row sum. The old rough "~2^33" agrees.
- **Capture–recapture at 283, ρ = 0.30 [M]** (`logs/k283/capture_r0.30_pooled.txt`): four independent samples.
  Two are 60 min × 1 core, with 5,801 / 5,933 orbits at row sum 19; two are 60 min × 4 cores, with 22,284 / 22,801.
  - Pooled pairwise Lincoln–Petersen (Σ nᵢnⱼ / Σ mᵢⱼ):
    - row sum 19: 3 recaptures, **about 2^35.0 blocks** (95%: 2^33.7–2^37.6);
    - row sum 7: 2 recaptures, **about 2^34.2** (2^32.7–2^36.2).
  - The estimate converges, though only to about ±2 bits with so few recaptures. Tabu walks catch some orbits more
    easily, which biases it low.
  - This is about 2.4 bits above the universal curve. The ρ ≤ 0.30 flat set at 283 is larger than the other slices'
    curve suggests, as the old "flat sets larger than the estimate" option of the gap analysis allowed.
  - Caution for repeats: `flat_gen_w.py` seeds workers with 1000k + ρ·1000 + 7919·seed, so equal seeds at different
    lengths give nested, not independent, samples.
- **Existence:** a core occupies 94 rotations × 2³ block reversals × 6 slot assignments ≈ 2^12.1 tuples. With the
  measured sizes, the whole ρ ≤ 0.30 flat set holds 2^(3·35.0 + 34.2 − 122.8 − 12.1) ≈ 2^4.3, about 20 cores (about
  0.04–2 on the universal curve). Alpöge's 283 blocks have ρ = 0.291, 0.299, 0.287, 0.290.
- **Plain k-tree cost** (`logs/ladder/cost_inputs_283b.json`, wide key, k = 15): over the whole flat set in
  reversal-canonical pools (2^34 / 2^34 / 2^34 / 2^33.2), about 2^16.6 buckets and **about 2^53 operations per
  core**. Pools of 2^30 hold 2^−2.8 tuples.
- **Conclusion [E]:** a plain flat-pool k-tree could have produced Alpöge's 283, but only at cluster scale. The search
  needs essentially the whole ρ ≤ 0.30 flat set (about 2^35 blocks per slot, tens of TB of prep data) and about 10^16
  element operations per core (about 1,000 GPU-hours at 10^10 per second). That is the same scale as 859/13 (§6).
  It is far beyond one machine, but plausible for a well-resourced group. A substantially better join than a plain
  4-list k-tree would lower it.
