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
- **Rung 1 (311): the method works at scale; the rate is slow [M/E].**
  - The pools are the whole reachable ρ ≤ 0.32 flat set.
  - A production-scale planted test recovers Alpöge's core.
  - The measured and modelled rate is one core per about 1,200 k = 9 buckets, about 5 hours of join per core on 4 cores.
  - Production result: see §3.
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

(Sections 3–7 below are filled in as the runs finish.)

## 1. What changed in the method

| change | effect | where |
|---|---|---|
| Full rotation closures in all four slots (no symmetry breaking) | each core present in all e rotations (31–37 distinct buckets) | `build_slot_pools*.py --closure` |
| Reversal-canonical slots | one block per reversal pair: pools ÷ 2, pairs ÷ 4, cores per bucket unchanged | `--revcanon` |
| 6-bit (k ≤ 10) and 5-bit (k ≤ 12) packed keys | keyed classes beyond 8 (k = 8 was OOM at 311: 2^31 pairs per bucket) | `kt_join_slots*.py`, `kt_join_disk.py` |
| Partner-key order independent of τ [P] | A and C sorted once, not per bucket: 18 s → 6 s per bucket at 311 | same |
| Disk-capable join | presorted memory-mapped pools, keys computed in the merge, shared files for equal pools | `kt_join_disk.py` |
| Skip histogram merge when one pass fits | one third less pair work | `kt_join_slots.py` (fast mode) |

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
| run B: reversal-canonical pools (2^23.69 / 2^23.69 / 2^23.42 / 2^22.62), k = 9 | 74M AB and 33M CD pairs per bucket, 15–20 s; model one core per about 1,200 buckets (about 5 h) |

Results of run B: (pending).

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

Results: (pending).

## 5. Rung 3: v = 491, d = 7, profile (29, 27, 15, 13) — stopped [E]

| ρ | flat set per row sum [M] | n′ Gaussian | n′ sphere | distinct cores in the whole flat set |
|---|---|---|---|---|
| 0.34 | 2^30.5–2^30.7 | 106.35 | 109.7 | ~11 |
| 0.36 | 2^33.0–2^33.3 | 109.01 | 112.3 | ~2,400 |

- Pools of 2^27 (about 4 GB per slot of prep on disk; the session has about 30 GB) hold 2^−1.7 solution tuples at
  ρ = 0.34.
- The whole flat set needs about 1,500 h and 2^46.4 operations per core.
- Alpöge's 491 blocks have ρ = 0.355, 0.338, 0.356, 0.311, just above the existence edge (about 0.34).

## 6. 859/13 (open): estimate and decision

(ρ = 0.38 measurement and the Kaggle plan pending.)

| ρ | flat set (row sums 1 / 25 / 53) [M] | n′ Gaussian | n′ sphere | distinct cores in the whole flat set [E] |
|---|---|---|---|---|
| 0.34 | 2^28.7 / 2^28.75 / 2^27.55 | 112.53 | 115.7 | 2^−12 |
| 0.36 | 2^31.1 / 2^31.3 / 2^30.2 | 115.13 | 118.3 | 2^−4.6 |

## 7. Rung 4: v = 283 (measurement only)

(pending: ρ = 0.30 capture–recapture with long generations.)
