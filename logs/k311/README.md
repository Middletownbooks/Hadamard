# v = 311, b = 0, d = 5: flat-pool k-tree, rung 1 of the calibration ladder

Purpose: stretch the flat-pool k-tree (H(3532) at 883, two new cores at 487) to a published order whose slice is
larger than anything we have solved (n = 139.4 untwisted constraint bits, against 130.9 at 487 and 111.9 at 883),
and measure how far the method reaches (`report/prompt_flatpool_stage1.md`, `report/alpoge_gap_analysis.md` §2).
Alpöge's 311 blocks sit at mean PAF²/v 0.094–0.105, i.e. flatness ratio rho ≈ 0.31–0.32.

- Slice: U = order-5 subgroup of F₃₁₁^*, 63 orbits (r = 63: fits the 64-bit masks), 31 lag classes, e = 62.
- Profiles (row sums ≡ ±1 mod 10, Σ s² = 1244): (21, 21, 19, 1) (Alpöge's), (21, 19, 19, 9), (31, 11, 9, 9).
  Runs use (21, 21, 19, 1).

## New pool machinery

- **Multiplier closure [P, M].** x ↦ x(g^s y) rotates the e = 62 orbit bits, keeps the row sum and permutes the lag
  classes (all of size 2d), so it keeps F = Σ_c PAF(c)²: flat sets are unions of rotation orbits. Checked
  numerically (`scripts/flat_probe.py`). Uses: (i) coverage of a pool = fraction of rotations already present;
  (ii) every generated block gives up to 62 flat blocks for free.
- **Capture–recapture** of the flat-set size from two independently seeded generations, at the level of multiplier
  orbits (`scripts/flat_capture.py`; a lower estimate, since tabu walks reach some orbits more easily).
- **Slot pools** (`scripts/build_slot_pools.py`) and a **slot join** (`scripts/kt_join_slots.py`, the kaggle join
  with one pool per slot). Symmetry breaking (slot A = one block per orbit) is supported but only pays when samples
  cover most of each closure; with low coverage, equal closure samples in all slots give more tuples per unit of work.

## Probes (1 min on 4 cores each, `probe_rho*_s*.log`)

| rho | blocks per row sum per minute | closure coverage | Gaussian n′ | partial joins (k ≤ 10) | flat set (capture–recapture) |
|---|---|---|---|---|---|
| 0.40 | 3.8 M | 1.6% | 89.9 (3.0 per class) | +0.035 bits per class | far larger than any sample |
| 0.36 | 0.6 M | 1.7% | 86.0 (2.87 per class) | +0.031 | 2^28.8–2^29.0 blocks |
| 0.32 | 30–46 k | 1.7% | 81.0 (2.70 per class) | +0.037 | 2^23.5 / 2^24.3 / 2^24.6 (row sums 1 / 19 / 21) |

Choice: rho = 0.32. Its flat set (about 2^24 per row sum) fills pools of 2^24 per slot almost completely after a
few 10-minute chunks (each chunk finds 60–70% of the orbits); solution tuples in the pools ≈ 2^(95.5 − 82) ≈ 2^13
(each core in several rotated copies, landing in different buckets). rho = 0.36 with pools twice as large gives a
similar count with slower buckets and a less concentrated bucket distribution.

## Production (rho = 0.32)

- Generation: 4 chunks × 10 min × 4 cores (`gen_r0.32_c*.log`; chunks 5–6 were dropped once saturation was clear).
- Slot pools: all four slots uniform samples of 2^24 from the rotation closures (`slots_r0.32.log`).
- Join: `scripts/kt_join_slots.py 311 5 21,21,19,1 SLOTS OUT --k 8 --mem 3e9 --bits 2.75` (8 chunk classes keyed,
  1 histogram class, 21 classes hashed; per-bucket expectation uses 2.75 bits per remaining class, the partial-join
  value). Log: `join_r0.32_log.jsonl` (copied from the scratch run directory), solutions `solutions_v311.txt`.
- Runner: `runners/run_prod.sh RHO NCHUNKS L LA K BITS` (restart-safe: chunks, slots and buckets resume).
