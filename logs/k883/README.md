# v = 883, b = 0, d = 21, profile (43, 41, 1, 1): flat-pool k-tree run

- Slice: `report/b0_919_answer.md` §4 (n = 111.9, log₂E = 60.1; 43 orbits, 21 classes). The profile (43, 41, 1, 1) is
  the only admissible one (block sums ≡ ±1 mod 42, Σ s² = 3532).
- Pools (`scripts/s1_pools.py 883 21 0.39 180 1 7`, one core, 218 s): 1.27M / 1.12M / 1.10M distinct blocks with
  PAF ratio ≤ 0.39 at row sums 1 / 41 / 43 (`pools_rho0.39_1core_218s.log`).
- Flat-subspace constraint bits (`scripts/s1_nprime.py`):
  - Gaussian n′ = 73.0 (pools of 50,000 per row sum), 3.65 bits per class;
  - partial joins (pools of 4,096): 3.77 / 7.54 / 11.34 bits measured against 3.65 / 7.31 / 10.97 Gaussian for
    k = 1 / 2 / 3, i.e. +0.12 bits per class, so n′ ≈ 75 [M/E].
  - Compare v = 487 (d = 9): n′ ≈ 83, where the same pipeline found two new cores (`report/v487_new_solutions.md`).
- Production: `kaggle/kt_join_kaggle.py` with KT_V=883 KT_D=21 KT_PROFILE=43,41,1,1 KT_RHO=0.39 KT_K=9
  KT_LMAX=33554432 KT_GEN_MINUTES=40 KT_MEM_BUDGET_BYTES=3e9 (4 cores). Logs: `log_inst0*.jsonl` here.

## Run 1 (rho = 0.39): pools saturated, stopped

- The flat set at rho <= 0.39 is small: 40 min on 4 cores gave only 2.32 / 2.32 / 2.43 M distinct blocks (row sums
  41 / 43 / 1), and a single 5-minute generator run already reached 1.8–2.0 M (`gen_yield_5min_rho045_rho039.log`).
- 742 buckets ran (about 1 s each, AB ≈ 42k and CD ≈ 46k pairs per bucket, 0 candidates); log
  `log_inst0_rho039_L2e21.jsonl`. With L ≈ 2^21.15 the pools hold about L⁴/2^{n′} ≈ 2^9 solutions, too few.

## Run 2 (rho = 0.45, pools 2^25)

- Yield at rho = 0.45: 8.6–11.9 M distinct blocks per row sum in 5 minutes on one core, not yet saturating.
- n′(0.45) = 77.4 Gaussian (pools of 50,000). Partial joins give 3.97 / 7.92 / 11.89 bits measured against
  3.95 / 7.89 / 11.83, so n′ ≈ 77.8 (`nprime_rho0.45.log`).
- With L = 2^25 the pools hold about 2^{100 − 77.8} ≈ 2^22 solutions, against 2^9 in run 1.
- Production: the same pipeline with KT_RHO=0.45, KT_INSTANCE=1, KT_LMAX=33554432, KT_GEN_MINUTES=8. Generation
  is shortened from 40 minutes because at rho = 0.45 the raw output is about 6x larger, and pool assembly
  (concatenate + unique) would otherwise peak near the memory limit. Eight minutes on 4 cores gives about
  2^25 distinct blocks per row sum.

**Run 2 result:** 345 buckets at 12.8 s each, 0 candidates (`log_inst1_rho045.jsonl`). Pools were 2^24.05 to 2^24.44.
- The pipeline's own estimate assumes 3.42 bits per class (calibrated at 487). At 883 the measured value is about
  3.97, so the realistic expectation is about 2^{36.4 − 43.7} ≈ 0.006 per bucket, about 2 over the run.
- **Same-shape control [M]:** v = 127, d = 3 has the same 43 orbits and 21 classes (n = 80.7, so it is much easier).
  The same pipeline found 1.46 M verified cores in its first bucket. The mechanics are correct; rate calibration at
  883's difficulty needs the 883 run itself.

## Run 3 (rho = 0.45, pools 2^25.3)

- Pool size enters hits per bucket as L⁴, while bucket time grows about as L, so the pools were enlarged:
  KT_LMAX=41300000, KT_GEN_MINUTES=16, KT_INSTANCE=2.
- Expected: about 16–20× more hits per bucket for about 2× the bucket time.
