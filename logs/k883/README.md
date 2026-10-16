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
