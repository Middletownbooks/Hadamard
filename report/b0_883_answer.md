# H(3532), v = 883: a b = 0 core in the d = 21 multiplier slice, found by the flat-pool k-tree

Date: 2026-10-16. Branch `claude/focused-planck-j1645f`.
Scripts: `kaggle/kt_join_kaggle.py` (search), `scripts/verify_h3532.py` (verification), `scripts/k883_pjoin.py`
(partial joins). Logs: `logs/k883/` (run history in `logs/k883/README.md`). Bundle: `pub/H3532/`.

**Tags.** [M] exact computation (script named); [P] proved here; [E] Gaussian estimate (estimator and calibration
named); [L] literature recollection, not checked from this container.

---

## 0. Bottom line

1. **H(3532) is constructed and verified [M].** v = 883 was on the open list (`report/prompt_alpoge_design.md`), and
   3532 is "Unknown" in `data/sage_unknown_orders_le4000.json` (Sage `hadamard_matrix(n, existence=True)`).
   - **Core.** Four ±1 circulants A, B, C, D of length 883 with Σ_k PAF_k = 3532·δ and row sums (43, 41, 1, 1).
     File: `data/core_v883_b0_d21_0.txt`.
   - **Matrix.** The Goethals–Seidel array of the four circulants (`g` layout, back-identity R).
   - **Verification.**
     - `scripts/verify_h3532.py`: integer PAF, U-invariance of all four blocks, and H Hᵀ = 3532·I in int64
       (`logs/k883/verify_h3532_core0.json`).
     - The standalone `pub/H3532/verify.py` checks the sequences and the published normalized matrix
       (`logs/k883/pub_verify_h3532.txt`).
     - SHA-256 of the packed rows (`data/H3532_packed.npz`):
       `5d131c4fc91ad40c2e7b98a37554825cb2c65c0ba67e5cb94231841d4546790c`.
   - **The 32 tuples found are one core [M].** Each of the 32 lines of `data/k883_solutions_v883_run4.jsonl` passes the
     integer PAF and U-invariance checks (`logs/k883/verify_all32_paf.json`). Up to reversing any block (x ↦ −x) and
     swapping C and D, all 32 are the same core: 32 = 2⁴ · 2.
2. **The slice [P]** (§1).
   - U = ⟨729⟩ is the subgroup of order 21 of F₈₈₃^* (729 = 2⁴², and 2 is a primitive root mod 883).
   - All four blocks are U-invariant and untwisted.
   - Each block has 43 orbit values: {0} and 42 U-orbits. There are 21 lag classes (⟨U, −1⟩-orbits, 42 lags each).
   - Estimator: n = 111.9 bits, log₂E = 60.1 [E] (`report/b0_919_answer.md` §4).
   - Restricted to flat blocks (PAF ratio ≤ 0.45): n′ ≈ 78 [M/E] (`logs/k883/README.md`).
3. **Search [M]** (§2).
   - Method: the flat-pool k-tree that found the new 487 cores (`report/REPORT.md` §4–5, `report/v487_new_solutions.md`).
   - Flat pools: 2^24.3 to 2^24.6 blocks per row sum. Join: 7 chunk classes, then best-first τ-buckets.
   - Result: bucket 232 of run 4 returned 32 tuples, all verified.
   - Time: 79 minutes on 4 cores for run 4, reusing run 3's pools. All four runs together took about 5.6 hours of
     wall-clock time on the 4-core container.
4. **Calibration [M/E]** (§3).
   - The pipeline logs its own expectation, using 3.42 bits per remaining class (calibrated at 487).
   - At 883 the measured value is 3.97 bits per class. With that, the realistic expectation before the hit is about
     25 tuples.
   - Solutions arrive in clusters of up to 32 tuples, so 25 tuples is about one core.
   - One core was found.

## 1. The slice [P]

- **Multiplier.** U = ⟨u⟩ with u = 2⁴² ≡ 729 (mod 883), the subgroup of order 21 of F₈₈₃^* (882 = 2·3²·7²).
  Every block satisfies X(ux) = X(x).
- **Equations.** −1 ∉ U because 21 is odd. So the lag classes are the 882/42 = 21 orbits of ⟨U, −1⟩.
  ΣPAF = 4v·δ becomes 21 class equations Σ_X PAF_X(c) = 0.
  - One equation is implied by the others: Σ_c 42·PAF_X(c) = s_X² − v for each X, and Σ s_X² = 4v.
- **Row sums.** A U-invariant block has s = X(0) + 21k, with k the sum of 42 orbit signs (so k is even). Hence
  s ≡ ±1 (mod 42), and |s| ∈ {1, 41, 43} because 83² > 3532.
  - The only way to write 3532 as four such squares is 43² + 41² + 1² + 1².
  - So the profile is (43, 41, 1, 1) up to signs and order. It is the unique admissible profile.
- **Freedom.** 4 × 43 = 172 signs against 20 independent class equations. Gaussian n = 111.9 bits
  (`scripts/b0_open_screen.py`, untwisted). This equals the T-form cost of `report/tmat_routes.md` §10.

## 2. Search [M] (`kaggle/kt_join_kaggle.py`)

1. **Flat pools.** Tabu walks generate U-invariant blocks with PAF ratio sqrt(F/(mL)) ≤ ρ for each row sum, where F is
   the sum of squared class PAFs. ρ = 0.45 in the successful run.
   - Run 3 generated for 16 minutes on 4 cores. After deduplication the pools held 2^24.26 blocks (row sum 43),
     2^24.36 (41) and 2^24.61 (1).
   - C and D share the row-sum-1 pool.
2. **Coverage of the flat set [M/E].**
   - Reversal x ↦ −x maps U-orbits to U-orbits and keeps every PAF, so the flat set is reversal-closed.
   - The fraction of pool members whose reversal is also in the pool is 0.885 / 0.884 / 0.908 (row sums 43 / 41 / 1;
     200,000 samples each; `logs/k883/diag/reversal_closure.log`).
   - If a pool is a uniform sample of its flat set, that fraction is the coverage. So the pools hold about 90% of all
     blocks with ρ ≤ 0.45, and generating longer would not have grown them much (`logs/k883/README.md`, run 3).
3. **Join.**
   - Classes are taken in a fixed random order: 7 chunk classes, one histogram class t, and 13 rest classes hashed
     linearly to 64 bits. The 13 include the one implied class.
   - Buckets are τ-vectors on the chunk classes, most probable first. Bucket τ holds the AB pairs whose PAF sums equal
     τ on the chunk classes, and the CD pairs at −τ.
   - Class t is matched through histograms. The rest are matched by exact equality of the 64-bit keys.
   - Every candidate is checked on all 883 lags before it is written.
   - Bucket 232 (τ = (2, 2, 2, 6, 6, 2, −2)) held 2^23.8 AB pairs and 2^24.4 CD pairs, and gave 32 candidates, all
     solutions.
4. **Run 4 configuration** (`logs/k883/README.md`): KT_V=883 KT_D=21 KT_PROFILE=43,41,1,1 KT_RHO=0.45 KT_K=7
   KT_LMAX=41300000 KT_MEM_BUDGET_BYTES=3e9 KT_INSTANCE=2 KT_SEED=0. Pools come from run 3, which used the same
   settings with KT_K=9 and KT_GEN_MINUTES=16.
   - Preparation took 11.5 minutes; the join ran at 17.4 s per bucket.
   - The hit came at t = 4,761 s. The worker was stopped at bucket 252.

## 3. Calibration [M/E]

- **Per-class cost.** Partial joins on the production pools (`scripts/k883_pjoin.py`, 8,192 blocks per row sum, 2^52
  tuples; `logs/k883/diag/pjoin_883_8192.log`):

  | classes k | 1 | 4 | 8 | 12 | 13 |
  |---|---|---|---|---|---|
  | measured bits | 3.97 | 15.95 | 31.99 | 48.0 (16 matches) | 51.0 (2 matches) |

  That is 3.92–4.0 bits per class with no drift through 13 classes. The earlier 4,096-block measurements agree
  (`logs/k883/nprime_rho0.45*.log`).
- **Expected tuples against the outcome.**
  - The pipeline's estimate uses 3.42 bits for each remaining class. Replacing it by the measured 3.97 bits
    multiplies the estimate by 2^{−0.55·r}, where r is the number of remaining classes (11 for k = 9, 12 for k = 7).

  | run | buckets | pipeline's own estimate (tuples) | at 3.97 bits per class | found |
  |---|---|---|---|---|
  | 1 (ρ = 0.39, pools 2^21) | 742 | 0.45 | 0.03 (3.77 bits per class at ρ = 0.39) | 0 |
  | 2 (ρ = 0.45, k = 9) | 345 | 139 | 2.1 | 0 |
  | 3 (ρ = 0.45, k = 9) | 248 | 161 | 2.4 | 0 |
  | 4 (ρ = 0.45, k = 7), to bucket 232 | 233 | 2,011 | 20.7 | 32 (one core) |

  - Before the hit, the expectation adds up to about 25 tuples. Adding the +0.17 bits per class seen at 487 would give
    about 7.
- **Clustering.** Tuples are not independent events.
  - A core (A, B, C, D) gives up to 2⁴ tuples from reversing any of its blocks, all in the same bucket with identical
    class vectors. Swapping C and D doubles that, because they share a pool.
  - With about 90% coverage per block, a core present in the pools shows about 32·0.9⁴ ≈ 21 tuples on average.
  - So about 25 expected tuples is about 1.2 expected cores, and finding one core in bucket 232 is in line with that.
  - The zeros of runs 2 and 3 (about 4.5 tuples, i.e. about 0.2 cores) were not evidence against the model.
  - Lesson for the estimator: count expected clusters, not tuples, whenever the pools are closed under a symmetry
    that fixes the class vector.
- **Per core [E].** At the run-4 rate (about 0.089 tuples per bucket at 17.4 s per bucket), one core costs about 240
  buckets, roughly 70 minutes on 4 cores once the pools exist.

## 4. Notes

- **Inequivalent cores.** A second core would need another few hundred buckets. The run stopped after the first core.
  `kaggle/kt_join_kaggle.py` resumes from its cache and log; the pools are in the session scratchpad and are not
  committed, being about 1.6 GB.
- **Multiplier equivalence.** w ∈ F₈₈₃^* \ ±U maps the core to another U-invariant core with its classes permuted. Such
  images land in other buckets when all four images are in the pools. They are equivalent and were not searched for.
- **Literature [L].**
  - The repository's open list and Sage both have H(3532) unknown.
  - The construction is a plain multiplier-invariant Goethals–Seidel core, like the 487 cores of
    `report/v487_new_solutions.md`.
  - Whether H(3532) appears in work after Sage's tables cannot be checked from this container.

## Reproduce

```
python3 scripts/verify_h3532.py data/core_v883_b0_d21_0.txt          # exact check, writes data/H3532_packed.npz
cd pub/H3532 && python3 verify.py                                    # standalone bundle check
# search (pools are generated afresh; time to a core is random, about 1-2 h of join after ~40 min of pools + prep):
KT_V=883 KT_D=21 KT_PROFILE=43,41,1,1 KT_RHO=0.45 KT_K=7 KT_LMAX=41300000 KT_GEN_MINUTES=16 \
  KT_MEM_BUDGET_BYTES=3e9 KT_MAX_HOURS=60 KT_NBUCKETS=100000 KT_INSTANCE=2 KT_NINST=1 KT_SEED=0 \
  KT_OUT=OUTDIR KT_GEN_DIR=GENDIR python3 kaggle/kt_join_kaggle.py
```
