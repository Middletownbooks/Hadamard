# Kaggle: flat-pool k-tree join (v=523, d=9 by default)

`kt_join_kaggle.py` is self-contained (numpy + numba; both preinstalled on Kaggle, no internet needed).
Paste it into one notebook cell, or upload it and run `!python kt_join_kaggle.py`.

## What one session does
1. **Pool generation** (`gen_minutes`, all cores): tabu walks emit flat blocks (PAF ratio ≤ `rho`) at the profile's row sums.
2. **Prep** (~35–50 min at 2^26 per row sum): per-block PAF on the lag classes → chunk keys, bucket class, 64-bit hash. Cached.
3. **Join**: buckets τ in most-probable-first order; each bucket = Wagner stage 1 (8-class key), t-masked passes under
   `mem_budget_bytes`, parallel Bloom matching, exact verification. Stops at `max_hours`.

Outputs in `/kaggle/working/kt/`: `solutions_v523.txt` (JSON lines, every entry verified Σ PAF = 0),
`log_inst{I}.jsonl` (per bucket: stage-1 sizes, candidates, `exp_hits` = model estimate), pools/prep caches.
Check solutions independently with `verify_file('/kaggle/working/kt/solutions_v523.txt')` (builds H and tests H·Hᵀ = 4v·I).

## Several instances in parallel
Copy the notebook N times; in each set `INSTANCE=0..N-1` and `NINST=N` in `CONFIG`.
Each instance generates its **own** pools (seeded by INSTANCE) and takes buckets b ≡ INSTANCE (mod N), so the runs are independent.

## Continuing across sessions
Save the session output ("Save Version" → output kept), then in the next session add that output as an input dataset and set
`RESUME_FROM='/kaggle/input/<dataset>/kt'`. Pools and prep are reused (no regeneration); finished buckets are skipped.
(Keep the same INSTANCE for a given resume chain.)

## Settings for Kaggle's larger RAM (≈30 GB, 4 CPU)
Recommended: `Lmax=2**27`, `mem_budget_bytes=10e9`, `gen_minutes=90`. Larger pools raise hits per bucket roughly as L⁴ while bucket time grows ~L,
so this is the main lever. Watch the first bucket's `AB`/`CD` and `passes` in the log.

## Expectations (estimates, not measurements)
At v=487 (validated): ~0.034 expected hits per bucket at 2^26 pools, 2.8 min/bucket on 4 cores; 2 new solutions in ~40 buckets.
At v=523 the flat-subspace constraint is ~6–7 bits larger, so expect ~50–100× fewer hits per bucket at the same pool size:
order 10⁻³–10⁻⁴ per bucket, i.e. many instance-sessions per solution. Replace this estimate with the logged `exp_hits` after the first buckets.

## Ready-made v=523 files (Lmax = 2^27, 30 GB Kaggle CPU sessions)
| file | profile | pools (row sums) | gen_minutes | output dir |
|---|---|---|---|---|
| `kt523_A_35-17-17-17.py` | (35,17,17,17) | 35, 17 | 90 | `/kaggle/working/kt_A_35-17-17-17` |
| `kt523_B_37-19-19-1.py` | (37,19,19,1) | 37, 19, 1 | 120 | `/kaggle/working/kt_B_37-19-19-1` |

These are the two admissible d=9 profile classes for v=523; they are independent shots (no prior favours either).
Each file is a complete paste-into-one-cell program; only edit `INSTANCE`/`NINST` (and `RESUME_FROM` when continuing).
Raw generator chunks go to `/tmp/kt_gen` and are deleted once pools are built, so `/kaggle/working` holds only the caches
(≈3.5 GB per row sum at 2^27: A ≈ 7 GB, B ≈ 10.5 GB) plus logs and solutions — under Kaggle's 20 GB output limit.
Expect per bucket ≈4× the stage-1 work of the 2^26 local run (≈10–15 min/bucket on 4 cores) for ≈16× the hit rate [estimate];
the local 2^26 run on profile A measured ~3.4 min/bucket and an estimated ~3×10⁻⁴ hits/bucket (unvalidated per-class bits at v=523).
