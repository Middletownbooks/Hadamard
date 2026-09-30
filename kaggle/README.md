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


## Fix (bucket assignment) — read this if an instance ended right after `"stage": "ready"`

Earlier versions gave bucket `b` to an instance only when `b % NINST == INSTANCE`. With `NINST=1` and `INSTANCE >= 1`
that matched **no** bucket, so the run finished generation and prep and then ended immediately (`"stage": "end"`,
0 buckets). Fixed: an instance now takes `b % NINST == INSTANCE % NINST`.

- **Independent parallel runs (recommended):** `NINST=1` and a different `INSTANCE` per notebook. Each instance
  builds its own pools, seeded by INSTANCE, and works through all buckets best-first on those pools.
- **Split one pool set:** `NINST=N` and `INSTANCE=0..N-1`. This only makes sense if the instances share pools,
  which they do not by default.

**Recovering a run that ended at "ready":** its pools and prep are saved in `OUT/cache_inst{INSTANCE}`.
1. Add that committed output as an input dataset.
2. Set `RESUME_FROM` to its `kt_*` directory and keep the **same INSTANCE**.
3. Run the fixed notebook. The cache is read in place, so it skips the multi-hour generation and prep and starts on
   buckets within minutes.

## v=515 exact two-sided lift, compressions 5, 12 and 13 (`kt515_u161_sol5_lift.py`, `kt515_u161_sol12_lift.py`, `kt515_u161_sol13_lift.py`)

A self-contained program for one notebook cell (numpy + numba, no internet, no input files). It embeds the solver
modules and data. It runs the exact lift on the 56 yc classes of u=161 compression 5; the classes of compressions 3
and 4 run locally.

1. **Paste the whole file into one cell.** Edit only `INSTANCE` and `NINST` at the top: run `NINST` notebooks with
   `INSTANCE = 0 … NINST−1`, and the same `NINST` in each. Notebook i takes classes c with c % NINST == i. A bad
   INSTANCE/NINST pair aborts immediately with `ASSIGNMENT CHECK FAILED`.
2. **To pick classes by hand,** set `CLASS_LIST = [...]`. This overrides INSTANCE/NINST.
3. **Self-test.** With `SELFTEST = True` (the default), the first minute checks the machinery on the known v=347
   matrix. It must print `SELFTEST PASS`.
4. **What to watch in the cell output:**
   - `ASSIGNMENT CHECK` shows your class list.
   - For each class: `class c: lift list sizes …`, then `CLASS RESULT {...}`.
   - `SOLUTION {...}` lines carry the sequences and `hadamard_verified`, checked on the full 2060 × 2060 matrix.
5. **Timing and memory.** The estimate is roughly 3–6 hours per class on 4 cores [E, not yet measured on Kaggle].
   After `MAX_HOURS = 11` it stops *starting* new classes; a class cut off mid-join is simply redone next time. Class 0
   stores 5.0·10⁸ pairs, one batch of about 8 GB at `JOIN_BATCH = 1e9`. Lower `JOIN_BATCH` if memory is short: this
   costs time, not correctness.
6. **Continuing.** Save the version, add its output as an input dataset, set `RESUME_FROM` to the `lift515_c<N>` folder
   inside it, and keep the same INSTANCE/NINST. Finished classes are skipped. The output is tiny (a log line per class).
- `kt515_u161_sol12_lift.py` is the same program for compression 12, which has only **4 classes**. Use 4 notebooks,
  `NINST = 4` and `INSTANCE = 0..3`, so each runs one class (about 3–6 h [E]). Its output folder is `lift515_c12`.
- `kt515_u161_sol5_lift.py`: 56 classes, output folder `lift515_c5`.
- `kt515_u161_sol13_lift.py`: 36 classes, output folder `lift515_c13` (heavy lists, about 3–6 h per class [E]).

## v=523 profile A with the INSTANCE fix (`kt523_A2_35-17-17-17.py`)

This is the same program as `kt523_B2_37-19-19-1.py`, with profile A (35,17,17,17), 90 min of generation, and output
in `/kaggle/working/kt_A_35-17-17-17`. Edit only `INSTANCE` (any value 0, 1, 2, …; each INSTANCE gets its own
independent pools). Caches go to `/tmp` by default, and `SOLUTION` lines print in the cell output.

Profile A generates and preps two pools (row sums 35 and 17) instead of three, so buckets start about 45–60 min
earlier than with B. Tested locally at INSTANCE=2 (small pools): the assignment check passed, both pools were built,
and 398 buckets ran with non-zero A+B and C+D pair counts.

## Williamson matrices of order 65 (`wil65.py`)

Exact search, pre-registered in `williamson/PREREGISTRATION.md`. The total is about 8.2×10^13 pair operations
(roughly 300–500 core-hours), split into 5 row-sum cases × K = 32 shards = 160 shards.

1. Paste `wil65.py` into a CPU notebook (no GPU or internet needed). Edit the settings block:
   - `CASE` (0–4);
   - `K = 32` (the same for every session);
   - `SHARDS`: suggested 8 per session, e.g. `[0..7]`, `[8..15]`, `[16..23]`, `[24..31]`. That is 4 sessions
     per case and 20 sessions in total.
   - `NPROC = 4`.
2. Run All. It first runs the n = 45 self-test (about 2–5 min); the session aborts if the self-test fails.
3. Each finished shard prints `shard i: done ledger=... work=... sols=...`. Any hit prints a line starting with
   `SOLUTION`.
4. Download `/kaggle/working/wil65/wil_n65_c<CASE>_s<i>of32.json` for every shard (the `.log` files are
   optional). If a session dies, rerun only the missing shards: finished shard files are skipped.
5. Put all the JSONs in one folder and run `PYTHONPATH=scripts python3 scripts/wil_collect65.py <folder> 32`. It
   checks that every shard is present and that the ledger and work equal the exact totals in
   `logs/wil_n65_cost.json`, and it re-verifies any solution.

Shards differ in size: each takes every 32nd Z5 representative. Expect roughly 2–4 hours per shard per core.
