# Prompt (Priority 1, high effort): scale the b = 0 flat-pool k-tree up the calibration ladder, then aim at an open order

*Self-contained, for a fresh Claude Code session on `middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`.*

- **Read first, in this order:**
  1. this file;
  2. `report/alpoge_gap_analysis.md` §2 (why);
  3. `report/prompt_flatpool_stage1.md` (the original task list, with its progress log);
  4. `logs/k311/README.md` (where rung 1 stopped);
  5. `logs/k883/README.md` (the run that built H(3532): what worked, how pools saturate, how the per-bucket rate was
     calibrated).
- **A separate max-effort session works on the b ≥ 3 stage-1 item** (`report/prompt_p2_bge3_stage1.md`). Do not edit
  its files (`theory/ab_s1_z*`, `logs/s1z/`).

## Rules

- **Pushing:** other sessions push to this branch. Before every push run
  `git fetch -q origin claude/focused-planck-j1645f && git merge -q --no-edit origin/claude/focused-planck-j1645f`,
  then `git push -u origin HEAD:claude/focused-planck-j1645f`.
- **Commits:** trailers follow the existing commits (`Co-Authored-By:` and `Claude-Session:`); no pull requests.
- **Processes:** never kill processes you did not start. Never `pkill -f` with a pattern that matches your own
  shell; kill by PID.
- **Restart safety:** write runners that resume from disk. Containers restart and usage ceilings interrupt sessions.
- **Check-ins:** whenever workers run, arm a check-in (`send_later`, about 1 hour) so the session never sits idle; an
  idle session was reclaimed once and killed production workers.
- **Commit as you go:** commit logs and small JSON as results arrive. Pools and caches stay out of git: keep them in
  the scratch directory and check `.gitignore`.
- **Verification:** a construction counts only if verified exactly (integer PAF, H Hᵀ = 4v·I) by a script sharing no
  code with the search.
- **Tags:** separate exact results [M], proofs [P] and estimates [E].
- **Progress log:** append to the progress log at the end of `report/prompt_flatpool_stage1.md` as you work, so the
  next session can pick up.
- **Compute:** 4 cores and about 15 GB of RAM. Install `numba` with pip if missing.

## Background (what you need to know)

**The array.** b = 0 Goethals–Seidel: four ±1 circulants of length v with Σ_k PAF_k = 4v·δ, rows Σ s_k² = 4v.
- Slice: every block is invariant under the order-d subgroup U of F_v^* (v prime). That gives 1 + (v−1)/d orbits,
  and lag classes = ⟨U, −1⟩-orbits.
- Constraint bits n per slice: `scripts/b0_open_screen.py`.

**Method (ours, and by the forensics also Alpöge's).** Build pools of *flat* blocks, then join four lists.
- A flat block is one whose PAF energy satisfies F ≤ ρ²·m·L; the blocks are local minima from per-block tabu
  descent.
- Restricting to flat blocks lowers the constraint count from n to n′. The measured gain is about 1.7–1.9 bits per
  lag class.
- The join is a 4-list k-tree (Wagner): k chunk classes keyed, best-first τ-buckets, a histogram class, a 64-bit hash
  of the rest, a bloom match, then exact verification.
- A 4-list join needs pools of about 2^(n′/3) blocks per row sum.

**Results so far:**
- **H(3532), v = 883, d = 21:** n = 111.9, n′ ≈ 78, pools 2^24.3–2^24.6, 4 runs, about 5.6 hours on 4 cores.
  Per-bucket calibration: 3.42 bits per remaining class was calibrated at 487; 3.97 was measured at 883.
- **Two new cores at 487** (d = 9; n = 130.9, n′ ≈ 83).
- **Failed:** 523 runs at n′ ≈ 91 (now moot, since H(2092) was built differently).

**Alpöge's b = 0 cores are deeper and his slices larger.** Mean PAF²/v is 0.082–0.105 at 223, 283 and 311 (ρ ≈ 0.29–0.32), against our best of 0.11 or more.

**The ladder** (published orders as rungs; each rung measures how far the method reaches):

| rung | v, d | orbits / classes | n | flat n′ [E] | lists 2^(n′/3) | state |
|---|---|---|---|---|---|---|
| 1 | 311, 5 | 63 / 31 | 139.4 | ≈ 81–84 | 2^27–28 | pools generated (ρ = 0.32, 3 chunks), **join not run** |
| 2 | 223, 3 | 75 / 37 | 156.6 | ≈ 90 | 2^30 | two-word join validated by a planted test |
| 3 | 491, 7 | 71 / 35 | 169.1 | ≈ 106 | 2^35 | needs disk-backed pools |
| 4 | 283, 3 | 95 / 47 | 207.0 | ≈ 119 (measured) | 2^40 | flat-set size at ρ ≤ 0.30 unknown (±5 bits) |

**The open-order target the ladder gates:** **859, d = 13, untwisted** (67 orbits, 33 classes, n = 173.4,
row sums (−25, 1, 1, 53); `logs/b919/b0_open_screen.jsonl`).
- With 1.7–1.9 bits per class its flat n′ ≈ 110–117 [E]: between rungs 3 and 4.
- Next in line: 859/11 (n = 204.4), 739/9 (209.7), 659/7 (236.6).
- 991/15 (n = 176.4) is 3-cycle-twisted, so its lists are not independent; it needs the 2-list approach of 919 and is
  not part of this prompt.

## State of rung 1 (311) when the last session stopped

- **Profile:** (21, 21, 19, 1) (Alpöge's). Others: (21, 19, 19, 9), (31, 11, 9, 9).
- **Flat-set depth ρ = 0.32 was chosen.** The flat set saturates at about 2^23.5 / 2^24.3 / 2^24.6 blocks per row sum
  1 / 19 / 21 (capture–recapture, lower estimates; `logs/k311/capture_rho0.32_c1c2.txt`).
- **Expected solution tuples in pools ≈ 2^(95.5 − 82) ≈ 2^13.** Each core lands in several rotated copies, in
  different buckets.
- **Generation chunks 1–3 logged** (`logs/k311/gen_r0.32_c*.log`). The pools lived in the scratch directory
  (`SP=/tmp/k311`) and were probably lost with the container. **Check first; regenerate if missing** (4 chunks × 10
  min × 4 cores; each chunk finds 60–70% of the orbits).
- **Production runner** (restart-safe; generation → slot pools → slot join, resuming from its log):
  `logs/k311/runners/run_prod.sh RHO NCHUNKS L LA K BITS`, e.g. `run_prod.sh 0.32 4 <L> <LA> 8 2.75`. Copy the runners
  into the scratch dir first, as the script expects.
  - The join command inside is
    `scripts/kt_join_slots.py 311 5 21,21,19,1 SLOTS OUT --k 8 --mem 3e9 --bits 2.75`.
  - 2.75 bits per remaining class is the partial-join value at 311.

## Tools

| tool | purpose |
|---|---|
| `scripts/flat_probe.py` | generation of flat orbit representatives, partial joins, n′ |
| `scripts/flat_capture.py` | capture–recapture flat-set size |
| `scripts/build_slot_pools.py`, `scripts/kt_join_slots.py` | 64-bit masks, r ≤ 64 orbits: 311 |
| `scripts/flatw.py`, `flat_gen_w.py`, `build_slot_pools_w.py`, `kt_join_slots_w.py` | two-word masks, r ≤ 128: 223, 491, 283, 859. Selftests pass. |
| `scripts/s1_nprime.py`, `scripts/k883_pjoin.py`, `scripts/flat_floor.py`, `scripts/exact_flat_count.py` | n′ and floor measurements |
| `kaggle/kt_join_kaggle.py` | the production k-tree that built H(3532) (env-configured) |
| `scripts/verify_h3532.py`, `pub/H3532/` | verifier and bundle patterns |
| published cores (`data/core_v223_b0.txt`, `data/seq_v283_g.txt`, `forensics/originals/`) | for planted tests and for comparing found cores |

**Multiplier closure [P]:** x ↦ x(g^s·y) rotates the e = (v − 1)/d orbit bits, keeps the row sum and F. So flat sets
are unions of rotation orbits, and one generated block gives up to e flat blocks.

## Tasks, in order

1. **Rung 1 (311): finish it.**
   - Regenerate pools if needed, build the slot pools, run the join, and log per-bucket progress and hits.
   - Any solution: verify exactly, then test it for inequivalence to Alpöge's 311 core (multipliers, signs, block
     swaps, reversal).
   - Compare hits per bucket with the logged expectation, and record the measured bits per remaining class (as was
     done at 883).
   - **If no hit after about 2× the expected work:** diagnose with planted tests. Plant Alpöge's 311 blocks in the
     pools and process their bucket (`--only-tau`, as done at 223). Decide whether the miss is pools (coverage) or
     join (bucket order) before spending more.
2. **Rung 2 (223).**
   - The planted 223 core was **not among the first 352 best-first buckets** on small ρ = 0.36 pools
     (`logs/k223/`). Measure where real solutions rank under the current bucket score.
   - If they rank low, fix the ordering before production: score by the expected tuple count in the bucket, using the
     calibrated per-class bits, not by key frequency alone.
   - Then measure ρ, flat-set size, n′ and pool size as at 311, and run.
3. **Rung 3 (491): disk-backed pools and join passes** for 2^27–2^30 blocks per row sum (`np.memmap`, partitioned
   passes, as the 523 kernels did for joins).
   - Measure n′(ρ) and the flat-set size before running.
   - If the 491 cost extrapolated from rungs 1–2 is more than about 2^47 operations on this machine, stop and report
     instead of running.
4. **Decision on 859/13.**
   - Using the rungs' measured n′(ρ), flat-set sizes and per-class bits, estimate 859/13's cost: flat-set size at the
     ρ that a solution needs, pool sizes, buckets, time.
   - **Run it only if it is within reach** (≤ about 100 core-hours here). Otherwise write down exactly what it
     needs: pool size, disk, CPU or GPU hours, a Kaggle plan. The repo has Kaggle notebooks in `kaggle/`.
5. **Rung 4 (283),** only as a measurement, unless it is cheap.
   - A converging flat-set count at ρ ≤ 0.30 (the old Wang–Landau runs did not converge below ρ ≈ 0.34;
     `report/REPORT.md` §0 item 5).
   - That decides whether a plain k-tree could have produced Alpöge's 283 at all, or whether his join was better.

## Deliverables

- **`report/b0_flatpool_ladder.md`:**
  - per rung: ρ, flat-set size, n′ (Gaussian and partial-join), pool sizes, buckets, time, hits, measured per-class
    bits, and verified cores with their inequivalence check;
  - the extrapolation to 859/13 and the go/no-go.
- **Any new open-order matrix:** an independent verifier `scripts/verify_h<N>.py`, a bundle `pub/H<N>/` (sequences,
  normalized matrix, standalone `verify.py`, as in `pub/H3532/`), and a commit at once.
