# Task list for a new instance: bigger, deeper flat pools (b = 0), then a stage-1 sampler (b ≥ 3)

Written 2026-10-17 for a fresh Claude Code session on this repo (the previous container may have been reset by usage
ceilings). Branch: `claude/focused-planck-j1645f`. Read this file first, then `report/alpoge_gap_analysis.md` (why
these two items), then the progress log at the end of this file, which the previous session kept up to date.

## Rules (from the standing research prompt; keep them)

- Other sessions push to this branch. Before every push: `git fetch -q origin claude/focused-planck-j1645f &&
  git merge -q --no-edit origin/claude/focused-planck-j1645f`, then `git push -u origin claude/focused-planck-j1645f`.
- Never kill processes you did not start; never `pkill -f` with a pattern that matches your own shell; kill by PID.
- Write restart-safe runners (resume from files on disk). Arm a check-in (`send_later`) whenever workers run, so the
  session never idles while they run. Commit results as they arrive; containers restart.
- A construction counts only if verified exactly. Separate theorems [P], exact computations [M], estimates [E].
  Calibrate every Gaussian against exact counts.
- Commit trailers follow the existing commits (`Co-Authored-By:` and `Claude-Session:` lines); no pull requests.
- Large regenerable files (pools, caches) stay out of git (see `.gitignore`); commit logs, small JSON and scripts.

## Priority 1: bigger, deeper flat pools for the b = 0 k-tree

Why: Alpöge's b = 0 cores (223, 283, 311, 487, 491) are local minima of their own PAF energy, like the blocks of
our flat-pool k-tree, so his b = 0 method is the same family; his pools were deeper (mean PAF²/v ≈ 0.08–0.10 at
223/283/311, i.e. flatness ratio ρ ≈ 0.29–0.32) and his slices larger (139–207 constraint bits, against our best 131)
(`report/alpoge_gap_analysis.md` §2). The goal is to stretch our method and measure how far it reaches; the
calibration ladder uses published orders as rungs, then targets an open b = 0 slice.

Existing machinery:
- `kaggle/kt_join_kaggle.py`: production k-tree (pools by tabu walks with F ≤ ρ²·m·L, best-first τ-buckets on k
  chunk classes, histogram class, 64-bit hash of the rest, bloom match, exact verification). Env overrides
  `KT_V, KT_D, KT_PROFILE, KT_RHO, KT_K, KT_LMAX, KT_GEN_MINUTES, KT_MEM_BUDGET_BYTES, KT_INSTANCE, KT_OUT`.
  Restart-safe (pools and prep cached under KT_OUT, finished buckets skipped).
- `scripts/s1_pools.py`, `scripts/s1_nprime.py` (Gaussian and partial-join n′ for a pool), `scripts/k883_pjoin.py`
  (partial joins on production pools), `logs/k883/README.md` (the H(3532) run: what worked, how pools saturate,
  how the per-bucket rate was calibrated).
- Constraint bits per slice: `scripts/b0_open_screen.py V1,V2 [N]`.

Ladder (all b = 0, multiplier slice U of order d, profile = row sums):

| rung | v, d | orbits / classes | n (untwisted) | flat n′ [E] | lists 2^(n′/3) | status |
|---|---|---|---|---|---|---|
| done | 883, 21 | 43 / 21 | 111.9 | ≈ 78 | 2^26 | H(3532) |
| done | 487, 9 | 55 / 27 | 130.9 | ≈ 83 | 2^28 | 2 new cores |
| 1 | 311, 5 | 63 / 31 | 139.4 | ≈ 84 | 2^28 | see progress log |
| 2 | 223, 3 | 75 / 37 | 156.6 | ≈ 90 | 2^30 | our 523 runs (≈ 91) failed |
| 3 | 491, 7 | 71 / 35 | 169.1 | ≈ 106 | 2^35 | needs disk-backed pools |
| 4 | 283, 3 | 95 / 47 | 207.0 | ≈ 119 (measured) | 2^40 | Alpöge's hardest b = 0 |

Steps for each rung:
1. Admissible profiles: row sums s ≡ ±1 (mod 2d), Σ s² = 4v (311: (21, 21, 19, 1) [Alpöge's] and (21, 19, 19, 9)).
2. Flat-set size versus ρ: generation yield and saturation (`logs/k883/README.md`: pools saturate; reversal closure
   x(y) ↦ x(−y) measures coverage), at ρ ∈ {0.32, 0.36, 0.40, 0.45}.
3. n′(ρ) by partial joins (`scripts/s1_nprime.py`, `scripts/k883_pjoin.py`), against the Gaussian.
4. Choose ρ maximising expected solutions in pools, L⁴·2^(−n′(ρ)), subject to memory (15 GB here: about 2^25–2^26
   blocks per row sum with the current join).
5. Run `kaggle/kt_join_kaggle.py` under a restart-safe runner; log to `logs/k<v>/`; commit logs as they arrive.
6. Any solution: verify exactly (integer PAF, H Hᵀ = 4v·I), compare with the published core (inequivalence under
   multipliers, sign, block swaps), record the cost against the ladder.

Method improvements to make, in this order:
- (a) Multiplier closure of pools: the flat set is invariant under the e = (v − 1)/d multipliers (they permute
  orbits and, uniformly, lag classes) and keeps the row sum, so each generated block gives up to e flat blocks.
  This fills a flat set to saturation at a fraction of the generation time, which matters most for deep (small ρ)
  pools, whose walks rarely hit.
- (b) Disk-backed pools and join passes for 2^27–2^30 blocks per row sum (rungs 3–4).
- (c) A converging measurement of the flat-set size at 283 (the Wang–Landau runs of `REPORT.md` §0 did not converge
  below ρ ≈ 0.34), which decides whether rung 4 is reachable at all.

## Priority 2: a stage-1 sampler for b ≥ 3 configurations with 30–60 G'-orbits

Why: at the open orders every b ≥ 3 configuration with ≤ 26 G'-orbits that our exact batches finished is empty, and
117 of 270 configurations (more than 26 orbits) were never examined, because exhaustive stage 1 does not scale
(`report/alpoge_gap_analysis.md` §4, `report/alpoge_design_answer.md` §11). **Read `report/b599_answer.md` first:**
the 599 c = 6 lead (H' = Z₉, −1 on Z₉, 30 G'-orbits) is dead: conditioned on its 420 admissible Q-compression
quadruples its stage 1 holds about 2^4.9 ordered solutions (whole-slice "+48" was an artefact of energy-only
sampling); its neighbour H' = Z₃² with −1 on one axis has a rich stage 1 (2^45.7) but stage 2 at 2^58–2^61 per
core. Whole-slice and typical-y estimates of stage 1 are not to be trusted; count conditioned on z.

Existing machinery (check `git log` before writing new files: several of these already do part of the job):
- `theory/ab_s1_cpsat.py`: CP-SAT stage 1 for these slices (orbit variables, pair products, implied block-sum and
  Q-character constraints, symmetry breaking); `--decide` and `--enumerate N` (every solution verified exactly).
- `theory/ab_s1_sample.py`: random stage-1 solutions for a prescribed Q-compression quadruple (random α-triples
  completed by a hash lookup of the fourth, then the exact S/D joins); used for the monomial c = 6 pipeline.
- `theory/ab_sweep.py` (slice build), `theory/ab_s1sd.py` (exact stage 1), `theory/ab_s1_chars.py` (c = 4 only),
  `theory/ab_stage2.py` and `theory/ab_prod347.py` (stage 2, validated on 347), `theory/ab_typical_levels.py` and
  `ab_gauss.stage2_typical(sl, Y)` (stage-2 Gaussian at a given y), `theory/ab_borderfit.py` (border fit of a y),
  `scripts/cyclo_stage1_cpsat.py` (CP-SAT stage 1 for cyclic GP slices; a 523 y in 323 s).
- 599 c = 6 configuration (dead, see above; useful as a test case): `{"v":599,"c":6,"b":5,"h":99,"m":594,"Q":[2,3],"H1":[9],"q":11,"s":2,"f":5,
  "mods":[2,3,9,11],"nQ":2,"U":[[2,2,8]]}` (30 G'-orbits; `logs/ab/p1/batch/levels_batch_p*.jsonl`).

Steps:
1. Sample stage-1 solutions y with the two tools above (CP-SAT with random objectives for diversity; or
   `ab_s1_sample.py` per border-class compression quadruple); validate on 347's slice (the published y must satisfy
   the model up to the model's symmetry breaking: positive first letter per block, non-decreasing block energies).
2. Before sampling anywhere: a z-conditioned stage-1 count (`theory/b599_flevels.py`) for each unexamined
   configuration (117 with > 26 G'-orbits; `logs/ab/` batch logs list them), to find rich stage 1 with cheap stage 2.
3. Rank sampled y by the stage-2 Gaussian yield (`ab_gauss.stage2_typical` at the actual y), as the 419 gadgets were
   ranked; calibrate the ranking on 347 (exact stage 2 per y is affordable there); run stage 2 on the best.

## Progress log (newest last; the previous session appended here as it worked)

- 2026-10-17: task list written.
- 2026-10-17, rung 1 (311) probes [M] (`scripts/flat_probe.py`, logs `logs/k311/probe_rho*_s*.log`):
  - rho = 0.40: 1 min on 4 cores gives 3.8 M blocks per row sum, nearly all in distinct multiplier orbits (coverage of
    the rotation closure 1.6%, reversal closure 0%): the flat set is far larger than any sample. Gaussian n' = 89.9
    (3.0 bits per class); partial joins to k = 10 sit 0.035 bits per class above it, so n' ≈ 91.
  - rho = 0.36: 0.6 M per row sum per minute, closure coverage 1.7%; n' = 86.0 Gaussian, ≈ 87 with the partial-join
    offset (+0.031 per class).
  - New tools: multiplier closure (rotations of the e = 62 orbit bits keep F and the row sum [P]; checked
    numerically), `scripts/build_slot_pools.py` (slot A = one block per multiplier orbit, slots B, C, D = samples of the
    closures: symmetry breaking, so hits are distinct cores), `scripts/kt_join_slots.py` (the kaggle join with one pool
    per slot), `scripts/flat_capture.py` (capture-recapture flat-set size from two seeds).
  - Choice of rho: distinct cores in the slot pools = |A| L^3 2^(-n'), |A| growing with generation time; rho = 0.36
    beats 0.40 by about 2^5 at equal time. The production join on 64-bit masks covers r <= 64 orbits only (311 has
    63); rungs 2-4 (r = 71, 75, 95) need two-word masks.
- 2026-10-17: caution for the next session: `theory/ab_s1_cpsat.py` and `theory/ab_s1_sample.py` already existed;
  a draft of mine briefly overwrote the first and was reverted from git (`git checkout`). Check for existing files.
- 2026-10-17, rung 1 production at rho = 0.32 (6 chunks x 10 min, `logs/k311/gen_r0.32_c*.log`):
  - The flat set saturates at this depth: chunk 1 found 135k-238k orbits per row sum, chunk 2 overlaps 60-70% of them.
    Capture-recapture (chunks 1, 2; `logs/k311/capture_rho0.32_c1c2.txt`): 196k / 334k / 401k orbits for row sums
    1 / 19 / 21, i.e. 2^23.5 / 2^24.3 / 2^24.6 blocks after closure (lower estimates).
  - So at rho = 0.32 pools of 2^24 per slot are essentially the whole flat set. Solution tuples in pools ≈
    2^(95.5 - 82) ≈ 2^13 (each core appears in several rotated copies, which land in different buckets). rho = 0.36
    with pools of 2^25 gives a similar count at twice the pool size; 0.32 was kept.
- 2026-10-17, item 2: `theory/ab_s1_zcount.py` (new): stage-1 count conditioned on the Q-compressions for any
  configuration of `theory/ab_sweep.py` with nQ >= 1 (generalises `theory/b599_flevels.py`): exact per-block counts by
  DP over Q-coset components, compressions within the Parseval bound around each admissible row-sum mean,
  z-quadruples by a 2+2 join on (PAF_Q, |z|^2), exact y-quadruple count per z-quad (Q-level and lag 0), lattice
  Gaussian of the remaining lags from uniform samples given (z_k, E_k). Smoke test on 347's slice
  (`logs/s1z/zcount_347_smoke.log`): 222 compressions per block, 34,560 ordered z-quadruples, log2 N = 32.4 and
  about +2.3 bits per z-quad after the 6 remaining lags. NOT yet calibrated: next step is an exact count per
  z-quadruple (S1SD with zallow = that quadruple's vectors, or `theory/ab_s1_sample.py`) on small configurations
  (the c = 6 controls v = 35, 47; the c = 4 calibration slices in `logs/ab/paley6/calib/`), and random rather than
  first-found z-quadruples in the total.
- 2026-10-17, item 2 calibration [M vs E] (`theory/ab_s1_zcount.py --calibrate 12`, `logs/s1z/zcount_347_calib.log`):
  on 347's slice, 12 random z-quadruples counted exactly (every block completion enumerated, energy-split ejoin4):
  7 have 0 solutions where the Gaussian predicts 2^4.6-2^5.3 (exact obstructions it cannot see); 5 have 2^5.2-2^10.7
  where it predicts 2^5.3-2^7.9 (0-3 bits low). Averaged over the 12 the slice total is about 2.3 bits low. Use the
  totals to screen configurations (rich versus poor stage 1, to within a few bits); never trust a single z.
- 2026-10-17, rungs 2-4 tooling (two-word masks, r <= 128): `scripts/flatw.py` (kernels; selftest against the 64-bit
  path and direct FFT PAFs), `scripts/flat_gen_w.py` (generation to multiplier-orbit representatives),
  `scripts/build_slot_pools_w.py`, `scripts/kt_join_slots_w.py`. Smoke test at 223 (d = 3, 75 orbits, rho = 0.36,
  profile (17, 17, 17, 5)): blocks flat and canonical, join runs end to end. Before any 223 production: a planted test
  (put Alpöge's 223 blocks, `data/core_v223_b0.txt`, profile (23, 11, 11, 11), into the pools and process its bucket).
- 2026-10-17, planted test of the two-word join [M]: Alpöge's 223 core (natural coordinates via
  `forensics/near_structure.load_all`, slots (23, 11, 11, 11)) placed in small rho = 0.36 pools; processing only its
  bucket (`scripts/kt_join_slots_w.py --only-tau -2,6,-6,2,-6,-2,6,6`) returns it as a verified solution
  (`logs/k223/planted_core_*`). Note: that bucket was not among the first 352 in best-first order on these pools.

