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
(`report/alpoge_gap_analysis.md` §4, `report/alpoge_design_answer.md` §11). The best costed lead is 599 at c = 6
(q = 11, s = 2, H' = Z₉, −1 on Z₉; 30 G'-orbits; stage 2 about 2^43 per core, D_def 21.5).

Existing machinery: `theory/ab_sweep.py` (slice build), `theory/ab_s1sd.py` (exact stage 1: α-quadruples →
S-join → D-join), `theory/ab_s1_chars.py` (stage 1 by the character decomposition over Q = Z₂², c = 4 only),
`theory/ab_stage2.py` and `theory/ab_prod347.py` (stage 2, validated on 347), `theory/ab_typical_levels.py`
(stage-2 budget per y), `scripts/cyclo_stage1_cpsat.py` (CP-SAT stage 1; found a 523 y in 323 s).

Steps:
1. A sampler of stage-1 solutions y (not an enumerator): CP-SAT per row quadruple and α-congruence class, or a
   randomised hashed join that keeps only a sample of each list; validate it on 347's slice (its published y must be
   reachable; compare the sampled y distribution with the exact stage 1 there).
2. Rank sampled y by the stage-2 Gaussian yield (`ab_typical_levels` evaluated at the actual y), as the 419 gadgets
   were ranked; calibrate the ranking on 347 (exact stage 2 per y is affordable there).
3. Apply to 599 c = 6, then to the unexamined configurations with the best whole-slice budgets.

## Progress log (newest last; the previous session appended here as it worked)

- 2026-10-17: task list written.
