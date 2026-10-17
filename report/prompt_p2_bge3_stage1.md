# Prompt (Priority 2, maximum effort): stage-1 counting and sampling for the large b ≥ 3 configurations, then stage 2 where it is rich and cheap

*Self-contained, for a fresh maximum-reasoning Claude Code session on `middletownbooks/hadamard`, branch
`claude/focused-planck-j1645f`.*

- **Read first, in this order:**
  1. this file;
  2. `report/alpoge_gap_analysis.md` §4;
  3. `report/b599_answer.md` (how a "+48 bit" stage 1 turned out to be 2^4.9);
  4. `report/paley6_answer.md` (exact compressed stage-2 screens that killed 120 million stage-1 solutions);
  5. `report/alpoge_design_answer.md` §1–3 and §11 (the staged b ≥ 3 systems, borders, cost model);
  6. the progress log at the end of `report/prompt_flatpool_stage1.md` (item 2 entries).
- **A separate high-effort session works on the b = 0 flat pools** (`report/prompt_p1_flatpools.md`). Do not edit
  its files (`scripts/flat*`, `scripts/*slot*`, `logs/k311/`, `logs/k223/`).

## Rules

- **Pushing:** other sessions push to this branch. Before every push run
  `git fetch -q origin claude/focused-planck-j1645f && git merge -q --no-edit origin/claude/focused-planck-j1645f`,
  then `git push -u origin HEAD:claude/focused-planck-j1645f`.
- **Commits:** trailers follow the existing commits; no pull requests.
- **Processes:** never kill processes you did not start. Never `pkill -f` with a pattern matching your own shell;
  kill by PID.
- **Restart safety and check-ins:** write runners that resume from disk. Arm a check-in (`send_later`) whenever
  workers run. Commit as results arrive.
- **Check before writing:** `theory/ab_s1_cpsat.py` and `theory/ab_s1_sample.py` already exist, and a draft once
  overwrote one of them. Check `git log` and the existing files before creating new ones.
- **Tags:** [P] proofs, [M] exact computations (name the script), [E] estimates (name the estimator and its
  calibration).
- **Verification:** a construction counts only if verified exactly by an independent script.
- **Compute:** 4 cores and about 15 GB of RAM. Install `ortools`, `numba` and `sympy` with pip if missing.

## The problem

**Goal:** a Hadamard matrix H(4v) for an open v via a b ≥ 3 core, or exact obstructions that close the remaining b ≥ 3
configurations.

**Open v:**
- below 3000: 537, 571, 599, 643, 659, 669, 719, 739;
- 3000–4000: 751, 789, 823, 839, 859, 863, 907, 947, 955, 971, 991.

The repo has built 2092, 2292, 2588, 2884, 3532, 3668, 3676 and 3732. 919 and 883 are built; drop them as targets.

**b ≥ 3 construction** (Alpöge's `u`/`v`/`w` builders; `report/bge3_construction.md`):
- c = b + 1 is even with c | v + 1, h = (v + 1)/c − 1, G is abelian of order m = c·h, H ≤ G has index c, and
  Q = G/H.
- **Core:** Σ_k PAF_k = 4m at 0, −4b on H \ 0, +4 off H.
- **Border:** a frame of an H(4c) inflated by h, with the forced corner C = −K·Cb(z)ᵀ·Kcᵀ/(4c) required to be ±1.
  - c = 4, 6, 8: published frames transfer to every open order with the same c, and Singer and monomial classes
    exist.
  - c ≥ 10: frames exist (via H₂ ⊗ H(2c)) but no flat corner is known.

**GP-type staged slices** (`report/alpoge_design_answer.md` §1):
- G = G′ × Z_q with q prime, Z_q ≤ H, H = H′ × Z_q, and a slice group U × M (M of index s on F_q^*).
- **Stage 1** is the Z_q compression y = α + f·σ on G′, with
  Σ_k PAF_{G′}(y_k) = 4m − 4b(q − 1) at 0, −4bq on H′ \ 0, and +4q off H′.
- **Stage 0** is its Q-compression z, with Σ_k |ẑ_k(ψ)|² = 4b for ψ ≠ 1 and Σ a_k² = 4v.
- **Stage 2** is the cyclotomic norm system, solved by the hierarchical S/D solver (`theory/ab_prod347.py`,
  `ab_prod_generic.py`).
- **Positive control:** a new H(1388) found in 347's own slice in 28 minutes (D_def 16.2).

**What is known [M]:**
- Every configuration with ≤ 26 G′-orbits that the exact batches finished is empty, except 659's large-index
  single-stage slices, whose stage 2 is out of reach.
- All 61 Paley-type c = 4 configurations with ≤ 16 G′-orbits are empty.
- **117 configurations with more than 26 G′-orbits were skipped** (`logs/ab/c4_batch.jsonl`,
  `c6_bordered.jsonl` and `c8_bordered.jsonl`, entries with `skipped: "L1 > …"`). By order and c:

  | order, c | configurations |
  |---|---|
  | 643, c = 4 | 29 |
  | 839, c = 8 | 28 |
  | 599, c = 6 | 25 (two analysed in `b599_answer.md`) |
  | 719, c = 6 | 19 |
  | 919, c = 8 | 8 (order already built; low priority) |
  | 823, c = 8 | 6 |
  | 739, c = 4 | 2 |

**Why estimates cannot be trusted here [M].**
- Whole-slice and typical-y stage-1 estimates sample y with the lag-0 energy only. At 599 that evaluated a far lower
  tail: the Q-compression energy had mean 3,588 against a target of 416. "+48 bits" became 2^4.9.
- The z-conditioned estimator (`theory/ab_s1_zcount.py`) counts y exactly per z-quadruple for the Q-level and lag 0,
  and uses a lattice Gaussian for the remaining lags. On 347 it is about 2.3 bits low in total. **Per z it is
  unreliable**: 7 of 12 random z had exactly 0 solutions where it predicted 2^4.6–2^5.3. Where solutions are scarce it
  undercounts by 1.3–5.5 bits, because they are mostly symmetric blocks (`b599_answer.md` §3).
- **Stage-2 Gaussians are unreliable too.** At 991 ⟨4⟩ one predicted 2^76.7 lifts at the Z₁₃ level where the exact
  count was 0. The exact compressed stage-2 screens of `paley6_answer.md` §3 decide this (the Z_q level per
  Q-character, joint Z₂ × Z_q, and Q × Z_q). Always run them on a y before any stage-2 search.
- **Cost model for stage 2** (calibrated at 347 and 523): a core needs about 2^{D_def} S-solutions, where
  D_def = C(D|S) − h_D, plus the S-join and D-check costs (`report/paley_523_answer.md` §3.6;
  `report/alpoge_design_answer.md` §11). Budget: ≤ about 2^45 per core.
  - Example of a rich-but-too-expensive survivor: 599 H′ = Z₃², −1 on one axis, with stage 1 2^45.7 but stage 2 at
    2^58–2^61 per core.

## Tools

| tool | purpose |
|---|---|
| `theory/ab_s1_zscreen.py OUT.jsonl [--K 12] [--N 400] [--secs 900]` | **written, never run**: batch z-conditioned stage-1 screen over the 117 skipped configurations, restart-safe |
| `theory/ab_s1_zcount.py`, `theory/b599_flevels.py` | z-conditioned counts (`--calibrate N`) |
| `theory/ab_s1_cpsat.py` (`--decide`, `--enumerate N`), `theory/ab_s1_sample.py` (random stage-1 solutions for a prescribed z-quadruple), `theory/ab_s1_chars.py` (exact stage 1 by the character decomposition over Q, c = 4) | stage-1 solvers |
| `theory/ab_s1sd.py`, `theory/ab_sweep.py`, `theory/ab_cbatch.py`, `theory/ab_exact_small.py` | exact stage 1, slice build, batches, brute-force validation |
| `theory/ab_s2_zq_screen.py`, `ab_s2_z2zq_screen.py`, `ab_s2_qzq_screen.py`, `ab_s2_hzq.py`, `ab_full_chars.py` | exact compressed stage-2 screens of `paley6_answer.md` |
| `theory/ab_gauss.py` (`stage2_typical(sl, Y)`), `ab_typical_levels.py`, `ab_levels*.py` | stage-2 Gaussians and level budgets (D_def) |
| `theory/ab_prod347.py`, `ab_prod_generic.py`, `ab_stage2.py`, `ab_small_exact.py`, `ab_s2_stream.py` | stage-2 solvers (347-validated) |
| `theory/ab_borderfit.py`, `ab_assemble.py`, `ab_mono_border.py`, `ab_singer.py` | border fit and assembly; bundle writer from commit d7c7c6a; `scripts/verify_packed_h.py` |
| `theory/b599_symblock.py` | exact search restricted to symmetric blocks (scarce solutions tend to be symmetric) |
| `logs/s1z/` | the 347 smoke test and calibration logs |

## Tasks, in order

1. **Calibrate the z-conditioned count where exact counts are possible.**
   - Use 347's slice (all z, not 12), the small c = 6 controls (v = 35, 47), the c = 4 calibration slices in
     `logs/ab/paley6/calib/`, and same-shape rungs of the large configurations (smaller q, same G′ and U).
   - Measure the error of the total, not per z, and how it depends on the number of z-quadruples and on the share of
     symmetric solutions.
   - Fix the sampling: random rather than first-found z-quadruples.
   - Decide what margin a "rich" verdict needs.
2. **Run `ab_s1_zscreen.py` over the 117 configurations** with the calibrated margin.
   - Record per configuration: exact emptiness (no compression or no z-quadruple) [M], or the total [E] with its
     margin.
   - Add the symmetric-block exact search (`b599_symblock.py`) for configurations whose total is small but nonzero:
     scarce solutions are usually symmetric.
3. **For each rich configuration:**
   - typical-y stage-2 budget (D_def, cost per core);
   - drop those over 2^47;
   - for the rest, sample bordered stage-1 solutions (`ab_s1_sample.py` per bordered z-quadruple, or CP-SAT
     `--enumerate` with random objectives);
   - run the exact compressed stage-2 screens on every y;
   - rank the survivors by stage-2 yield at the actual y.
4. **Stage 2 on the best:**
   - calibrate S-solution counts and D-success on the first y's, as at 523 and 347;
   - cost it, then run production, restart-safe;
   - any core: assemble with the fitted border and forced corner, verify independently, write the bundle
     `pub/H<N>/`, and commit at once.
5. **If nothing survives:**
   - write the cleanest obstructions;
   - say what the remaining b ≥ 3 frontier is: more than 60 G′-orbits, c ≥ 10 borders (the corner-flatness
     condition of `alpoge_design_answer.md` §3.1), or noncyclic Q and non-split H as in 359.

## Deliverables

**`report/bge3_stage1_answer.md`:**
- the calibration of the z-conditioned count, with numbers;
- a table of all 117 configurations (exact or estimated stage 1, stage-2 budget, verdict);
- the sampled y's, screened, with stage-2 results;
- any core (verified, bundled), or the obstructions and the frontier.

Logs go in `logs/s1z/`, scripts in `theory/ab_*` or `theory/s1z_*`.
