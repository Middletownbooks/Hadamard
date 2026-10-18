# Prompt (maximum effort): a multi-level stage-2 solver, calibrated on 523 and 347, then costed and run on 599 c = 6 ×4 on Z₉

*Self-contained, for a fresh maximum-reasoning Claude Code session on `middletownbooks/hadamard`, branch
`claude/focused-planck-j1645f`.*

- **Read first, in this order:**
  1. this file;
  2. `report/bge3_stage1_answer.md` §0, §1, §6 and §8 (why this slice is the best b ≥ 3 option, and its numbers);
  3. `report/b599_answer.md` §5 (the nested-level cost formula, worked out on the Z₃²-a neighbour);
  4. `report/paley_523_answer.md` §3 (the two-level Wagner tree, D-check and cost model that built H(2092));
  5. `report/alpoge_design_answer.md` §1, §3 and §5 (b ≥ 3 staged systems, borders, the 347 control that built a new
     H(1388)).
- **Another session may be running the b = 0 ladder** (`logs/k283/`, `scripts/flat*`). Do not touch its files or
  processes.

## Rules

- **Pushing:** other sessions push to this branch. Before every push run
  `git fetch -q origin claude/focused-planck-j1645f && git merge -q --no-edit origin/claude/focused-planck-j1645f`,
  then `git push -u origin HEAD:claude/focused-planck-j1645f`.
- **Commits:** trailers follow the existing commits; no pull requests.
- **Processes:** never kill processes you did not start. Never `pkill -f` with a pattern matching your own shell;
  kill by PID.
- **Restart safety and check-ins:** write runners that resume from disk. Arm a check-in (`send_later`) whenever
  workers run; an idle session was reclaimed once and killed production workers. Commit results as they arrive.
- **Check before writing:** look at `git log` and the existing files first. Several "missing" tools already exist,
  and a draft once overwrote one.
- **Tags:** [P] proofs, [M] exact computations (name the script), [E] estimates (name the estimator and its
  calibration).
- **Verification:** a construction counts only if verified exactly by a script sharing no code with the search.
- **Compute:** 4 cores, about 15 GB of RAM, about 30 GB of scratch disk. Install `numba`, `ortools` and `sympy` with
  pip if missing.

## 0. Goal and deliverables

H(2396) (v = 599, open) through the b = 5, c = 6 Goethals–Seidel construction with an inflated border, in the slice
below.
- Its stage 1 is rich: about 2^30.8 solutions per symmetry class [E].
- About 2^15–2^19 cores per symmetry class are expected [E].
- What is missing is a stage-2 solver about 3–6 bits cheaper than the current model: **2^49.9–2^52.6 per core against
  a 2^45–2^47 budget.**

Deliverables, in order of value:
- **(A)** A verified H(2396): an independent verifier, plus a bundle `pub/H2396/` (core on G, `border.json` with
  layout, K, Kc and the corner, the normalized matrix, a standalone `verify.py`; the b ≥ 3 bundle writer is from commit
  d7c7c6a).
- **(B)** A multi-level stage-2 solver with its cost model validated against exact counts and wall-clock on 523 and
  347, then a calibrated cost per core at 599 ×4, and production if it is within about 2^47 (about 200 core-hours).
- **(C)** If it does not fit: the measured floor, with what it would take (compute, or which level is missing).

Write to `report/b599x4_answer.md`, with scripts in `theory/m2_*.py` (solver) and `theory/b599x4_*.py`, and logs in
`logs/m2/` and `logs/b599x4/`.

## 1. The slice (exact)

**Parameters:** v = 599, b = 5, c = 6, h = 99, m = 594.
- G = Q × H′ × Z₁₁ with Q = Z₂ × Z₃ and H′ = Z₉; H = H′ × Z₁₁; G′ = Q × H′ has order 54.

**Slice group:**
- Paley(11) on Z₁₁ (q = 11, s = 2, f = 5);
- U = ×4 on Z₉ (order 3; on Z₉ its orbits are {0}, {3}, {6}, {1, 4, 7}, {2, 8, 5});
- 30 G′-orbits and 90 orbits per block.
- Configuration key `[599, 6, 11, 2, [9], [[2, 3, 4]]]` (`theory/ab_sweep.py` format; mods [2, 3, 9, 11]; the
  U-entry acts on the Z₉ axis by ×4).

**Core condition** (b ≥ 3): Σ_k PAF_k = 4m at 0, −4b = −20 on H \ 0, +4 off H. The border is an inflated frame of an
H(24) with forced corner C = −K·Cb(z)ᵀ·Kcᵀ/(4c), which must be ±1.
- c = 6 borders: the 347 template class and the monomial/Singer class.
- At least 17 of this slice's compression multisets are bordered (`theory/s1z_border.py`).

**Stage 1** (Z₁₁ compression on G′): letters y = α + 5(p + n) ∈ {±1, ±9, ±11}, with
Σ_k PAF_{G′}(y_k) = 2176 at 0, −220 on H′ \ 0 and +44 off H′.
- Q-level: Σ_k |ẑ_k(ψ)|² = 20 for ψ ≠ 1, and Σ a_k² = 2396.
- Estimates: 2^61.1 ordered solutions, 2^+30.8 per symmetry class [E] (`theory/s1z_screen.py`; calibrated within
  ±3.6 bits where solutions exist).

**Stage 2** (index 2, q = 11 ≡ 3 mod 4) [P] (`report/alpoge_design_answer.md` §11.3). With u = p − n and
a = α − σ/2, given y it is a PAF problem on G′ for the ±2-valued u_k on the free points, plus a linear system
Σ_k CC_odd(a_k, u_k) = 0. The S/D split on Q's Z₂ applies, with an S level and a D level on Z₃ × Z₉.

**Typical-y budget [E]** (six y's; `logs/s1z/budget_seeds_599rich.log`):

| quantity | value |
|---|---|
| H₂ | 98–100 |
| C₂ | 111.4–112.4 (−10.4 to −14.4 bits per y) |
| C_S | 63.2–65.6 |
| H_S | ≈ 78.5 |
| h_D | 21–21.5 |
| D_def | 24.3–27.5 |
| **cost 2^{D_def + C_S − H_S/2}** | **2^49.9–2^52.6 per core** |

**Exact screens so far [M]:** all 85 compression multisets pass the Z₁₁ and Z₂ × Z₁₁ compression screens. The joint
Z₃ × Z₁₁ and Q × Z₁₁ screens did not fit in memory per compression quadruple, but **are cheap per y**.

## 2. The cost model and why a nested level should help

- **Model** (calibrated at 347 and 523): a core needs about 2^{D_def} S-solutions. Each costs
  2^{C_S − H_S/2} in a plain 2+2 S-join, plus a D-check of about 2^{h_D/2}.
- **523 (H(2092)):**
  - the realised S-step was 2^{F + 1 − n} ≈ 2^18.1 pair visits per S-solution, with a two-level Wagner tree (Z₃
    middle level on parents, then the order-9 equations);
  - the D-check took 1.5–1.8 ms;
  - D_def 26.1 and 2^26.1 S-solutions per core;
  - about 2^44.2 pair visits and about 80 core-hours per core, measured.
  - The ideal Wagner cost was 2^40.7, out of reach because a Z₃ middle value had only a few hundred parent pairs
    (§3.3 there).
- **347 (new H(1388)):** D_def 16.2, 11 middle values, 2^14.2 S-solutions, 28 minutes.
- **Nested level** (`report/b599_answer.md` §5): with a second level M, the cost per S-solution is

      2^(C_S − H_S + H_M/2) + 2^((C_S − C_M) − (H_S − H_M)/2),

  which helps exactly when M is constraint-dense (C_M > H_M/2).
- **599's S level lives on Z₃ × Z₉, with several U-stable candidate levels:**
  - Q's Z₃ compression;
  - the Z₉ → Z₃ quotient (compression along 3Z₉);
  - Z₃ × Z₃;
  - the ×4-orbit structure;
  - and their combinations, possibly a three-level tree.
  Measure (H_M, C_M) for each on real y's before choosing.

## 3. Plan

1. **Stage-1 supply.** No 599 ×4 y exists yet. On the similar 719 slice the α-list sampler and CP-SAT both failed
   (`report/bge3_stage1_answer.md` §7).
   - Build an S/D stage-1 generator per bordered compression multiset: S-lists of σ per block subsampled to fit
     memory, a hash join on the S-level equations, then a D-join. That is the design costed in `b599_answer.md` §5 at
     about 2^21–2^32 per y [E].
   - Alternatively, CP-SAT per prescribed compression quadruple (`theory/s1z_cpsat.py`), or a residue-aware 4-list
     join as at 523 (`theory/q523_ymitm.py`).
   - Validate on 347's slice (its published y must be producible) and on small rungs with exact counts
     (`theory/s1z_calib*.py`).
   - Each y found: run the joint Z₃ × Z₁₁ and Q × Z₁₁ screens at that y (cheap per y), and keep only survivors.
2. **Multi-level stage-2 solver**, general over G′ and the levels: a tree of middle levels, radix-partitioned
   cache-resident joins (`mjoin_rad2` in `theory/q523_wagner.py`), a compiled D-check (`theory/q523_dcheck_nb.py`),
   checkpoints and `--resume`. Start from `theory/ab_prod_generic.py`, `ab_stage2.py` and `ab_prod347.py`.
3. **Calibrate on 523 and 347 before touching 599:**
   - **523:** on the production y (`logs/q523/y523_best.json`), reproduce the record's S-solutions per middle value
     (K = 32 lists, §3.3 and §6 there) and the D-success (2^−26.0 per S-solution). Then measure what a nested level adds
     there: its pair visits per S-solution against 2^18.1.
   - **347:** its own slice (`logs/ab/cfg_347_s6.json`). Re-find cores with the new solver, and compare S-solutions per
     core (model 2^16.2) and wall-clock with the 28-minute run.
   - Also planted tests at full size on both (as `theory/q523_planted.py` and the 347 plants did).
   - **Record measured against predicted for every quantity.** The model is the acceptance criterion: if it misses
     523 or 347 by more than about 2 bits, fix the model before costing 599.
4. **Cost 599 ×4.**
   - On the first real, screened y's: measure (H_M, C_M) per candidate level, S-solutions per middle value, D-success
     per S-solution, and time per pair visit.
   - Pick the level tree and give the cost per core with its spread over y.
   - Also try **y-selection**: cost varies about 2.7 bits over typical y's, so rank y's by their own level budget
     before running.
5. **Production if ≤ about 2^47** (about 200 core-hours on 4 cores): restart-safe, with a tally of expected cores.
   - Any core: assemble with the fitted border and forced corner (`theory/ab_borderfit.py`, `ab_assemble.py`), check
     H Hᵀ = 2396·I independently, commit at once, then bundle.
   - Do not stop a calibrated run early: H(2092) took about 3 expected cores, a 5% wait.
6. **If it does not fit:** report the measured floor per level and what would close the gap (more cores, a GPU join,
   another level), and whether the Z₃² variants (6–10 bits dearer) or a different U on Z₉ could be cheaper with the new
   solver.

## 4. Tools (read before writing new code)

| tool | purpose |
|---|---|
| `theory/s1z_lib.py`, `s1z_screen.py`, `s1z_cpsat.py`, `s1z_border.py` | compression multisets, z-screen, CP-SAT per compression quadruple, border classes |
| `theory/s1z_qzq.py`, `s1z_typescreen.py`, `s1z_zqscreen.py`, `s1z_char_batch.py` | exact compression screens of stage 2 (per type, aggregated, per Q-character) |
| `theory/s1z_budget_batch.py`, `theory/ab_typical_levels.py`, `theory/ab_levels.py`, `theory/b599_cond.py` (`budget_sd`) | level budgets, D_def, the S/D budget of stage 1 |
| `theory/ab_s1sd.py`, `ab_s1_sample.py`, `ab_s1_cpsat.py`, `theory/b599_symblock.py` | stage-1 solvers |
| `theory/ab_prod347.py`, `ab_prod_generic.py`, `ab_stage2.py`, `ab_small_exact.py`, `ab_s2_stream.py` | the existing single-middle-level stage 2 (347-validated) |
| `theory/q523_wagner.py`, `q523_wagner2.py`, `q523_dcheck_nb.py`, `q523_dcheck2.py`, `q523_prod.py`, `q523_dfine_calib.py`, `q523_dsuccess_real.py`, `q523_planted.py` | 523 kernels and calibrations. Read only; copy into `theory/m2_*`. |
| `theory/ab_borderfit.py`, `ab_assemble.py`, `ab_template.py`, `ab_mono_border.py` | border fit and assembly |
| `scripts/verify_packed_h.py` | independent verifier; b ≥ 3 bundle writer from commit d7c7c6a |
