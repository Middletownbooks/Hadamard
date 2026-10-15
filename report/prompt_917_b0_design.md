# Prompt: H(3668) at v = 917 by the b = 0 Goethals–Seidel array: a split dual level for the index-26 slice (self-contained)

*For a fresh Claude Code session on `middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`. Budget:
about 150 core-hours in all, on the session's 4 cores.*

- **Tags.** [M] exact computation in the repo (name the script); [P] proved (give the proof); [E] Gaussian estimate
  (name the estimator and its calibration); [L] literature recollection (unchecked: arXiv and OEIS are blocked).
- **Branch.** Other sessions push to it; another session is working on 917 by a different route (c = 6, b = 5) at
  the same time. Fetch and merge before every push. Never kill processes you did not start, and never `pkill -f`
  (or `pgrep -f` in a wait loop) with a pattern that matches your own shell; kill by PID.
- **The container is reclaimed when the session is idle,** which kills every background process (it happened
  twice). While workers run, keep the session active (blocking waits of at most 10 minutes) and also arm a check-in
  (`send_later`); write restart-safe runners that resume from their logs. Cap memory per process
  (`resource.setrlimit(RLIMIT_AS, ...)`, about 3.5 GB each for 4 workers): an uncapped worker was OOM-killed at 14 GB.
- Commit and push as results arrive. Commit trailers follow the repository's existing commits; no pull requests.

## 0. Goal and deliverables

v = 917 = 7·131 is open (4v = 3668). Build H(3668) as the Goethals–Seidel array (b = 0) of four ±1 circulants of
length 917 with Σ_k PAF_k = 3668·δ, in the slice below, where the Gaussians say cores are abundant but no level
decomposition found so far is cheap.

Deliverables, in order of value:
- **(A)** A verified H(3668): H Hᵀ = 3668·I exactly, checked by an independent script that shares no code with the
  search (pattern: `scripts/verify_h3732.py`, which also checks the slice structure with `--slice q s`), and a bundle
  `pub/H3668/` made like `pub/H3732/` (`scripts/make_bundle_b0.py`).
- **(B)** A method validated end to end on the same shape at a smaller order (§3), costed per core, with the run itself
  if the cost is within budget.
- **(C)** Exact obstructions if the slice (or a sub-slice you choose) is empty.

Write `report/b0_917_answer.md`; scripts `theory/b917_*.py`; logs `logs/b917/`.

## 1. Where the project stands

**Built here, each verified three ways:** H(2092) (v = 523, b = 1), H(2292) (573, b = 1), H(2588) (647, b = 1,
doubly staged), H(2884) (721, Cooper–Wallis), and **H(3732) (v = 933, b = 0, this method's ancestor)**:
`report/alpoge_design_answer.md` §0, §4; `pub/H3732/`.

**The b = 0 staged method that built 933** (`theory/ab_b0z3.py`, `report/alpoge_design_answer.md` §4, §6):
- G = Z_F × Z_q (cyclic by CRT). The slice: each x_k is constant on {g} × C_j (or on U-orbit × C_j), where the C_j
  are the s cosets of the index-s subgroup M of F_q^* (|M| = f), plus the point 0 of Z_q.
- **Stage 1** (compression along Z_q onto Z_F): y_k(g) = α_k(g) + f·σ_k(g), Σ_k PAF_{Z_F}(y_k) = 4v·δ. Exact and
  small (`theory/ab_b0_stage1.py`).
- **Dual** (compression along Z_F onto Z_q): the class values T_{k,j} = Σ_g x_k(g, C_j), Σ_k PAF_{Z_q} = 4v·δ, one
  equation per lag class.
- **Cross**: every character nontrivial on both factors.
- At 933 (F = 3, q = 311, s = 10): the dual vectors (4^10 per block) are listed and joined; then each dual solution's
  cross level is an exact 2 + 2 join over its fibres.
- Cost rule [E, checked at 933]: a dual solution with fibre product 2^F yields about 2^(F − C_cross) cores, and its cross
  join costs at least 2^(F/2). At 933: C_cross ≈ 42, F ≈ 34.6, cost 2^24.7 per core (measured 2^25.7, 3 minutes).
- The cross-level Gaussian (lattice-exact, `theory/q523_gauss.gauss_bits`, `theory/ab_b0z3_calib.py`) is within about 1.5
  bits of exact core counts at q = 71 and 311. b = 0 Gaussians have been reliable here because the dual compression
  is a sum of several signs per class, close to continuous. Doubly staged b = 1 slices, where they failed, are a
  different situation (`report/next_design_answer.md` §0.5).

**Why 917 is next** (`report/alpoge_design_answer.md` §6; `logs/ab/b0_rank2.log`):

| slice | stage 1 [M] | stage 2 per y [E] | remark |
|---|---|---|---|
| q = 131, s = 26 (f = 5), U = ⟨2⟩ on Z₇ (order 3) | **38,380 y's** | **H₂ 278, C₂ 225: +52 bits per y** | the target |
| q = 131, s = 13, U = ⟨2⟩ | 96 | −27.8 | dead [E] |
| q = 131, s = 26, U = Z₇^* | 5 | −411 | dead [E] |
| q = 131, s = 13, U = −1 | **0** | – | empty [M] |
| q = 131, s = 2, 5; s = 10 with U ≠ 1 | screen A fails | – | empty [M] |
| q = 131, s = 10 or 13, U = 1 | not enumerated | fibres ≈ 2^50 per block | cross join ≈ 2^90 per dual solution |

So the target slice holds about 2^15.2 · 2^52 ≈ 2^67 cores [E]. The problem is finding one, not existence.

The cfg (JSON, as used by `theory/ab_sweep.build` and `theory/ab_s1sd.cfg_from_json`):
`{"v":917,"c":1,"b":0,"h":null,"m":917,"Q":[],"H1":[7],"q":131,"s":26,"f":5,"mods":[7,131],"nQ":0,"U":[[0,3,2]]}`
(`U = [[axis, order, generator]]`: the multiplier 2 of order 3 on axis 0 = Z₇).

## 2. The structure of the target slice [P unless marked]

Rows: the U-orbits of Z₇ are R₀ = {0}, R₁ = {1, 2, 4}, R₂ = {3, 5, 6} (weights 1, 3, 3). Columns: {0} and the 26
cosets C_j = γ^j M of the order-5 subgroup M of F₁₃₁^* (γ primitive). Since f = 5 is odd, −1 ∉ M and
−C_j = C_{j+13}. Each block is an 81-cell ±1 array x_k(R_i, C_j) (with the 0 column).

Characters: χ trivial, or in the orbit {χ₁, χ₂, χ₄} (conjugate orbit {χ₃, χ₅, χ₆}); ψ trivial or one of 26 classes.
With η = Σ_{g∈R₁} ζ₇^g = (−1 + i√7)/2:
- **stage 1** (χ ≠ 1, ψ = 1): row sums;
- **dual** (χ = 1): per column j, w_j = x(R₀,C_j) + 3x(R₁,C_j) + 3x(R₂,C_j); 13 real equations (t and t + 13 are
  conjugate);
- **cross** (χ in the R₁ orbit): u_j = x(R₀,C_j) + η x(R₁,C_j) + η̄ x(R₂,C_j); 26 real equations;
- in each, the value at character ψ_t is the column-0 term plus Σ_j v_j P_{t,j}, with P the Gauss periods of index 26
  and P_{t,j+13} = conj(P_{t,j}).

**Column states.** A column (x₀, x₁, x₂) ∈ ±1³ is one of two types:
- *sym* (x₁ = x₂ = ε): w = x₀ + 6ε ∈ {±5, ±7}, u = x₀ − ε ∈ {0, ±2} (real);
- *anti* (x₁ = −x₂ = ε): w = x₀ ∈ {±1}, u = x₀ + iε√7.

Given w_j, the column is unique for w = ±5, ±7 and has two choices (the sign of ε) for w = ±1. So the cross level's
freedom given the dual is at most one bit per ±1-column, F ≤ 104 bits. The cross level has 26 equations
(C_cross ≈ 140 bits [E, rough; measure it]). The 933 order (dual first, then cross) therefore yields about
2^(F − C_cross) ≤ 2^−36 cores per dual solution at a 2^52 join each: **not viable** [E]. The +52 bits per y come from the
huge number of dual solutions (6^26 vectors per block before the equations).

**The half-degree subfield.** With s_j = w_j + w_{j+13} and d_j = w_j − w_{j+13} (j < 13), each dual equation reads
(Σ_j s_j Re P_{t,j} + column-0 term)² + (Σ_j d_j Im P_{t,j})² summed over the blocks. Check this identity numerically
first. So the dual splits into a "real" part (13 pair sums) and an "imaginary" part (13 differences), coupled only
through the sum of their two energies per equation and through the per-column letter constraint (w_j, w_{j+13}) ↔
(s_j, d_j). The cross level has the same real/imaginary structure in u.

## 3. What to do

1. **Reproduce the numbers.**
   - Stage 1: `python3 theory/ab_b0_stage1.py '<cfg>'` should give 38,380 y's (up to symmetry).
   - Stage-2 H₂/C₂ on a few y's: `theory/ab_gauss.stage2_sampled`, with the DP fibre sampler.
2. **Same-shape control first: v = 553 = 7·79**, s = 26 (f = 3, odd: the same C/−C pairing) with U = ⟨2⟩ on Z₇. Also
   v = 7·157 = 1099 (f = 6, even: the classes are self-paired) if useful. Enumerate stage 1 there, and decide small
   instances exactly (brute force or exhaustive joins over the 81-cell arrays at a reduced size) so that every
   Gaussian you use is calibrated against exact counts. The cross-level Gaussian at 933 was calibrated this way
   (`theory/ab_b0z3_calib.py`).
3. **Design a level order that is cheap here.** Candidates, to be costed with exact counts on the control before any
   production:
   - *Subfield-split dual:* list real parts (s) and imaginary parts (d) per block under the energy (spectral) bound, join
     them across blocks against the 13 equations as an 8-list problem, then the cross level.
   - *Type pattern first:* the sym/anti pattern of each block's 26 columns fixes |u_j| and the parity structure of w.
     It may let the cross level be solved before the dual, or jointly with it.
   - *Cross-first:* the anti columns carry all of the cross level's imaginary part.
   - *Spend slack on structure:* with +52 bits per y you can impose further constraints that cut the search, e.g. a
     prescribed sym/anti pattern, column-0 values, or a symmetry of some blocks, as long as the calibrated estimate
     stays positive. The natural sub-slices s = 13 and U = Z₇^* are dead (table above), so the restrictions must be finer
     than multiplier slices.
   - The same ideas apply to **955 = 5·191** (s = 38 = 2·19 with U = −1 on Z₅: 274 y's, +113 bits per y [E]; s = 19
     with U = 1: 3.27·10⁶ y's, +35.9 [E]). Treat it as the fallback target.
4. **Cost it** in core-hours per expected core with the calibrated figures; run it if within budget; verify, bundle,
   report.

## 4. Tools in the repo

- `theory/ab_sweep.py` (`build(cfg)`: orbits, letters, targets), `theory/ab_s1sd.py`, `theory/ab_b0_stage1.py`,
  `theory/ab_gauss.py` (Gaussians, fibre sampler);
- `theory/ab_b0z3.py` (the 933 solver: stage 1, dual, cross join), `theory/ab_b0z3_calib.py`,
  `theory/q523_gauss.py` (lattice-exact Gaussian);
- `scripts/verify_h3732.py` (the independent verifier with the slice check; adapt for 917),
  `scripts/make_bundle_b0.py`, `pub/H3732/` (bundle format);
- the report sections to read first: `report/alpoge_design_answer.md` §1, §4, §6 and §9.2; then
  `report/next_design_answer.md` §0 on when Gaussians fail.
