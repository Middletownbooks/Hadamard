# Prompt (maximum effort): bordered Goethals–Seidel arrays with a character set T, the next push at the open orders below 4000

*Self-contained, for a fresh maximum-reasoning Claude Code session on `middletownbooks/hadamard`, branch
`claude/focused-planck-j1645f`.*

**Read first, in this order:**
1. this file;
2. `README.md` (status of our orders);
3. `report/b599x4_answer.md` §0–3 (a general obstruction, and why at-mean estimates fail);
4. `report/bge3_stage1_answer.md` §0–2 and §8 (our b ≥ 3 coverage and its gaps);
5. `report/alpoge_design_answer.md` §1, §3 and §5 (our staged b ≥ 3 systems, borders, the 347 control).

**Then clone Kumar's repository** (Apache-2.0; read-only here). Call the `add_repo` tool for `vatskishan03/hadamard`
with read access, then run
`GIT_LFS_SKIP_SMUDGE=1 git clone --depth 50 https://github.com/vatskishan03/hadamard /home/user/vatskishan03/hadamard`.
Read in it:
- `search/NOTES.md` §18–19 (the general construction and the ℤ₂ split that built 2148);
- `hadamard/3964/README.md`, `hadamard/2636/README.md` and `hadamard/2148/README.md`;
- `METHOD.md`, `FRONTIER.md` and `search/README.md` (the "Bordered arrays" and compression-chain tools).

## Rules

- **Pushing:** other sessions push to this branch. Before every push run
  `git fetch -q origin claude/focused-planck-j1645f && git merge -q --no-edit origin/claude/focused-planck-j1645f`,
  then `git push -u origin HEAD:claude/focused-planck-j1645f`.
- **No `main` and no pull requests:** never push to `main` (the owner fast-forwards it) and open no pull requests.
- **Commits:** trailers follow the existing commits.
- **Processes:** never kill processes you did not start. Never `pkill -f` with a pattern that matches your own shell;
  kill by PID.
- **Restart safety:** write runners that resume from disk. Arm a check-in (`send_later`) whenever workers run; an
  idle session was reclaimed once and killed production workers. Commit results as they arrive.
- **Check before writing:** look at `git log` and the existing files first; several sessions have rebuilt tools that
  already existed.
- **Tags:** [P] proofs, [M] exact computations (name the script), [E] estimates (name the estimator, its calibration,
  and whether it is read at the true target).
- **Verification:** a construction counts only if verified exactly by a script that shares no code with the search.
- **Credit:** reusing Kumar's code or ideas is fine (Apache-2.0); credit it in code headers and in the report.
- **Compute:** 4 cores, about 15 GB of RAM, about 30 GB of scratch disk. Install `numba`, `ortools` and `sympy` with
  pip if missing; a C compiler with OpenMP is available.

## 0. Goal and deliverables

A Hadamard matrix at one of the open orders below 4000. By Table 5 of Egan–Kotsireas–London–Ó Catháin (2026, cited
in `README.md` as [EKLO]) these are 4n for

**n = 599, 643, 669, 719, 739, 751, 839, 859, 863, 907, 947, 955, 971.**

The search uses the bordered Goethals–Seidel construction with a character set (§1). All recent b ≥ 3 successes are
instances of it:
- Alpöge et al.: 1388, 1436, 1916;
- Kumar: 2148, 2636, 3964;
- [EKLO]: 5788, 8476.

Deliverables, in order of value:
- **(A)** A verified H(4n) at an open n: an independent verifier, and a bundle `pub/H<4n>/` (sequences on G, the border
  and corner written out, the normalized matrix, a standalone `verify.py`).
- **(B)** A ranked census of the structures in §3 at every open n. For each: exact gate status [M], a stratified count
  [E] read at the true target, the border-problem status, and a cost. The best ones go into a calibrated production
  run.
- **(C)** Exact obstructions closing families.

Write to `report/b3_push_answer.md`, with scripts in `theory/bt_*.py` and logs in `logs/bt/`.

## 1. The construction

Kumar `search/NOTES.md` §18; [EKLO] Theorem 7 is the case where T is every nontrivial character.

- **Data:**
  - G abelian; K ≤ G of index q; Q = G/K;
  - T a set of t characters of Q, closed under Galois action;
  - Γ = 4 Σ_{χ∈T} χ, a function on Q.
- **Sequences:** four ±1 sequences on G with

      Σᵢ PAFᵢ(g) = −Γ(g + K)   for g ≠ 0.

- **Completion:** a border of 4t rows and columns, constant on the cosets of K, completes the Goethals–Seidel array
  over G to a Hadamard matrix of order 4n with n = |G| + t. This needs the **border problem** to be solvable:
  - X is Hadamard of order 4t;
  - R̂ is 4t × 4q with orthogonal rows;
  - Ĉ Ĉᵀ = I₄ ⊗ Γ;
  - X Ĉᵀ + R̂ Mcᵀ = 0, where Mc is the Goethals–Seidel array over Q of the compressions to Q.
- **Consequences:** Σᵢ |ŷᵢ(χ)|² = 4t for χ ∈ T and 4n otherwise, so the compressions must solve a small exact
  problem.

**The families that matter below 4000:**

| family | data | which open n | instances |
|---|---|---|---|
| Alpöge / [EKLO] Theorem 7, k = 4 | Q = ℤ₂², T = the three nontrivial characters, t = 3, \|G\| = n − 3 | every n ≡ 3 (mod 4): all open n except 669 | 1916, 3964 (G = ℤ₂² × odd), 2636 (G = ℤ₄ × ℤ₄ × ℤ₄₁, K = 2ℤ₄ × 2ℤ₄ × ℤ₄₁, non-split) |
| Theorem 7, k = 6 or 8 | t = 5 or 7 | requires k \| n − t: k = 6 fits 599, 719, 839, 863, 947, 971; k = 8 fits 599, 719, 751, 839, 863 and others | 1388, 1436 |
| ℤ₃, all characters | Q = ℤ₃, T = all three, t = 3: targets −12 on K \ 0 and 0 off K; Σ r² = 12 ("two spikes") | 3 \| n: **669** among the open orders | 2148 |
| product family | Q = ℤ₂² × ℤ₃, t = 9 | n ≡ 9 (mod 12): 669 | Kumar's estimates are below 10^−10 for 669 (§18); recheck exactly |

**Factorisations of n − 3** for the k = 4 family, which determine the groups G with K of index 4:

| n | n − 3 |
|---|---|
| 599 | 2²·149 |
| 643 | 2⁷·5 |
| 719 | 2²·179 |
| 739 | 2⁵·23 |
| 751 | 2²·11·17 |
| 839 | 2²·11·19 |
| 859 | 2³·107 |
| 863 | 2²·5·43 |
| 907 | 2³·113 |
| 947 | 2⁴·59 |
| 955 | 2³·7·17 |
| 971 | 2³·11² |

The 2-part can be organised many ways: ℤ₂² × ℤ_m, ℤ₄ × ℤ_m with Q a quotient, ℤ₄ × ℤ₄ × …, ℤ₂³ × …, ℤ₈ × ….

**The symmetry groups S that have worked** (all are subgroups of Aut(G) that preserve K):
- **3964:** S has order 27, generated by σ = (α on ℤ₂², ×66 on ℤ₂₄₇) and τ = ×35 on ℤ₂₄₇, where α is the order-3
  automorphism of ℤ₂². 46 orbits, 23 equations.
- **1916:** a split structure, with different groups for different blocks: X, W invariant under ⟨(α, 16), (1, 13)⟩
  (order 12, 46 orbits); Y, Z under ⟨(1, 13)⟩ (order 4, 128 orbits).
- **2636:** S has order 15: β of order 3 on ℤ₄ × ℤ₄ acting trivially on ℤ₄₁, and ×10 (order 5) on ℤ₄₁. 54 orbits
  (sizes 1, 3 five times, 5 eight times, 15 forty times), 27 equations.
  - Over ℤ₂² the same shape is empty; over ℤ₄ × ℤ₄ it is not, because β-invariant bent patterns meet the energy budget.
- **2148:** S = the eighth powers of GF(89) acting on ℤ₆ × GF(89) (order 11, 54 free signs per sequence). Solved by
  the ℤ₂ split into U and W parts: 284,292 U-solutions gave 1 completion, at the predicted rate.

**What our earlier b ≥ 3 sweeps did not cover** (`report/bge3_stage1_answer.md` §8.3):
- automorphisms of Q coupled with multipliers on the odd part (3964, 1916);
- non-split K (2636);
- the ℤ₃ family (2148);
- split structures with different groups per block (1916).

## 2. What is known and dead

**Exact obstructions to run before anything else:**
- **Coset-Parseval [P, M]** (`theory/b599x4_cosetscreen.py`, `report/b599x4_answer.md` §1–2). Let N ≤ K be a subgroup
  of prime order, and C_N the union of the S-orbits that are not unions of N-cosets. The structure has no core when
  |C_N|(1 + 1/|N|) < n (|N| odd), or |C_N| < (1 − 1/|N|)·n (|N| even).
  - It killed 599 ×4 and 259 configurations in all.
  - It bites whenever S contains an element acting on one axis ℤ_{p^e} (e ≥ 2) as u ≡ 1 (mod p^{e−1}), u ≠ 1, and
    trivially on the other axes.
  - Adapt it to the T-targets of §1 (the Parseval budget per character class changes with Γ).
- **Compression gates.** For every S-invariant subgroup N of G, the sums over N-cosets satisfy the same kind of
  equation on G/N (with Γ pushed forward). Lists of compressed vectors per block are joined 2+2.
  - Kumar `chain.py`, `cchain4.py`, `comp_scan.py`, `gcomp.py`; ours `theory/s1z_*`, `theory/ab_s1sd.py`.
  - The 3964 search was this chain to the end: the ℤ₁₉-level join gave 304 solutions, then the ℤ₁₃ level, with a DP
    listing only the 70–140 distinct compressed vectors per fibre.
- **The border problem.** Solvable for every compression pattern at Q = ℤ₂² (2048 of 2048 in our tests,
  `report/alpoge_design_answer.md` §3) and at Q = ℤ₃ (864 of 864, Kumar `border3.c`). Frames exist for every even k
  (`report/alpoge_design_answer.md` §3.1). For k ≥ 10 no flat corner is known. Check border solvability per surviving
  compression class before any deep search.

**Results of a Sonnet session** (the owner's sandbox; its scripts may appear later under `sonnet/` in this repo: use
them if present, re-derive them otherwise).

Notation there: b is the border width divided by 4 (= t).

Exact negatives:
- 955, b = 1, ℤ₁₈ × ℤ₅₃ (f = 13): empty. Only two sets of non-±1 entries are possible, and CP-SAT refuted all 2,059
  placements. This agrees with our `report/slice_screens_20261013.md`.
- 947, b = 5, ℤ₆ × ℤ₁₅₇: empty.
- 863 in ℤ₂² × ℤ₂₁₅ with an odd, α-paired symmetry: no solutions at the ℤ₂² gate.
- Cyclic ℤ_{n−3} with a ℤ₄ quotient and odd K: dead for 599, 643, 719, 751, 859, 907, 947, 955. The runs for 739, 839,
  863, 971 died: unknown.

Census of ℤ₂² × ℤ_κ (all symmetry orders, one common group; Gaussian E, not read at the true target):

| n | orbits | log₁₀ E | note |
|---|---|---|---|
| 839 | 78 | 11.7 | generator (α, 49) |
| 839 | 110 | 30 | |
| 863 | 84 | 6.0 | ⟨(1, 4)⟩, order 14 |
| 863 | 140 | 40 | |
| 955 | 92–100 | 1.5–7.6 | generators (α, 47), and (α, 67) with (1, 183) |
| 751 | 132 | 8.5 | |

- All other hosts were dead.
- ℤ₄² × ℤ_κ: nothing alive for 739 or 947; 643 was killed for memory.
- 971 in ℤ₂² × ℤ₂ × ℤ₁₁² with GL(2, 11): best E ≈ 10^10.7 at 112 orbits.

Controls that reproduced:
- 2636: 54 orbits, 27 classes, exactly 2,954 compressed vectors per sequence.
- 1916: 16 level-1 solutions at row sums (6, 6, 38, 20).
- 3964: 304 solutions at the ℤ₁₉ level.
- 5788 (ℤ₂² × F₃₆₁): E ≈ 10^19.9 at 86 orbits.

Not finished:
- second-level gates for 839, 863, 955 (839's first level: 6,720 level-1 quadruples, 20 row-sum patterns);
- a mixed-group estimator for split structures;
- border-completion checks.

**Calibration lessons [M]:**
- **Gaussian estimates fail in both directions here.**
  - 3964: about 20 predicted, and the structure was rich.
  - 1916: 7 predicted, more than 7,800 found in one branch.
  - Slices at +46 to +102 bits were exactly empty more than fifteen times in this repo.
- **Estimates read at the sample mean are not costs.** `report/b599x4_answer.md` §0.2:
  - on a benign slice they understate D_def by about 2 bits;
  - on a rigid one they are unbounded, because the target can lie outside the reachable set.
- **Rank only by:** exact gates, then stratified counts (exact for the coarse compressions, Gaussian for the rest),
  read at the true target, calibrated on controls with exact counts.

## 3. Plan

1. **Controls first** (each must be reproduced exactly before any open-order claim):
   - 3964: the ℤ₁₉ join gives 304 solutions; the ℤ₁₃ level then gives a core.
   - 1916: the split structure, 16 level-1 solutions at (6, 6, 38, 20).
   - 2636: the structure's counts (2,954 compressed vectors per sequence).
   - Optionally 2148 via the ℤ₂ split, or its n = 177 analogue: 26,095,424 solutions exactly in its stage-1 class.
2. **One structure engine** for the construction of §1. It takes:
   - G as a product of cyclic groups (any 2-part), K, Q and T;
   - S generated by automorphisms (multipliers per axis, automorphisms of the 2-part such as α and β, coupled
     generators), optionally different S per block;
   - and produces orbits, class equations with the Γ-targets, the coset-Parseval screen, and the list of S-invariant
     subgroups N for compression gates.
   Validate it on the controls.
3. **Census at every open n below 4000:** families k = 4 (all n ≡ 3 mod 4), k = 6 and 8 where k | n − t, and the ℤ₃
   and product families for 669.
   - Run the cheap exact gates first (coset-Parseval, the Q-level, the first compression level), then the stratified
     count at the true target.
   - Include split structures and non-split K.
   - Start from Sonnet's live hosts: 839 (78 orbits), 863 (84), 955 (92–100) and 751 (132), plus the unknown cyclic
     cases 739, 839, 863 and 971.
4. **Second-level and deeper gates** on the survivors, as compiled code (C or numba): DP per fibre over orbits to list
   distinct compressed vectors, exact 2+2 joins, pruning of row-sum patterns level by level.
   - Run them to the end where feasible (that is how 3964 was found).
   - Otherwise switch to the staged solvers: our multi-level stage-2 solver (`theory/m2_*`, validated at 347 and 523,
     `report/b599x4_answer.md`), or the ℤ₂ U/W split (Kumar `bgs_uw.py`).
5. **Border check** for each surviving compression class: Kumar `border3.c` and `border_probe.py`; ours
   `theory/ab_borderfit.py`, `ab_assemble.py`, `ab_template.py`.
6. **Production** on the best calibrated structure, restart-safe, with a tally of expected cores; do not stop early.
   - Any core: assemble the border, check H Hᵀ = 4n·I independently, commit at once, then bundle.
7. **If nothing reaches:** obstructions per family, and the cheapest remaining structure per order with its measured
   cost.

## 4. Our tools (read before writing new code)

| tool | purpose |
|---|---|
| `theory/b599x4_parseval.py`, `theory/b599x4_cosetscreen.py` | the coset-Parseval obstruction and its screen |
| `theory/s1z_lib.py`, `s1z_screen.py`, `s1z_cpsat.py`, `s1z_border.py`, `s1z_qzq.py`, `s1z_typescreen.py`, `s1z_zqscreen.py` | compression multisets, exact compression screens of stage 2, border classes |
| `theory/ab_lib.py`, `ab_sweep.py`, `ab_s1sd.py`, `ab_qaut.py`, `ab_assemble.py`, `ab_template.py`, `ab_borderfit.py`, `ab_border_c4.py` | staged b ≥ 3 systems, Q-automorphism slices, assembly and borders |
| `theory/m2_*`, `theory/ab_prod347.py`, `ab_prod_generic.py` | multi-level stage-2 solver and production (347: 23 cores at 2^15.7 S-solutions per core, model 2^16.1) |
| `theory/nd_joint_exhaust.py`, `nd_cross.py`, `nd_ejoin.py` | doubly staged joint tests, energy-split joins (647 ×353) |
| `scripts/verify_packed_h.py`, `scripts/verify_b0_bundle_rebuild.py`, `pub/*/verify.py` | verification; the b ≥ 3 bundle writer is from commit d7c7c6a |
