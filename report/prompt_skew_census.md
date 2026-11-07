# Prompt (high effort, bounded): new skew-Hadamard orders 4v from the Đoković–Seberry census, with this repository's staged tools

*Self-contained, for a fresh Claude Code session on `middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`.*
**Time box:** about 30 hours of wall-clock time on 4 cores. Commit after each step. The goal is **one or more
skew-Hadamard matrices of orders 4v not known before**, verified and published, or an honest map of why the cheap
candidates fail.

**Read first:**
1. this file;
2. `data/skew_open_census.json` (the open list and its provenance; §1);
3. in Kumar's repository `https://github.com/vatskishan03/hadamard` (Apache-2.0), at commit **d885022** (6 November
   2026) or later:
   - `skew/README.md` (every skew order he has, with structure and status);
   - `search/NOTES.md` §10, §12, §14, §20, §22, §23, §26, §35, §36, §48, §66, §73, §87–§89 and §110–§115 (his skew
     methods, calibrations, negative results and census scheduler);
   - `NOVELTY.md` (the literature he searched);
4. in this repository:
   - `report/open_screen_answer.md` §1–§4 (slice screens, the calibrated cost model, the stage-1 tool);
   - `report/h4892_answer.md` §0 and §5 (stage-1 and ℤ_q-level gates, the energy identity);
   - `report/h2396_plugin_answer.md` §0 (plug-in calibration);
   - `report/b0_883_answer.md` (the b = 0 flat-pool k-tree that built H(3532)).

## Rules

- **Kumar's repository:** clone it read-only. Call the `add_repo` tool for `vatskishan03/hadamard` with read access,
  then run `GIT_LFS_SKIP_SMUDGE=1 git clone --depth 50 https://github.com/vatskishan03/hadamard /home/user/vatskishan03/hadamard`.
  - Run `git log -1`. If it is newer than d885022, recompute §1 against it before anything else.
  - Never push to it. If you adapt his code (for example `search/cp_lowe.py`, `cp_enum.c`, `kt4.c`,
    `skew_certify.py`, `spin_ls2.c`), say so in the file header: Apache-2.0, credit K. Kumar.
- **Pushing:** before every push run
  `git fetch -q origin claude/focused-planck-j1645f && git merge -q --no-edit origin/claude/focused-planck-j1645f`,
  then `git push -u origin HEAD:claude/focused-planck-j1645f`.
- **Branches:** never push to `main` or to `codex/order-3452-compression`. Open no pull requests.
- **Commits:** trailers follow the existing commits. Scripts go in `theory/skew_*.py` (or `.c`), logs in `logs/skew/`,
  bundles in `pub/skew/<order>/`.
- **Processes:** never kill processes you did not start. Never `pkill -f` with a pattern matching your own shell;
  kill by PID. Never `rm` a path built from an unchecked shell variable.
- **Restart safety:** runners resume from their output files. Arm one `send_later` check-in while workers run.
- **Tags:** mark claims [P] proved, [M] exact computation (script named), [E] estimate. A cap or a timeout is never a
  zero.
- **Compute:** 4 cores and about 15 GB of RAM. Install `numpy`, `numba` and `ortools` with pip if they are missing;
  `gcc` is available.
- **Posting:** do not post anywhere outside the repository. The owner decides that.

## 1. What is open (computed 2026-11-07; recompute and correct if needed)

`data/skew_open_census.json`:
- **The census** (Đoković and Seberry, May 2025, as transcribed in Kumar's `search/data/skew_census_undecided.txt`)
  lists **737** odd v < 2500 with no skew-Hadamard matrix of order 4v.
- **Settled since:** 107 of them, all by Kumar except Schumer's 919 (order 3676, Zenodo, September 2026).
  - Kumar's orders of exactly 4v with v odd come from Goethals–Seidel difference families, Seberry's theorem with new
    ingredients, and the Paley-skew theorem of `skew/PALEY_SKEW.md`.
  - His orders 8v and 16v (negacyclic arrays, Wallis products, the 8v rows of the Paley-skew table) do **not** settle
    the census value v. They are excluded from the 107.
- **Open now: 630 values,** 227 of them prime and 82 at most 600. The smallest are 89, 101, 107, 149, 153, 167, 177,
  179, 191, 193, 201, 223, 225, 229, 233, 235, 257, 265, 269, 277, ….
- **Before building on this list, also check:**
  - EKLO v2 (arXiv:2610.27175v2; its Table 3 gives 356 = 4·89 as the smallest open skew order);
  - Schumer's Zenodo record;
  - a quick web or arXiv search for skew-Hadamard constructions dated 2025–2026, if the network allows.

  Record what you checked in the answer.

## 2. The construction (exact)

**The array.** A skew Goethals–Seidel array over an abelian group G of order v (v odd).
- Four ±1 functions A, B, C and D on G with Σ PAF = 0 off 0 (so Σ PAF(0) = 4v).
- A is of skew type: A(0) = 1 and A(−g) = −A(g) for g ≠ 0.
- Then the Goethals–Seidel array, with group-developed blocks and R the matrix of g ↦ −g, is a skew-Hadamard matrix
  of order 4v: H Hᵀ = 4v·I and H + Hᵀ = 2I (Goethals–Seidel 1970; Wallis).
- Row sums: A has row sum 1, so 4v = 1 + b² + c² + d² with b, c and d the other row sums.

**The slices.** As in Kumar's skew families and this repository's b = 0 screens:
- A, B, C and D are constant on the orbits of a multiplier group H ≤ Aut(G) of **odd** order, so −1 ∉ H and negation
  pairs the nonzero orbits.
- A is then determined by one sign per pair of orbits.
- Variants he used, which you may use: Paley products over ℤ_m × F_q with w-positions; cyclotomic products; spin
  structures with an extra μ ∉ H; non-cyclic G (ℤ₅² × ℤ₁₁, GF(343), and so on).

**Verification of any matrix:** a script sharing no code with the search must check integer PAF on every lag, the
skew first block, the H-invariance, H Hᵀ = 4v·I and H + Hᵀ = 2I. Copy the pattern of `scripts/verify_h2092.py`, and
check also with Kumar's `verify/verify.py` format if you write his certificate format.

## 3. What Kumar has already done, and where this repository's tools differ

**His coverage is large and active.** Between 1 and 5 November a census scheduler on a 32-core, 246 GB machine found
2300, 6628 and 7204. It ran heuristic searches on 661 and 883 to budget (30 h and 44 h) without a solution (§115).
His tools:
- low-energy pools with a ranked k-tree (`cp_lowe.py`, `cp_enum.c`, `kt4.c`; lists up to 2·10⁸);
- tabu search on spin structures (`spin_ls2.c`, hundreds of core-hours);
- an exact composite-census screen at 8–50 orbits (`census_screen.py`: 903 of 1643 structure-patterns survive, over
  131 values; §14);
- prime ranks with lists of 1.2·10⁷ (`skew_rank.py`).

**Do not race him on targets he is visibly running or has listed as leads.** These include the spin primes 1951, 2269
and 2437 (§14), the values named in §110–§115, and anything his notes mark as queued.

**This repository's tools that his notes do not use, or use differently:**
- **Exact staged decisions.** Compress along a large prime factor, decide stage 1 exactly (`theory/open_screen_stage1_exact.py`,
  method B, strong pruning), then the ℤ_q-level energy identity of `report/h4892_answer.md` §5. Then cost stage 2
  with the calibrated deficit model (`report/open_screen_answer.md` §2). This kills empty slices in minutes, before
  any long search. Kumar's skew runs rely on heuristic estimates, which he says are "blind to exact emptiness".
- **Composite v with 50–130 orbits.** His exact composite screen stopped at 50. Our staged pipeline handled 54–140
  (523, 933, 4892).
- **The b = 0 flat-pool k-tree** (`report/b0_883_answer.md`, `scripts/kt_*.py`) and the twisted slices (H(3676)), with
  the skew constraint added to block A.
- **The Theorem-B style use of exact normal forms** at coarse quotients (`report/openai_4p3_review.md` §1), where they
  apply.

## 4. Steps

**Step 1. The coverage map (≤ 3 h).**
- For every open v (§1), record from Kumar's notes whether and how he tried it (structure, orbits, method, budget,
  result), and whether he lists it as a lead.
- Write `logs/skew/coverage.jsonl`, and commit.

**Step 2. Our screen of the open values (≤ 8 h).**
- **Slices to list,** for each open v not marked as his lead, smallest v first: odd-order multiplier slices on every
  abelian G of order v (cyclic, and non-cyclic for square factors), with 24–130 orbits per block. Include his variant
  shapes where he did not cover them.
- **Cheap exact gates first:**
  - the skew row-sum condition 4v = 1 + b² + c² + d² with each sum in the orbit-weight residue classes;
  - lag-0 and row feasibility;
  - the ℤ₂ screen, where relevant;
  - coset-Parseval at each prime divisor;
  - Theorem α and Lemma R analogues, where they apply.
- **Then** the calibrated estimate, read at the true target (expected solutions, the stage-2 deficit for staged
  slices, and the cost per solution).
- **Rank** by estimated cost on 4 cores and 15 GB, excluding slices Kumar decided or exhausted.

**Step 3. Calibration (≤ 3 h), before any production.** Reproduce blind two or three of Kumar's skew orders with the
same tools and budgets as in Step 4, using only his structure, never his signs:
- **Suggested:** 476 (v = 119, H of order 3), 1004 (v = 251 prime, H of order 5), 1412 (v = 353, H of order 11) or
  3052 (v = 763 = 7·109, H of order 27).
- **Report:** time to first solution, and found against expected.
- **Kill rule:** if a chosen method cannot reproduce a control of comparable estimated cost, do not use it in
  production.

**Step 4. Exact stage 1, then search, for the top candidates (the rest of the box).**
- For each candidate in cost order, decide stage 1 exactly (two methods, planted controls) and apply the ℤ_q-level
  gates before running any long search.
- Then search with the cheapest calibrated method: an exact staged lift, the flat-pool k-tree, or Kumar's ranked
  k-tree on our hardware. Run at most 4 workers, each restart-safe, with a budget per candidate fixed in advance from
  the calibration (for example at most 3× the predicted time per solution, and never more than 12 hours).
- **On a solution:**
  - verify it (§2);
  - check novelty again (census, Kumar's current `skew/README.md`, EKLO v2, any 2025–2026 source found in §1);
  - commit at once;
  - write `pub/skew/<order>/` with the sequences, the group and H, a certificate, `verify.py`, and a README with the
    construction, search, verification and novelty status;
  - add a line to the main `README.md`.

## 5. Deliverable

Write `report/skew_census_answer.md` with:
- the bottom line (new orders found, verified, published; or none);
- the coverage map summary;
- the ranked screen table (top 30), with gates, estimates and the next step for each;
- the calibration results;
- Step 4 decisions and runs, with exact results;
- any matrices, with verification and novelty;
- the next computation.

Commit, push, and stop.

## 6. Calibration warnings (from this repository and Kumar's)

- Estimates are not evidence of solutions. In this repository, every positive estimate decided exactly at the ordinary
  open orders was zero, mostly at the coarsest quotient. Kumar's estimator is "blind to exact emptiness" (440 of 527
  empty small cases had finite estimates). Decide stage 1 exactly first.
- On prime v, his low-energy Gaussian was good to a factor of a few (353); on composite v it was off by orders of
  magnitude in both directions (763).
- **Concurrency is a real risk.** Kumar is working this census on far more hardware. Prefer values and structure
  classes his notes show as untried. Record the date and time of every hit.
