# Prompt (maximum effort, targeted): order 2396 (v = 599) by Kumar's Paley plug-in, plus two combinations of Kumar's methods with ours

*Self-contained, for a fresh Claude Code session on `middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`.*

This is a targeted reasoning task with go/no-go gates. It is not an open-ended search.
- **Part A** is the main question.
- **Parts B and C** are bounded checks of two places where a published block in Kumar's censuses might be pushed
  through with this repository's tools.

Stop each part at its first kill. Write the answer file, and stop. **Total compute:** about 36 hours of wall-clock
time on 4 cores. Commit after every gate.

**Read first:**
1. this file;
2. in Kumar's repository `https://github.com/vatskishan03/hadamard` (Apache-2.0):
   - `hadamard/2572/README.md`;
   - `search/NOTES.md` §18 (bordered arrays with a character set T), §78 and §80–§85 (two-level sampling and Paley
     plug-ins), §102 (Theorem 7 with larger index, order 3356), §103 (order 3820 and "alternatives that fail") and
     §104 (the remaining complex Yamada–Kiyasu targets);
3. in this repository:
   - `report/open_screen_answer.md` §0–§4 (the screen; the cost model; the stage-1 tool);
   - `report/h4892_answer.md` §0 and §5 (stage-1 and ℤ_q-level gates; the energy identity);
   - `report/paley_523_answer.md` §0 and §3 (a staged cyclotomic lift that worked);
   - `report/b3_push_answer.md` §7 (borders) and `report/openai_4p3_review.md` §2 (Theorem B, a universal border by
     normal form).

## Rules

- **Kumar's repository:** clone it read-only. Call the `add_repo` tool for `vatskishan03/hadamard` with read access,
  then run `GIT_LFS_SKIP_SMUDGE=1 git clone --depth 50 https://github.com/vatskishan03/hadamard /home/user/vatskishan03/hadamard`.
  Never push to it. If you adapt his code (`search/plug7_gen.c`, `search/plug7_l2.c`, `search/plug7_sol.py`,
  `search/tbcore/`, `search/bgschain/`), say so in file headers: Apache-2.0, credit K. Kumar.
- **Pushing:** before every push run
  `git fetch -q origin claude/focused-planck-j1645f && git merge -q --no-edit origin/claude/focused-planck-j1645f`,
  then `git push -u origin HEAD:claude/focused-planck-j1645f`.
- **Branches:** never push to `main` or to `codex/order-3452-compression`. Open no pull requests.
- **Commits:** trailers follow the existing commits. Scripts go in `theory/h2396_*.py` (or `.c`), logs in `logs/h2396/`.
- **Processes:** never kill processes you did not start. Never `pkill -f` with a pattern matching your own shell;
  kill by PID. Never `rm` a path built from an unchecked shell variable. Arm one `send_later` check-in while workers
  run.
- **Tags:** mark claims [P] proved, [M] exact computation (script named), [E] estimate. A cap or a timeout is never a
  zero.
- **Matrices:** a matrix counts only if H Hᵀ = n·I is verified exactly by a script sharing no code with the search, as
  in `scripts/verify_h2092.py`. Every core is also rechecked by integer PAF on all lags.
- **Compute:** 4 cores and about 15 GB of RAM. Install `numpy`, `numba` and `ortools` with pip if they are missing;
  `gcc` is available.

## Part A. 599 by a Paley plug-in Legendre quadruple (main)

### A.1 The idea, and why it is worth one session

- **Kumar's 643.** Kumar settled 643 (order 2572) with a *Paley plug-in* core on G = ℤ₂₃ × F₇: each block is
  x₀ ⊗ δ + x₁ ⊗ χ, with χ the Legendre symbol of F₇, except at a few *defect* positions where x₁ = x₂. An exact identity
  forces 12–22 defects.
- **Why plug-ins are special.** In his calibration (§81) a small plug-in class (ℤ₃² × ℤ₁₇) had tens of thousands of
  solutions where the Gaussian predicted 5, and every one had this defect structure. Solutions concentrate near the
  algebraic Paley point.
  - That is the opposite of this repository's repeated "positive estimate, exact zero" at generic slices.
  - A generic two-level sampler on 643's group gave about 10⁻⁷ completions per instance. The plug-in sampler found the
    core in about 13.5 CPU-hours.
- **The gap.** He used plug-ins only inside quaternary (Szöllősi) cores, for the complex Yamada–Kiyasu route. His
  Legendre-quadruple ("Spence") censuses covered only n ≡ 1 (mod 4); 599 ≡ 3 (mod 4).
- **The host.** The plug-in host for real Legendre quadruples is the b = 1 slice with no multiplier on K and the
  squares on F_q. This repository's screen tested 57 such slices at open orders (`logs/open_screen/screens_s*.jsonl`,
  `logs/open_screen/z2gate.jsonl`):
  - 50 die at screen A, 6 more at the ℤ₂ screen or the w-screen;
  - **exactly one passes every exact screen: v = 599, key `["b1", 599, [], [23, 2], 13, 2, []]`.**
- **That host:**
  - G = ℤ₄₆ × F₁₃ (ℤ₅₉₈ by the CRT; ℤ₄₆ = ℤ₂₃ × ℤ₂), q = 13, |K| = 46, r = 138 orbits per block;
  - ℤ₂ screen passed with 3 compositions, w-screen passed;
  - generic estimate about 2^107 solutions in the slice, but a staged cost of about 2^72 (`logs/open_screen/deficit_s*.jsonl`).
  - Its |K|/q ≈ 3.5 is close to 643's 23/7 ≈ 3.3.
- **The risk** (Kumar §81, "existence trend"): "the plain plug-in runs out of freedom". At q = 17 he found solutions
  for |K| = 9, 11 and 13 and none for 15, 17 and 19. Whether the defect version still has solutions at |K| = 46 is the
  question this part decides first.

**The target.** Theorem 5 (EKLO): four ±1 sequences on G of even order m = 598, with Σ PAF = −4 off 0 and row sums
(±2, 0, 0, 0), give H(4m + 4) = **H(2396)**, the smallest open order. Assembly and verification as for H(2092)
(`pub/H2092/`, `theory/q523_finish.py`, `scripts/verify_h2092.py`).

### A.2 Gates

**A1. Derive the plug-in system for Legendre quadruples [P] (≤ 3 h).**
- **Set-up.** G = K × F_q, q ≡ 1 (mod 4) (q = 13, so χ is even and the Paley matrix Q is symmetric), H = {1} × S with
  S the squares of F_q^*. Each block is x₀ on K × {0}, x₁ on K × S, x₂ on K × N.
- **Variables:** s = x₀ − e, t = (x₁ − x₂)/2, e = (x₁ + x₂)/2, so X = sI + tQ + eJ.
- **Derive the K-level system** for the target Σ PAF = (4|G| + 4)δ − 4, the analogues of Kumar's (I), (J) and (Q) in
  §81.
- **Derive:**
  - the **defect identity**: the counts of defects of each type forced by the trivial-character equation;
  - row-sum lemmas;
  - every exact filter you can find, including mod-small congruences;
  - the ℤ_q-level energy identity of `report/h4892_answer.md` §5, which killed every non-empty stage 1 there.
- **Check** the derivation against brute force on small K, both for solutions and non-solutions.
- **Kill** if the identity has no solution at K = ℤ₄₆, q = 13.

**A2. Calibrate the existence trend [M] (≤ 8 h). This is the go/no-go.**
- **Analogues:** K × F₁₃ with |K| = 2k, k odd, and v − 1 = 13|K| (|K| = 6, 10, 14, 18, 22, 26, …): v = 79, 131, 183,
  235, 287, 339, …, all known orders. Also |K| ≡ 0 (mod 4) if the defect identity allows it.
- **Count exactly,** for each |K| that fits: the plug-in solutions (all defect allocations), the Gaussian (or tilted
  Gaussian, as in §82) prediction, and the defect counts.
- **Where exact counts stop being affordable,** run the two-level sampler (A3) at a fixed budget and record completions
  per instance, as Kumar did at p = 13 and 17.
- **Fit** log(count or rate) against |K| and extrapolate to 46, with an honest uncertainty.
- **Kill** if the solutions vanish at an |K| well below 46, or the extrapolated cost per core at |K| = 46 exceeds about
  2^46 operations (four cores for a few days), with at least one bit of margin.

**A3. The sampler.**
- **Design.** Adapt Kumar's two-level design (§82–§83: `plug7_gen.c` level 1 and `plug7_l2.c` level 2), or this
  repository's staged join, to the real Legendre-quadruple plug-in on K × F₁₃.
  - Level 1 samples x₀ and the defect allocation, with the exact filters.
  - Level 2 solves the t signs by a meet-in-the-middle join on the remaining K-level equations.
- **Correctness first.** Validate on the calibration analogues: brute-force equality on complete sets where possible,
  and planted recovery at the target size.
- **Measure** the time per instance at |K| = 46.

**A4. Production (only if A2 passes).** Run with restart-safe workers, at most the remaining time box. Log instances,
completions and the measured rate against the prediction.

On a core:
- verify it;
- assemble H(2396) and verify it independently;
- commit at once, and write `pub/H2396/` in the format of `pub/H2092/`.

**Do not post anywhere outside the repository;** the owner decides that.

## Part B. 599 over ℤ₄² × ℤ₃₇: push through Kumar's border block (≤ 8 h)

- **The block.** In `search/NOTES.md` §103 ("Bordered arrays over ℤ₄ʳ with partial Teichmüller character sets") Kumar
  records: *599 over ℤ₄² × ℤ₃₇ (t = 7, 70 orbits, 10^24.6): no border in any symmetric class tested.*
  - This is the family of §18, with Q = ℤ₄² and T = the trivial character plus orbits of order-4 characters,
    |T| = t = 7, |G| = 592 = 599 − 7, border of size 4t = 28.
  - 70 orbits with an estimated 10^24.6 solutions is the richest moderate-size structure anywhere in his notes.
  - Its 3356-style chain (exact quotients, heavy-vector last level) was not run, because the border failed first.
  - He did not search for borders without symmetry; at 478 he notes "borders without symmetry were not searched".
- **B1. Reconstruct the structure exactly** from §18, §102 and §103:
  - G, K, Q, T and H;
  - the target −Γ(g + K) with Γ = 4 Σ_{χ∈T} χ;
  - the Q-level compressions and their consequences: Σ|ŷ(χ)|² = 4t on T and 4n off T.

  Recompute his 10^24.6 with this repository's estimator, as a consistency check. If the structure cannot be pinned
  down uniquely from the notes, state the alternatives and treat each.
- **B2. The border problem, without symmetry.**
  - **What to find:** X (Hadamard of order 28), R̂ (28 × 64, orthogonal rows) and Ĉ (64 × 28, Ĉ Ĉᵀ = I₄ ⊗ Γ), with
    X Ĉᵀ + R̂ Mcᵀ = 0, for each feasible Q-level compression class.
  - **First:** enumerate the Q-level compression classes exactly. Ask whether a normal form like the one in
    `report/openai_4p3_review.md` §1 reduces them to a few classes, and whether a universal certificate argument like
    Theorem B applies.
  - **Then decide existence per class:** CP-SAT on the ±1 border variables with exact constraints, frame search from
    `theory/bt_border.py`, or a proof of nonexistence (for example from the support condition: rows of R̂ in the span of
    the T-characters).
  - **Kill** if no class has a border [P or M]. Then Part B is closed and the block is real.
- **B3. Only if a border exists:** run the cheapest exact quotient gates on the core.
  - The ℤ₃₇ level first, then the ℤ₄² levels.
  - Use Kumar's rule: test the Q-level join, and the quotient equations modulo m where H has a kernel of order m ≥ 9.
  - Estimate the 3356-style chain cost (heavy-vector last level, mod-8 span test, exact join).
  - **Go** to production only if a non-empty finer quotient exists and the estimated cost is ≤ 2^44. Otherwise record
    it and stop.

## Part C. The index-4 Szöllősi plug-ins at 1291 and 2083: is Kumar's "search" block our kind of problem? (≤ 6 h)

- **The block.** Kumar (§82, "Index-r plug-ins") estimates by a tilted Gaussian about 10^8.8 solutions for the complex
  Yamada–Kiyasu core of n = 1291 (order 5164), with K = ℤ₁₉, q = 17, r = 4, and 10^14.8 for the bordered n = 2083
  (order 8332), with K = ℤ₄₀, q = 13, r = 4. He writes "the obstruction is search, not existence":
  - the blocks have 95 bits;
  - the s/t decoupling of r = 2 has no analogue;
  - annealing stalls;
  - exact joins reach only |K| ≤ 5.
- **Why it might be ours.** An index-r plug-in is constant on the index-r cyclotomic classes of F_q in every K-fibre.
  That is exactly a generalised-Paley slice: compress along F_q to K (stage 1 with letters x₀ + L·σ), then a cyclotomic
  stage 2 with fixed class sums. This repository's staged solvers (523, 4892) do this.
- **C1. Formulate** the 1291 core (zero-diagonal Szöllősi array: A and B symmetric, {0, ±1}, vanishing at 0; C and D
  ±1; Σ PAF = 0 off 0 on N = 323) as a staged slice. Kumar's verifier (`verify/verify.py`, `_szollosi_zero_diagonal_core`)
  and his 955/1051/1699 certificates give the conventions and positive controls.
- **C2. Decide stage 1 exactly on K = ℤ₁₉** (A, B symmetric), and apply the ℤ_q-level energy identity to every stage-1
  solution, as in `report/h4892_answer.md` §5.
  - **Kill** if stage 1 is empty, or every stage-1 solution fails the ℤ_q level. This is likely, given 4892 and 7572.
- **C3. If alive:** cost stage 2 with the 523/4892 model at real stage-1 solutions. **Go** only if the cost is ≤ 2^44.
  Then search stage 2 exactly for every stage-1 solution, as for 4892.
  - A core gives a quaternary conference matrix of order 646, then QH(2582) by complex Yamada–Kiyasu
    (q = 1289 ≡ 1 mod 8), then H(5164) by Turyn doubling. Use Kumar's builder; verify independently.
- **2083:** run C1–C2 only if 1291 is decided early.

**Set aside after review** (record in the answer, one line each; no compute):
- 1187 (k4 over ℤ₄² × ℤ₂ × ℤ₃₇, work 10^26);
- 1987;
- 1223 over GR(4, 3) × ℤ₁₉ (heavy-leaf fibres of 10^7.2–7.6 per block, yield 10^−4.1 per leaf: about 2^63 per
  solution by his figures; the S/D splits measured here do not bridge that);
- 719 at index 16 (work 10^21.8);
- 1879 and 2215 (t = 1, about 0.03 inequivalent solutions).

## Deliverable

Write `report/h2396_plugin_answer.md` with:
- the bottom line, one line per part with its verdict and the gate reached;
- per gate: methods, scripts, logs, planted and brute-force controls, and results with tags;
- for Part A, the calibration table (|K|, exact count or rate, prediction) and the extrapolation;
- any matrix, with its verification;
- the next computation.

Commit, push, and stop.

## Calibration warnings (from this repository and Kumar's)

- Estimates are not evidence of solutions. At the open orders every exactly decided positive estimate in this
  repository has been zero, and Kumar's estimator "is blind to exact emptiness" (440 of 527 empty small cases had
  finite estimates).
- Plug-in classes are the documented exception: there the Gaussian undercounts. Measure; don't assume either way.
- Exact quotient tests, energy identities and congruences decide; run them before any long search.
- A planted test checks correctness, not cost.
