# Prompt (OpenAI model; ideas and proofs, short sessions): beat the deficit floor of staged Hadamard-core searches

*Self-contained, for a session with access to the public repository `https://github.com/middletownbooks/hadamard`.*
- **Push target:** you may push only to the branch `codex/order-3452-compression`.
- **The task:** come up with new algorithms or algebraic ideas, prove what you can, and estimate the rest. Each idea
  must leave behind a precise proposal with a cheap test.
- **Testing:** long searches are not your job. Claude Code instances with exact-gate tooling will test your proposals
  afterwards. Write them so that such an instance can implement each test without asking you anything.

## 0. Working rules (sessions may time out: checkpoint early and often)

- **Setup:**
  - `git fetch origin codex/order-3452-compression claude/focused-planck-j1645f`.
  - Check out `codex/order-3452-compression`.
  - Merge `origin/claude/focused-planck-j1645f` into it. The two branches share the Codex files unchanged, so the
    merge should be clean.
  - If it is not clean, abort the merge and read the Claude branch with `git show origin/claude/focused-planck-j1645f:PATH`.
- **Pushing:** push only to `codex/order-3452-compression`. Never push to `main` or to `claude/focused-planck-j1645f`.
  Open no pull requests.
- **Checkpoints:**
  - In your first ten minutes, create `report/openai_deficit_answer.md` with the section skeleton of §5, and push it.
  - Push again after every finished proposal, or every 30 minutes, whichever comes first.
  - If you run out of time, the file must say what is finished and what is not.
- **Runs:** keep any computation under about 10 minutes per run. Record the command, seed and output under
  `logs/openai_deficit/`.
- **Tags:** mark every claim [P] proved (proof in the text), [M] exact computation (script named) or [E] estimate
  (method named). Give an estimate with its assumptions. A "should work" with no cost is not a proposal.
- **Matrices:** never claim a Hadamard matrix without an exact integer verification of H Hᵀ = N·I.
- **Credit:** the bordered construction is Kumar's form of Theorem 7 of Egan–Kotsireas–London–Ó Catháin [EKLO]. The
  863 set-up and its first certificate are by OpenAI Codex. The 523 solver and the 863 lift analysis are by Claude Code.

## 1. The problem in one paragraph

**Hadamard matrices from four sequences.** Four ±1 sequences x₁, …, x₄ on a finite abelian group G whose periodic
autocorrelations sum to a prescribed function F give, through a Goethals–Seidel array plus a small border, a Hadamard
matrix of order 4n.

**Staged searches.** For the open orders we restrict the xᵢ to be constant on the orbits of a multiplier group S
(a slice), then solve in stages along a chain of quotients G → G/N₁ → … : first a compressed problem, then lifts.
At every open order we have, the last stage is the problem:
- the expected number of cores per stage-1 solution is tiny or huge, but the cost of finding one is governed by the
  **deficit** (§3);
- with every known split, the deficit is about 10–25 bits too large for four CPU cores.

**What we want:** an algorithm, or an algebraic construction, whose cost per core beats the deficit floor of §3. Or a
proof that a whole class of such algorithms cannot.

## 2. The concrete instance: v = 863, order 3452 (open)

Read `report/h863_lift_answer.md` §0, §1 and §3 first; the facts below are from there.

**The group, the target and the slice.**
- G = ℤ₂² × ℤ₅ × ℤ₄₃ (|G| = 860), K = ℤ₅ × ℤ₄₃, Q = G/K = ℤ₂².
- Target: Σᵢ PAF_{xᵢ} = 3440 at 0, −12 on K \ 0 and +4 off K. Equivalently, Σᵢ |x̂ᵢ(χ)|² = 3452 for every character
  χ nontrivial on K, and 12 for the three nontrivial characters of Q.
- A border of width 12 completing any such core is in `logs/h3452/border.json`.
- **Slice S7:** S = ⟨×4 on ℤ₄₃⟩, of order 7, acting trivially on H = ℤ₂² × ℤ₅.
  - Per block there are 20 cells h ∈ H. Each holds a singleton sign e(h) at (h, 0) and six class signs c_j(h) on
    h × C_j, where C_j are the six cosets of ⟨4⟩ in ℤ₄₃^*.
  - −1 ∉ ⟨4⟩; the classes pair under negation as (1, 2), (3, 6), (7, 9) by representatives.
  - There are 72 character classes, so 71 real equations.

**Stage 1 is solved [M].**
- The ℤ₄₃ compression y_i(h) = e_i(h) + 7 Σⱼ c_ij(h) is solved in seconds (`data/h3452_stage1.json`).
- Its fibre has 2^339.1 points and carries **2^63.5 expected cores** [E]: 59 remaining equations, 275.6 bits.

**Stage 2a is solved [M].**
- The G/ℤ₅ level is w_i(q, j) = Σ_u c_ij(q, u): the column sums of 16 tables of size 5 × 6 (4 blocks × 4 cells of
  ℤ₂²), whose row sums are fixed by y.
- It has 2^99 solutions per y [E]. 39 are verified in `data/h863_g5_solutions.json`.
- Producing one costs 2^9 pair visits (0.22 ms).

**Stage 2b is the bottleneck.**
- Given (y, w), the remaining unknowns are 16 sign tables of size 5 × 6 with **both margins fixed**, 2^172.7–2^175 of
  them.
- They must satisfy **47 real equations**, the characters nontrivial on both ℤ₅ and ℤ₄₃, at a cost of 209.8–210.3
  bits [E, calibrated to about 1 bit against exact hit counts].
- So each (y, w) holds about **2^−35.9 cores**, and a core needs about 2^35.9 pairs (y, w).
- To fit the budget, one pair (y, w) would have to be decided in about 2^10 operations. The best method found costs
  2^80.

**Everything already measured (do not re-derive; improve on it).** These are the splits:

| split | S-level → D-level | h_D (D-freedom, bits) | C(D \| S) (bits) | **D_def** |
|---|---|---|---|---|
| ℤ₅ Fourier (y → w → G) | w → mixed characters | 173.9 | 210.0 | **35.9** |
| ℤ₂ along an involution z of Q, given y (the 523 route) | x + x(·+z) → x − x(·+z) | 69.8 | 121.5 | **51.7** |
| ℤ₂ along z inside (y, w) | — | 14.6–16.8 | 91.7–92.1 | 74.9–77.4 |
| negation (symmetric / antisymmetric) | no S-level exists [P] | 84.0 | 275.6 | 191.6 |

What else has been tried:
- **Gauss periods.** The equations are, exactly, sums of four Hermitian norms in O_F[H], with F = ℚ(ζ₄₃)^⟨4⟩ of
  degree 6, a CM field. Only parity conditions follow [P].
- **Solver designs** (§3.3 of the answer): a 2 + 2 join per (y, w) costs 2^87; the inner ℤ₂ route costs 2^80; a
  Wagner tree over w-parents costs 2^101.
- **Other approaches:** annealing, tabu search and CP-SAT all fail or time out. Local search also fails on a control
  structure that has cores.
- **The S14 sub-slice** (an extra involution, deficit 22.7) is exactly empty at stage 1 [M].

## 3. The benchmark that worked: v = 523, order 2092 (solved here)

Read `report/paley_523_answer.md` §0, §3 and §5.

**The cost model.** A ℤ₂ = ℤ₁₈/ℤ₉ split gives an S-level (sums) and a D-level (differences), with:
- **D_def = C(D | S) − h_D = 26.1 bits**, so a core needs 2^26.1 S-solutions;
- each S-solution comes from a two-level Wagner join at about 2^18.1 pair visits, plus a 1.6 ms D-check;
- total: **log₂(pair visits per core) ≈ D_def + log₂(cost per S-solution) = 26.1 + 18.1 = 44.2**, about 80
  core-hours per core;
- the model matched production to within about 1 bit.

**The budget** is about 2^46 pair visits per core: four cores for a few days, at 22–30 ns per pair visit.

## 4. What to deliver: up to three proposals, ranked

**Each proposal must have:**
1. **Definition:** the algorithm, or the algebraic ansatz, precisely.
2. **Cost:** cost per expected core as a formula in measured quantities (h_D, C(D|S), fibre sizes, list sizes),
   evaluated at 863 stage 2b, or at the y-level, if it re-splits stages 2a and 2b.
   - Use the repository's numbers and scripts where possible: `theory/h863_split_deficits.py`,
     `theory/h863_yw_estimate.py`, `theory/h863_analog.py`, `theory/h863_lib.py` (`TableSampler`), and
     `theory/q523_gauss.py` (the lattice-Gaussian estimator).
   - Say which inputs are [M] and which are [E].
3. **Why it can beat the floor:** which assumption of the cost model it breaks, and why that is legitimate. Prove it
   if you can.
4. **A cheap test** a Claude instance can run in hours. Two test beds:
   - a planted instance (choose a random fibre point, set the targets to its autocorrelations, require recovery);
   - the known-order analogue v = 519 (G = ℤ₂² × ℤ₃ × ℤ₄₃, ×4; order 2076 is known). Its ℤ₃ split has D_def ≈ 21.7
     with a (y, w) fibre of 2^70.7 (`theory/h863_analog.py 3 43 4`), so it is in the 523 regime.
   - State the pass and kill criteria in advance.
5. **Transfer:** which other hosts it would help. The other open orders are the rich hosts in
   `report/b3_push_answer.md` §8.1 and the 863-shaped orders in `report/h863_lift_answer.md` §4, where the deficits
   are 38–90 bits.

**Directions worth considering** (you are not limited to these):
- **Representations (BJMM / May–Ozerov style).** Write each block as a sum of overlapping partial assignments, so
  that many representations of one solution exist. Can the number of S-solutions per merge grow faster than the
  deficit? The equations are quadratic, so say how the cross terms are handled.
- **Spend the surplus.** The y-fibre has 2^63.5 expected cores; only one is needed. Choose a sub-slice, a fixed
  sub-pattern, or extra algebraic structure that lowers C(D|S) faster than it lowers the expected count, ending with
  D_def ≲ 26 and E[cores per stage-1 solution] ≥ 2^10.
  - S14 did this but was empty: its stage-1 lag-0 equation has no solution in its alphabet.
  - Can you find such sub-slices that pass that lag-0 test? Check it exactly: `theory/h863_lag0_check.py` is the
    subset-sum test for the transfer hosts; adapt it to your slice.
- **Fourier side first.** Prescribe the Gauss-period values (the x̂ on the mixed characters), for instance from ideal
  factorisations in O_{F(ζ₅)} or from cyclotomic or relative-difference-set constructions. Then invert or round to
  signs with a provable or measurable success rate.
- **Linearisation or lattices.** Can an S/D split be found in which the D-level equations become linear (or affine)
  in the D-unknowns given S, so that D is found by lattice or linear algebra instead of a join?
- **A proof of impossibility.** If you believe no split of this group and slice can have D_def below about 30, a
  proof [P] for a clearly stated class of splits is a valid, useful answer.

## 5. The answer file: `report/openai_deficit_answer.md`

It has these sections:
- **Bottom line:** the best proposal, its claimed cost exponent at 863 and at 519, and confidence.
- **One section per proposal,** with items 1–5 of §4.
- **Ideas considered and rejected,** one line each with the reason, so they are not re-tried.
- **What a Claude max-effort instance should run first:** the exact commands or a precise specification, the expected
  outcome, and the kill criterion.
- **Unfinished items,** if the session ends early.

Commit with a message that names the proposals, and push to `codex/order-3452-compression`.
