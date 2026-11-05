# Prompt (high effort, gated): the block-exchanging ansatz (A, B, C, τC) at v = 863: exact stage 1, border, real-y costing, go/no-go

*Self-contained, for a fresh Claude Code session on `middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`.
This tests one proposal from the OpenAI model, `report/openai_deficit_answer.md` (branch `codex/order-3452-compression`,
commit 5b26dae).*

The work is a fixed sequence of exact gates with kill criteria (§3). Stop at the first kill, write the answer and push.
Do **not** build a production lift solver in this session. On a "go", write its specification (§3, Gate 4).

**Read first:**
1. this file;
2. `report/openai_deficit_answer.md`, Proposal 1 (§1–§4) and "What a Claude max-effort instance should run first";
3. `report/h863_lift_answer.md` §0, §1.1, §3.2–3.5. These give the S7 host, the deficit cost model, and the S14 sibling
   killed at stage 1;
4. `report/paley_523_answer.md` §0 and §3.6. This is the cost model that worked: D_def 26.1, about 2^44.2 pair visits
   per core.

## Rules

- **Setup:** run `git fetch origin claude/focused-planck-j1645f codex/order-3452-compression`, then
  `git merge --no-edit origin/codex/order-3452-compression` into `claude/focused-planck-j1645f`. A trial merge was
  clean. This brings in `theory/openai_deficit_*.py`, `logs/openai_deficit/` and the report.
- **Pushing:** before every push run
  `git fetch -q origin claude/focused-planck-j1645f && git merge -q --no-edit origin/claude/focused-planck-j1645f`,
  then `git push -u origin HEAD:claude/focused-planck-j1645f`.
- **Branches:** never push to `main` or to the codex branch. Open no pull requests.
- **Commits:** trailers follow the existing commits. Commit after each gate, with logs in `logs/blockswap/`.
- **Processes:** never kill processes you did not start. Never `pkill -f` with a pattern that matches your own shell;
  kill by PID. Arm one `send_later` check-in while workers run.
- **Credit:** the ansatz and its cost formula are by the OpenAI model (codex branch); the 863 set-up and certificate by
  OpenAI Codex; the bordered construction is Kumar's form of Theorem 7 of EKLO. Say so in the file headers.
- **Tags:** mark claims [P] proved, [M] exact computation (script named), [E] estimate. UNKNOWN from a solver decides
  nothing. A timeout is never an emptiness claim.
- **Verification:** a full core counts only if verified exactly by a script sharing no code with the search: integer
  PAF on all 860 lags, then the bordered H Hᵀ = 3452·I.
- **Compute:**
  - 4 cores and about 15 GB of RAM. Install `numpy`, `numba` and `ortools` with pip if they are missing.
  - Total time box: about 12 hours of wall-clock time. Individual CP-SAT calls get at most 300 s unless stated
    otherwise.

## 1. The set-up (exact)

**The family and the host.**
- k4 bordered family, t = 3. G = ℤ₂² × ℤ₅ × ℤ₄₃, K = ℤ₅ × ℤ₄₃, Q = ℤ₂²; elements are (q₁, q₂, u, t).
- Target: Σᵢ PAFᵢ = 3440 at 0, −12 on K \ 0 and +4 off K, which gives H(3452) with a border of width 12.
- S7 host: ×4 on ℤ₄₃, of order 7. Per block there are 20 cells h = (q, u) ∈ H = ℤ₂² × ℤ₅, each holding a singleton
  and six class signs.
- The H-level (ℤ₄₃-compressed) alphabet is y = e + 7s, with e = ±1 and s ∈ {−6, −4, …, 6}. That is
  y ∈ ±{1, 13, 15, 27, 29, 41, 43}.
- H-level targets: Σᵢ Σ_{(q,u)} yᵢ(q, u) yᵢ(q ⊕ a, u + b) is 2936 at (a, b) = (0, 0), −516 at a = 0, b ≠ 0, and 172
  otherwise.

**The ansatz (OpenAI model, Proposal 1).**
- Let τ(q₁, q₂, u, t) = (q₂, q₁, −u, t). This is an automorphism of G; it commutes with S and maps K to K. Write
  (τX)(g) = X(τg).
- The four blocks are **x₀ = A, x₁ = B, x₂ = C, x₃ = τC**, with τA = A and τB = B, and C free (S7-invariant only).
- [P] R_{τC}(g) = R_C(τg), so the target reads R_A + R_B + R_C + R_C∘τ = F. Both sides are τ-invariant.
- This is **not** the S14 sibling, which required all four blocks τ-invariant and is exactly empty
  (`report/h863_lift_answer.md` §3.5). Neither slice contains the other.
- **Stage 1** has 42 independent H-cells: 11 τ-orbits each for A and B, and 20 for C. In the lag-0 energy they have
  weights four 1s and thirty-eight 2s, with excess energy Σ(y² − 1)·weight = 2856.

**What has been run already [M]** (`theory/openai_deficit_symmetry.py`; reproduced on 2026-11-05 in 0.07 s):
- at 863: 48 Q-level tuples, of which 40 survive the joint Q-row-sum / energy gate;
- at 519 (G = ℤ₂² × ℤ₃ × ℤ₄₃): 24 tuples, 0 survive. The ansatz is dead at 519.
- A second implementation, `theory/openai_deficit_controls.py`, agrees.
- Planted switch-face controls recover their plants at m = 5 and m = 3. These check the algebra, not the cost.

**Synthetic estimate [E]** (`theory/openai_deficit_estimate.py`, read at the sample mean on synthetic y that satisfy
only lag-0):
- seed 863: D_def = 23.7, final S-merge F = 49.2, D-generator entropies 8.7 / 8.9 / 17.2, 2^32.9 expected cores per y,
  conditional cost exponent 45.1;
- seed 864: D_def = 24.4, F = 49.0, 2^32.3 expected cores per y, exponent 45.6.
- The cost formula is T_core = (1/p_D)·[P + 2^(F+1−n) + L_A·L_B + L_C], with sides 2^n assumed to be n = 29 and the
  parent cost P = 0. The 2^46 budget is four cores for a few days.
- **This is not evidence yet.** At 523, conditioning on a real y moved the estimates. At 519, a synthetic estimate of
  D_def 10 belonged to an empty slice.

## 2. Why this is worth one session

- **Within budget on paper.** This is the first proposal in the repository whose estimated cost at an open order
  is inside the 2^46 budget, even if only conditionally.
- **Cheap to decide.** Its first gate, an exact stage 1, takes minutes to hours with the method that decided S14.
- **Either outcome is clean.** Either we get real y certificates and real numbers, or a proof-quality kill.

## 3. Gates, in order (stop at the first kill)

**Gate 0. Independent re-derivation (≤ 1 h).** Write your own code, sharing nothing with the `openai_deficit_*`
scripts, that:
- checks the identity R_{τC}(g) = R_C(τg) on random sequences over G (integer PAFs);
- recomputes the 48 Q-level tuples and the 40 joint-gate survivors at 863 (and 0 of 24 at 519) from the definitions:
  - Q-level targets Σ PAF = 3452·δ;
  - τ-invariance of A and B on Q means a(01) = a(10);
  - τC's Q-compression is C's with q = 01 and q = 10 exchanged;
  - odd Q-cell sums;
  - energy reachable in the alphabet with the weights above.
- Normalize by block sign and the interchange of A and B, and say so.

Kill if the identity fails, or if the survivor set differs in a way that empties the gate. Otherwise record agreement
or the differences, and continue on the union of survivors.

**Gate 1. Exact stage 1 for every survivor (target ≤ 4 CPU-hours).**
- **Primary method: exhaustive, as `theory/h863_s14_stage1_exact.py` did for S14.** For each Q tuple:
  - list every τ-invariant A (11 orbit values) and B with the block's Q row sums, under the energy bound
    (block energy ≤ 2856 − the other blocks' minima);
  - list every free C (20 cells), whose energy counts twice;
  - join A + B against C on the exact 20-lag H-level vector. C's contribution is PAF_C + PAF_C∘τ. Do **not** treat τC
    as an independent fourth list;
  - re-check every hit directly on four explicit length-20 arrays.
  - If a list exceeds memory, split it by a cheap exact invariant (for example the lag-0 energy, or a Q-character
    sub-vector) and join per part.
- **Second method, sharing no code.** Either CP-SAT with exact AddMultiplicationEquality on all 20 lags (one worker,
  seed 863 + i for tuple i, 300 s per call), or a second exhaustive implementation. Use it on every tuple the primary
  method decides, or on a large random subset, and on every hit.
- **Planted test:** plant a random ansatz y, set the targets to its PAFs, and require recovery by both methods.
- **Kill:** every survivor exactly empty by the primary method, confirmed by the second. This is an exact theorem for
  the ansatz at 863 [M]. Record it and go to §4.
- **Pass:** at least one verified real y. Aim for 4 or more distinct y, across tuples if possible. Save them to
  `data/blockswap_y.json`.
- **Operational stop:** neither a y nor a complete decision within the time box. Record what was decided and what was
  not, and stop.

**Gate 2. Border (≤ 1 h).**
- **The need.** The k4 border depends only on the Q-compressions of the four blocks, in order (A, B, C, τC). For each
  passing y, find a border by `theory/bt_border.py`, or reuse `logs/h3452/border.json` if the Q-compressions are equal
  to those it was built for, up to an equivalence you prove.
- **Kill:** a Q tuple with no border is useless; drop its y. Kill the ansatz only if no passing y has a border.

**Gate 3. Real-y costing (target ≤ 4 CPU-hours).** Follow the OpenAI specification (its Proposal 1 §4, "Second
gate").
- **Estimator on real y.** Replace the synthetic E and S in `theory/openai_deficit_estimate.py` by E = y mod 14
  (representatives ±1) and S = (y − E)/7. Sample uniformly in the real y-fibre under the ansatz, with seeds 863, 864
  and 865 at 8000 samples each.
- **Read the levels at the true target** with `theory/q523_gauss.gauss_bits`:
  - all levels;
  - the S-level of the z = (1, 1) split (z is τ-fixed);
  - the parent level (also sum over ℤ₅).
- **Calibrate the tails.** Check the low-dimensional tails against exact hit counts, as `theory/h863_tail_calib.py`
  did.
- **Report:** log₂ E[cores | y], C(all), C(S), C(parent), F, h_D per generator, D_def, and the target-in-lattice test.
- **Measure, don't assume:**
  - real parent triples: join A-parent + B-parent against C-parent + τC-parent on exact signatures;
  - for a fixed matching middle signature, the distinct descendant S-candidates per side, deduplicated, and the visits
    spent generating them. This gives the achievable n;
  - actual D-list sizes L_A, L_B and L_C for sampled S-solutions, not the unconditional h_D;
  - the D-completion probability p_D, from the estimator.
- Then evaluate T_core = (1/p_D)·[P + 2^(F+1−n) + L_A·L_B + L_C] with measured P and n. Use
  (2^n_L + 2^n_R)/2^(n_L+n_R−F) if the sides are unbalanced.
- **Pass:** all of
  - log₂ E[cores | y] ≥ 10;
  - D_def ≤ 25;
  - tails within 1 bit where exact hits are numerous;
  - log₂ T_core ≤ 45, i.e. 2^46 with at least 1 bit of margin, **including** measured parent, list, D-check and memory
    costs.
- **Otherwise kill the solver design,** naming the binding number. The ansatz may stay mathematically alive.

**Gate 4 (only on a pass). The solver specification, not the solver.** Write `report/blockswap_solver_spec.md`:
- the three-generator Wagner tree (parents, middle buckets, the S-join, the A + B versus C D-check);
- memory layout, restart safety and planted-test design. Plants must not be forced into privileged middle buckets;
  report inclusion probability and recovery separately;
- the expected wall-clock time on 4 cores.

The production build is a separate session.

**Gate 5 (cheap, always). Transfer screen (≤ 2 h).**
- Run the Gate 0 Q/energy screen, with the analogous involution, at the 863-shaped hosts of
  `report/h863_lift_answer.md` §4: 583 as an analogue, then 1583, 2063, 1039, 1487 and the others in that table.
- Then the Gate 1 lag-0 subset-sum test (adapt `theory/h863_lag0_check.py`).
- Report the survivors per host. No deeper gates.

## 4. Deliverable

Write `report/blockswap_answer.md` with:
- the bottom line (the gate reached, the verdict, and the binding number);
- per-gate methods, scripts (`theory/blockswap_*.py`), logs and results with tags;
- any y certificates with their independent checks;
- the transfer table;
- the next computation.

Commit, push, and stop.

## 5. Calibration warnings (from this repository)

- Estimates are not evidence of solutions. Zero counts are common below about 2^30. At the open orders, every exactly
  decided positive estimate has so far been zero, up to 2^85, and S14 at 2^45 was decided exactly as zero at stage 1.
- Read estimates at the true target on real y, never at a synthetic mean, before using them for a decision.
- Sides of 2^29 distinct S-candidates per middle value must be measured. At 523, a middle value had only a few hundred
  parent pairs, which cost 3.5 bits against the ideal Wagner cost (`report/paley_523_answer.md` §3.3).
- A planted test checks correctness, not cost.
