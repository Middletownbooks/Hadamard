# Prompt (max effort, targeted): can the ℤ₄₃-compressed solutions of the 863 S7 host be lifted to a core? Exact gate, lift design, go/no-go

*Self-contained, for a fresh Claude Code session on `middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`.
This is a **targeted reasoning task, not a production run**:
- answer the three questions of §3, with exact gates and calibrated measurements;
- write the answer file;
- stop.

If the answer is "go", the production run is a separate, lower-effort session. You write its plan here; you do not
launch it.*

**Read first:**
1. this file;
2. the two Codex reports on the branch `codex/order-3452-compression`: `report/h3452_stage1_20261103.md` and
   `report/h3452_lift_20261103.md`;
3. `report/paley_523_answer.md` §0, §3 and §5: the staged solver that built H(2092), and its deficit and cost model;
4. `report/b3_push_answer.md` §0, §3, §6 and §8: the bordered family, the exact gates, the estimator calibration, and
   the cost of the rich hosts, including 863.

## Rules

- **Pushing:** before every push run
  `git fetch -q origin claude/focused-planck-j1645f && git merge -q --no-edit origin/claude/focused-planck-j1645f`,
  then `git push -u origin HEAD:claude/focused-planck-j1645f`.
- **Branches:** never push to `main` or to `codex/order-3452-compression`. Open no pull requests.
- **Commits:** trailers follow the existing commits. Commit small, with logs.
- **Processes:** never kill processes you did not start. Never `pkill -f` with a pattern that matches your own shell;
  kill by PID.
- **Credit:**
  - Codex's work: the construction set-up, the scripts `theory/h3452_*.py` and the certificates are by OpenAI Codex,
    directed by the repository owner. Say so in the headers of any file you derive from them.
  - Kumar and EKLO: the bordered construction (`theory/bt_engine.py`, `bt_border.py`) is based on their work; keep the
    existing credits.
- **Bringing in the Codex files:**
  - Run `git fetch origin codex/order-3452-compression`.
  - Then run
    `git checkout origin/codex/order-3452-compression -- theory/h3452_*.py data/h3452_*.json logs/h3452 report/h3452_*.md`.
  - Commit them unchanged, with a message that credits Codex. Then work on top of them.
- **Verification:** a full core counts only if verified exactly (integer PAF on all 860 lags, then the bordered
  H Hᵀ = 3452·I) by a script sharing no code with the search.
- **Compute:**
  - 4 cores and about 15 GB of RAM. Install `numpy`, `numba` and `ortools` with `pip` if they are missing.
  - Keep measured compute to about 12 hours in total.
  - Arm one `send_later` check-in while workers run.

**Tags:** mark claims [P] proved, [M] exact computation (script named) and [E] estimate. UNKNOWN from a solver decides
nothing.

## 1. The set-up (exact; from the Codex reports)

**The family and the host.**
- **Family:** k4 bordered construction (Theorem 7 of EKLO, Kumar's form), t = 3, border of width 12.
- **Group:** G = ℤ₂² × ℤ₅ × ℤ₄₃ (|G| = 860), K = ℤ₅ × ℤ₄₃, Q = ℤ₂².
- **Target:** four ±1 sequences on G with Σᵢ PAFᵢ = 3440 at 0, −12 on K \ 0 and +4 off K.
  - Equivalently Σᵢ |x̂ᵢ(χ)|² = 3452 for every character χ nontrivial on K, 12 on the three nontrivial characters of Q,
    and Σ rᵢ² = 3452 at the trivial character.
- **Host S7 (census key `2e93be9718324148`):** the multiplier is the identity on ℤ₂² × ℤ₅ and ×4 on ℤ₄₃, of order 7.
  - Per block there are 140 orbits: in each of the 20 cells h ∈ H = ℤ₂² × ℤ₅, a singleton (h, 0) and the six classes
    h × Cⱼ, where the Cⱼ are the cosets of ⟨4⟩ in ℤ₄₃^*.
  - −1 ∉ ⟨4⟩, so the classes pair under negation: (1, 2), (3, 6), (7, 9) by representatives.
  - There are 72 lag classes under S and negation, so 71 equations.
  - Census estimate log₂ E_Q = 136.9 [E].

**Codex's certificates.**
- **Level H = G/ℤ₄₃.** `data/h3452_stage1.json`, verified on all lags:
  - y_i(h) = e_i(h) + 7 Σⱼ c_ij(h), with c_ij(h) the class signs and e_i(h) the singleton sign;
  - e is fixed by y mod 14;
  - the H-level targets are 2936 at 0, −516 at (0, u ≠ 0) and 172 off ℤ₅;
  - rows (4, 6, 6, 58).
- **Border.** `logs/h3452/border.json`, checked exactly. It works for any core with these Q-compressions.
- **Necessary condition on the lags (h, 0).** Codex derived Σ_{i,j} PAF(c_ij)(h) = 480·δ on H.
- **The complementary projection along H, onto ℤ₄₃.** `data/h3452_prime_projection.json`, verified. Per block it is
  a 20 × 6 sign table with both margins fixed; Codex found compatible tables.
- **Lift heuristics.** Annealing, tabu search and rectangle switches all plateau: the best attempt fails 718 of 860
  lags. A full Boolean model returned UNKNOWN at 60 s.
- **The S14 sibling** (`fbf4f510b2566552`): fibres reported INFEASIBLE, but **not** independently confirmed.

**Measured on the branch since then (2026-11-03).**
- **Fibre estimate [E].** `theory/h863_fibre_estimate.py`, `logs/h863/fibre_estimate_codex_cert.jsonl`, using the
  repository's lattice Gaussian:
  - the fibre of the Codex y holds 2^339.1 sign assignments;
  - the remaining equations have rank 59, the target lies in their lattice, and the Mahalanobis term is 0.15;
  - they cost 275.6 bits, so the fibre holds **about 2^63.5 expected cores**.
  - So the heuristic plateau is no evidence of emptiness: annealing also failed on the 3964 control, which has cores.
  - But at the open orders, every positive estimate decided exactly so far was a zero. They were killed by exact gates
    at an intermediate level (`report/b3_push_answer.md` §0, §6.3).
- **The one such level left inside this fibre is G/ℤ₅ = ℤ₂² × ℤ₄₃.**
  - `theory/h863_g5_gate.py` is an independent CP-SAT model of it, joint with the fixed y. It returned UNKNOWN after
    1200 s on 4 workers (`logs/h863/g5_gate_codex_cert.json`).
  - The model: per block i and cell (q, u), six class signs with fixed sum; w_i(q, j) = Σ_u c_ij(q, u).
  - Its targets on ℤ₂² × ℤ₄₃ are 3392 at 0, −60 at (0, t ≠ 0) and 20 at q ≠ 0.
- **Cost of the 523-style route [E].** The 523 stage 2 had lattice rank 35 and C₂ = 146 bits, and cost about 2^44
  pair visits (about 80 core-hours per core).
  - Here the lift is 59 equations and 275.6 bits.
  - The ℤ₂ U/W-split scan (`report/b3_push_answer.md` §8.5) gives about 2^−61 completions per U-solution at this host,
    against about 2^−26 at 523.
  - So a direct 523-style S/D solver is out of reach. A deeper chain is needed.

## 2. Why this is worth one targeted session

- **No other host gets this far.** This is the only rich host at an open order with an intermediate level
  actually solved, a border in hand, and a large estimate per fibre.
- **Same H at other open orders.** The same H = ℤ₂² × ℤ₅ (with ℤ_p in place of ℤ₄₃) occurs at 1223 (p = 61),
  1423 (71), 1583 (79), 1783 (89) and 2063 (103). H = ℤ₂² × ℤ₇ occurs at 1039, 1319, 1487, 1879 and 2215. Their lifts
  are larger, so a solver must first work at 863.
- **Either answer is informative.**
  - A proved obstruction, or a measured cost far beyond reach, effectively closes the large-prime compression route
    for the k4 family.
  - A "go" gives a production plan.

## 3. Questions, in order

**Q1. Is the G/ℤ₅ level inside the Codex fibre empty? Decide it exactly.**
- **Methods.** Use whatever decides it:
  - a stronger CP-SAT formulation (symmetry breaking; the Z₄₃ projection of `data/h3452_prime_projection.json` as an
    optional added constraint, recorded separately because it restricts the fibre);
  - an exact enumeration. Per (block, q), the w-options are the column-sum vectors of 5 × 6 sign tables with the
    fixed row sums; join them on the 28 G/ℤ₅ equations, using the ℤ₂² character split;
  - a congruence argument [P].
- **If FEASIBLE:** verify the w on explicit length-172 sequences with independent code, and keep several solutions.
- **If INFEASIBLE:**
  - confirm with a second, independent method;
  - test whether the obstruction is structural, i.e. whether it follows from y mod a small modulus, or from the rows
    (4, 6, 6, 58);
  - generate further H-level solutions with the Codex pipeline (other seeds and the other row patterns with
    Σ r² = 3452) and test them;
  - try to prove the level empty for the whole host [P].

  An exact kill of the host's G/ℤ₅ level is a complete answer: record it and go to §4.

**Q2. The lift chain H → G/ℤ₅ → G: design and measure it [M sizes, E counts].**
- **Sizes and counts.** For a pair (y, w) from Q1, the remaining fibre is, per (block, q), a 5 × 6 sign table with
  both margins fixed. The singleton column is fixed by e. Compute:
  - log₂ of the fibre at each stage;
  - the rank and lattice-Gaussian bits of the remaining 48 or so mixed equations, i.e. the characters nontrivial on
    both ℤ₅ and ℤ₄₃;
  - the expected cores per (y, w);
  - the expected number of w per y.

  Use `bt_gates.stratified_bits` or `q523_gauss.gauss_bits`, with fibre samplers for the doubly fixed tables.
- **Solver design.** Design the final-stage solver and give its cost per expected core. Candidates:
  - the ℤ₂² character split (four q-tables per block, combined through the characters of ℤ₂²);
  - an S/D split on the negation pairs of classes;
  - a Fourier treatment along ℤ₅ per q;
  - the Gauss periods of the degree-6 field ℚ(ζ₄₃)^⟨4⟩, where the equations are norm conditions;
  - a Wagner tree on the resulting lists, as at 523.
- **Deficit.** State the deficit in the 523 form, D_def = C(D|S) − h_D, or its analogue, and the cost formula.
- **Validate the estimator and the solver on smaller analogues with the same multiplier structure:**
  - n = 519: G = ℤ₂² × ℤ₃ × ℤ₄₃, ×4 on ℤ₄₃. 4n = 2076 is a known order; this shape is the closest calibration;
  - n = 583: G = ℤ₂² × ℤ₅ × ℤ₂₉, ×16 of order 7 on ℤ₂₉;
  - planted solutions: plant a random full assignment, set the targets to its PAFs, and require recovery.

  Where the analogue slice has cores, compare cores found against expected [E vs M], as §3.5 of the 523 report does.

**Q3. Go or no-go.**
- **"Go"** requires all of:
  - a measured cost of at most about 2^46 pair visits per expected core, at measured speeds on this machine (523 was
    2^44, about 80 core-hours per core);
  - a calibration that held on an analogue or a planted test;
  - Q1 feasible for enough (y, w) pairs to cover the expected number of cores needed.
- **Otherwise "no-go"**, with the binding number: which deficit or list size exceeds what.
- **If "go":** write the production plan, i.e. scripts, restart safety and expected wall-clock time. Do not run it
  beyond a 1-hour smoke test.
- **Any core found along the way:** assemble it with the Codex border, verify it independently, commit it at once,
  and write `pub/H3452/` in the format of `pub/H3732/`, with the border written out.

## 4. Deliverable

`report/h863_lift_answer.md` contains:
- the bottom line (Q1 status [M or P]; Q2 numbers; Q3 verdict);
- the methods and scripts (`theory/h863_*.py`) and logs (`logs/h863/`);
- the calibrations;
- the transfer to the other orders of §2: for each, the analogous fibre and deficit estimate, if Q2 gives a formula;
- the next computation.

Commit and push it, then stop.

## 5. Calibration warnings (from this repository)

- Estimates are not evidence of solutions. Zero counts are common below about 2^30. At the open orders, every exactly
  decided positive estimate has so far been zero, up to 2^85.
- Lattice obstructions (congruences mod 2^k or mod |S|) are common. Test the quotient equations modulo small numbers
  before building large solvers (Kumar `search/NOTES.md` §102).
- Local search is not a substitute for an exact or staged solver (`report/b3_push_answer.md` §8.4).
- A planted-solution test checks correctness, not cost. Measure costs on real targets.
