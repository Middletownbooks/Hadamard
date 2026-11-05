# Prompt (high effort, bounded): screen all 85 open orders below 10000 for the cheapest staged slice

*Self-contained, for a fresh Claude Code session on `middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`.*
This is a ranking run with existing tools, followed by exact stage-1 decisions for the few best candidates. It is
not a production search.
- **Time box:** about 30 hours of wall-clock time.
- **Order:** commit per batch of orders, so that a cut-off session leaves a usable partial table.

**Read first:**
1. this file;
2. `README.md`;
3. `report/b3_push_answer.md` §0 and §6–§8 (the k4 census, calibration, and why rich hosts fail);
4. `report/paley_523_answer.md` §0, §3.3–3.6 and §5 (the deficit cost model that built H(2092));
5. `report/alpoge_design_answer.md` §4 (H(3732): a b = 0 GP slice of ℤ₃ × ℤ₃₁₁ with a dual-first deficit of 7.4 bits);
6. `report/h863_lift_answer.md` §0 and §5 (the regime a staged search needs);
7. `report/openai_4p3_review.md` §0 and §4 (Theorem B, and why seed families fail).

## Rules

- **Pushing:** before every push run
  `git fetch -q origin claude/focused-planck-j1645f && git merge -q --no-edit origin/claude/focused-planck-j1645f`,
  then `git push -u origin HEAD:claude/focused-planck-j1645f`.
- **Branches:** never push to `main` or to `codex/order-3452-compression`. Open no pull requests.
- **Commits:** trailers follow the existing commits. Logs go in `logs/open_screen/`, scripts in
  `theory/open_screen_*.py`. Don't commit raw worker output every few minutes.
- **Processes:** never kill processes you did not start. Never `pkill -f` with a pattern that matches your own
  shell; kill by PID.
- **Restart safety:** runners resume from their output files. Arm one `send_later` check-in while workers run.
- **Reuse:** use the repository's screens and estimators (§2) before writing anything new. Credit their authors
  (Kumar and EKLO for the bordered k4 and ℤ_q families; the OpenAI model for Theorem B and the normal form).
- **Tags:** mark claims [P] proved, [M] exact computation (script named), [E] estimate. A cap or a timeout is never a
  zero.
- **Compute:** 4 cores and about 15 GB of RAM. Install `numpy`, `numba` and `ortools` with pip if they are missing.

## 1. Why this screen

The repository's systematic screens covered the open orders **up to 4000** (599–971).
- Above 4000 there are 76 open orders, which have had only one-off looks:
  - 4412 and 4748;
  - the ℤ_q census of the composite orders (`report/zq_census_answer.md`);
  - the 863-shaped k4 transfer tables (`report/h863_lift_answer.md` §4, `report/blockswap_answer.md` §3);
  - the n = 4p + 3 algebraic work.
- Below 4000, every route is closed or out of reach.
- **Cheapness depends on factorisation, not size.** H(3732) came from a slice with deficit 7.4, against about 26 at
  523, which was barely feasible. Whether such a slice exists depends on how v, v − 1 and v − 3 factor, and on the
  cyclotomic subgroups available, not on the order's size. An order above 4000 may have a cheaper slice than anything
  at 599–971.
- **The goal is a ranked table.** For every open order, the best candidate slice per family, with exact gate status
  and a calibrated cost per core. Then exact stage-1 decisions for the top few.

## 2. Families, and the tools for each

The open orders 4v are in `data/open_orders_lt10000.json` (key `open`, 85 entries).

| family | group | applies when | screens and estimators already in the repository |
|---|---|---|---|
| **b = 0**: Goethals–Seidel core, Σ PAF = 4v·δ | abelian, order v | v odd | `scripts/b0_open_screen.py` and `scripts/gs_twist.py` (prime v); `theory/ab_sweep.py`, `theory/ab_s1_zscreen.py`, `theory/ab_s2_*_screen.py`, `theory/ab_gauss.py` (GP slices of composite v, as at 933); `theory/cost_sweep.py` |
| **b = 1**: Legendre quadruple, bordered GS | order v − 1 | always | `theory/p571_msweep.py` (all-subgroup sweep); `theory/gp_screen_candidates.py`, `theory/gp_z2screen.py`, `theory/gp_alpha_screen.py`, `theory/slice_compress_cpsat.py`; `theory/q523_ycomp.py` and `theory/q523_yscreen.py` (stage-2 deficit) |
| **k4**: Theorem 7, t = 3, border of width 12 | order v − 3, K of index 4 | v ≡ 3 (mod 4) | `theory/bt_census_run.py`, `theory/bt_census2.py`, `theory/bt_thm_alpha.py`, `theory/bt_thm_rigid.py`; `theory/h863_analog.py`, `theory/h863_split_deficits.py` (deficits for ℤ₂² × ℤ_m × ℤ_p); Theorem B gives the border for odd \|K\|, i.e. v ≡ 7 (mod 8) |
| **ℤ_q all-character** | order v − q, q \| v | v composite | `theory/zq_census.py`. Already censused: take its survivors from `report/zq_census_answer.md` §3 |
| **Cooper–Wallis** (T-matrices × Williamson) | t \| v, Williamson order w = v/t | v composite | `scripts/tmat_estimate.py` (with `rowsum_feasible`); `report/tmat_routes.md` |

If a tool is limited to a range or a group shape, extend it minimally, and say so.

## 3. Method

**Step 1. Enumerate candidate slices per order and family.**
- Multiplier or GP slices with about 40–160 orbits per block.
- For staged slices: those with a natural compression along a large prime factor (stage 1 on a small group, stage 2
  on cyclotomic classes), as at 523, 647 and 933.
- Cap the work per order; record any order where the cap binds.

**Step 2. Exact cheap gates, cheapest first.**
- lag-0 composition and row-sum feasibility;
- the ℤ₂ screen;
- α-congruence;
- coset-Parseval;
- Theorem α and Lemma R (k4);
- Q-level emptiness (k4; closed form by the normal form when |K| is odd);
- `rowsum_feasible` (T-matrices).

**Step 3. Calibrated cost for each survivor** [E], read at the true target, never at the mean.
- **Expected count:** E per stage-1 solution, or E_Q.
- **Stage-2 constraint:** C₂ in bits.
- **Best split:** the lowest D_def over the splits the tools support, and its D-fibre size.
- **Estimated cost per core:** log₂(cost per core) ≈ D_def + log₂(per-S-solution cost). The 523 values were 26.1 +
  18.1 = 44.2 bits. The budget is about 46 bits (four cores for a few days).
- **Regime test** (`report/h863_lift_answer.md` §5):
  - D_def ≲ 26;
  - D-fibre ≲ 2^40;
  - E per stage-1 solution ≥ 2^10;
  - C₂ ≲ 150 bits.

**Step 4. Rank and commit.** Write `logs/open_screen/table.jsonl`, one line per (order, family, best slice), with:
- the gate status;
- r, E, C₂, D_def, D-fibre and the cost exponent;
- the next exact gate and its cost.

**Step 5. Exact stage 1 for the top candidates** (at most 5, and only those in or near the regime).
- Decide stage 1 exactly, by exhaustive join or by CP-SAT with a second, independent method, as
  `theory/h863_s14_stage1_exact.py` and `theory/blockswap_stage1_exact.py` did.
- Planted recovery is required.
- Recent history: S14 at 863 (2^45) and the block-exchange slice (2^33 per stage-1 solution) both died here in
  minutes. So **no cost estimate counts until stage 1 is decided.**
- For a candidate whose stage 1 is non-empty, re-read its deficit at real stage-1 solutions.

**Step 6. Stop.** Do not build a stage-2 solver here. A candidate that passes Step 5 and stays in the regime gets a
separate targeted session.

## 4. Deliverable: `report/open_screen_answer.md`

It contains:
- **Bottom line:** whether any open order has a slice in the regime with a non-empty stage 1, and which.
- **The ranked table:** for every order, the best slice per family, with gates, estimates and cost, and caps marked.
- **Step 5 results:** exact decisions, with methods and planted controls.
- **The verdict on search:** either named targets for a staged solver, or the statement (with its evidence) that
  every open order's best slice is out of the regime or empty at stage 1.
- **The cut-off point** and what was not covered.

Commit, push, and stop.

## 5. Calibration warnings (from this repository)

- Zero counts are common below about 2^30. At the open orders, every exactly decided positive estimate has so far
  been zero, up to 2^85. Estimates rank; only exact gates decide.
- Small-orbit slices are the cheap ones, and they are the ones that most often die arithmetically. Decide stage 1
  exactly before believing a cost.
- A real stage-1 solution moves the estimates by up to about 1 bit (523 §5). A synthetic y read at the mean can be
  wildly wrong (519 in `report/openai_deficit_answer.md`: estimated deficit 10, exactly empty).
