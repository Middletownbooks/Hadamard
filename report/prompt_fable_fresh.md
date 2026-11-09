# Prompt (Claude Fable 5.2 in Claude Code): find something new at the open Hadamard or skew-Hadamard orders

*Self-contained, for a fresh Claude Code session on `middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`.*

## The goal

Construct a Hadamard matrix at an order where none is known, or a skew-Hadamard matrix of an order 4v that the
Đoković–Seberry census lists as undecided and nobody has settled since.

- **What counts as success:** a verified matrix. A **proved** obstruction that closes a natural family and redirects
  the field is the honest second prize. A new estimate or heuristic is not.
- **Where the room is:** this repository and its collaborators have pushed the existing toolkits about as far as they
  go on modest hardware (§2). What remains is room for a genuinely different idea. Spend real thought on that before
  spending compute, and test ideas cheaply and decisively before scaling any of them.

You choose the approach. The suggestions in §3 are not instructions.

## 1. The targets

- **Ordinary Hadamard.** 85 orders below 10000 are open: `data/open_orders_lt10000.json`, key `open`. The list is
  EKLO's (Egan, Kotsireas, London and Ó Catháin, arXiv:2610.27175v2, Table 5) with later settlements removed.
  - The smallest are 4·599 = 2396, then 2676 = 4·3·223, 2876, 2956, 3004, 3452, 3628, 3788 and 3884.
  - Every prime in the list is ≡ 3 (mod 4). Most entries are 4p; there are a few 4·3p, 4·5p, 4·7p, 4·9·191, 4·11·179,
    4·13p and 8p.
- **Skew-Hadamard.** 630 census values v < 2500 are open: `data/skew_open_census.json`, against K. Kumar's repository
  at commit d885022. Orders 8v and 16v do not settle a census value. The smallest are 89, 101, 107, 149, 153, ….
  1661 and 1963 may be settled by Schumer's 6644 and 7852; check before claiming either.
- **Check for newer work.** If the network allows, check Kumar's repository (`https://github.com/vatskishan03/hadamard`)
  and arXiv for anything newer than these lists before claiming novelty.

## 2. What is known, compressed

The detail is in the files named here. Read what you need.

**What has produced matrices.**
- **Goethals–Seidel cores** (four ±1 sequences with Σ PAF = 0) or their bordered forms:
  - Legendre quadruples (Theorem 5 of EKLO);
  - Theorem 7 (Σ PAF = −4(k−1) on a subgroup S of index k, +4 off S, with a width-4(k−1) border).
- **Restricted by multiplier groups** to slices of about 40–100 orbits per block, and found by:
  - staged compress-and-lift searches (this repository: 2092, 2292, 2588, 3732; `report/paley_523_answer.md`,
    `report/alpoge_design_answer.md`);
  - k-tree or Wagner joins on flat or low-energy pools (3532, 3676; Kumar's 2284, 2212, 4492);
  - exact quotient chains with a truncated heavy-vector last level (Kumar's 3964, 2636 and 3356, the last over the
    Galois ring GR(4, 3));
  - Paley plug-ins with a few defects, and a quaternary conference matrix route (Kumar's 643, 859, 955, 1051, 1699).
- **The August 2026 matrices of Alpöge, Voinov and Reynolds-Haertle:**
  - b = 0 cores with extremely flat spectra (Φ ≈ 0.1);
  - Legendre quadruples with no multiplier at 668 and 716, with a hidden half-shift template;
  - Theorem 7 at k = 4, 6 and 8.

  See `forensics/`, `report/alpoge_methods_forensics.md`, and the EKLO summary in `report/prompt_openai_fresh_idea.md`
  §1.
- **Skew matrices** come from the same arrays with block A skew and an odd-order multiplier group. Kumar has about 130
  of them, many in 2–14 days on a 32-core, 246 GB machine (his `skew/README.md` and `search/NOTES.md` §10–§115).

**What has been tried here and failed, with the measured reason.** Read `report/open_screen_answer.md`,
`report/h4892_answer.md`, `report/h2396_plugin_answer.md`, `report/skew_census_answer.md` and
`report/b3_push_answer.md` §0.
- **The central dichotomy.** At every open order screened (all 85 ordinary, all 630 skew), slices cheap enough to
  search turn out arithmetically empty, while slices rich enough to have solutions cost 2^54 or more per core.
- **Cheap means rigid.** Over 20 slices with positive Gaussian estimates were decided exactly at the ordinary orders,
  and all were empty, mostly at the coarsest quotient (the ℤ_q compression). Two fell to an exact energy identity.
- **Rich means expensive.** The best rich hosts (863: 2^63 expected cores per stage-1 solution; 599, 669, 751, 971)
  have no affordable final stage. Every split measured has a deficit of 36 bits or more, against about 26 that is
  feasible on 4 cores.
- **Rejected identities.** Algebraic identity classes proposed and rejected with proofs: tensor-Paley with any
  boundary, three-stratum products, elliptic and cyclotomic seeds with small repairs (`report/openai_4p3_review.md`,
  `report/openai_identity_answer.md` on branch `codex/order-3452-compression`).
- **Rejected searches.**
  - A projection (divide-and-concur) method failed its own blind calibration: 0 of 32 at p = 61
    (`logs/openai_fresh_gateA/`).
  - Paley plug-ins at 599 work at small K (a verified H(940) control), but cost at least 2^54 at the target.
  - The skew screen left two slices within budget, 361 and 325. Both were searched without a solution, 361 only to
    30% of its split values.
- **Kumar's own censuses**, which cover far more than ours, report the same picture. See `FRONTIER.md` and §102–§104
  for the ordinary orders, and §14 and §110–§115 for skew.

**Useful results to build on.**
- **Theorem B** (`report/openai_4p3_review.md` §2): one fixed border certificate completes every k = 4 Theorem-7 core
  with |K| odd. That covers 38 of the open orders, so a core alone suffices there.
- **The compression normal form** for those cores.
- **A tested toolkit** in `theory/` and `scripts/`:
  - exact stage-1 decision (`theory/open_screen_stage1_exact.py`) and calibrated Gaussian costs (`theory/q523_gauss.py`);
  - low-energy lists with a k-tree join (`theory/skew_*`, adapted from Kumar's `cp_lowe.py`, `cp_enum.c` and `kt4.c`);
  - independent verifiers (`scripts/verify_h2092.py`, `scripts/verify_skew.py`).

## 3. Where something new might come from (your call)

These are prompts for thought, not a plan:
- **The algebraic Paley point.** Every success came from structure that the Gaussian model undercounts or ignores.
  Kumar's plug-in calibration found tens of thousands of solutions where the model predicted 5, all near an algebraic
  point. Is there a family where an identity makes most of the equations hold automatically, and the remaining
  search is small and not rigid?
- **Cost exponents.** The bottleneck is a 4-list matching problem with quadratic equations, a deficit of 36 bits or
  more, and fibres of 2^35–2^174. Is there an algorithmic idea whose cost exponent is genuinely better? Representation
  methods with honest accounting, or algebraic structure in the Fourier domain (norm equations in cyclotomic fields)?
- **Other arrays, groups or recursions.** Different ones put an open order within reach: non-abelian or cocyclic
  matrices, new Theorem-7-type borders, or product theorems with a new ingredient. EKLO's Table 3 lists the smallest
  open ingredients: T-matrices 127, Williamson-type 71, TT(48), BS(45, 44), symmetric conference 146, skew 356.
- **Skew.** The census is large, and the ordinary-order rigidity may not bind the same way there. 200 of the 630 open
  values never appear in Kumar's notes.

## 4. Constraints

- **Git:**
  - Work on branch `claude/focused-planck-j1645f`. Before every push run
    `git fetch -q origin claude/focused-planck-j1645f && git merge -q --no-edit origin/claude/focused-planck-j1645f`,
    then `git push -u origin HEAD:claude/focused-planck-j1645f`.
  - Never push to `main` or to other branches. Open no pull requests. Commit trailers follow the existing commits.
  - Commit and push at least every few hours, so that an interrupted session leaves its state behind.
- **Kumar's repository:** clone it read-only if you need it. Call the `add_repo` tool for `vatskishan03/hadamard` with
  read access, then run `GIT_LFS_SKIP_SMUDGE=1 git clone --depth 50 https://github.com/vatskishan03/hadamard /home/user/vatskishan03/hadamard`.
  His code is Apache-2.0; credit it in any file you adapt.
- **Credit:** EKLO for the survey and Theorem 7, Alpöge, Voinov and Reynolds-Haertle for the twelve matrices, Kumar
  for his constructions, and the OpenAI model for Theorem B and the normal form.
- **Compute:** 4 cores and about 15 GB of RAM. Install `numpy`, `numba` and `ortools` with pip if they are missing;
  `gcc` is available.
  - Never kill processes you did not start. Never `pkill -f` with a pattern matching your own shell. Never `rm` a path
    built from an unchecked variable.
  - Make long runs restart-safe, and arm a `send_later` check-in while workers run.
- **Honesty and tags:**
  - Mark every claim [P] proved, [M] exact computation (script named) or [E] estimate.
  - A cap or a timeout is not a zero.
  - Calibrate any new search method on a solved order, blind, before trusting it on an open one, and fix its pass and
    kill criteria in advance.
- **Verification:**
  - A matrix counts only if a script sharing no code with its construction checks H Hᵀ = n·I exactly, and H + Hᵀ = 2I
    for skew.
  - Check novelty against the lists in §1 and any newer source you can reach.
  - Commit any verified matrix at once, with a bundle in `pub/` in the format of `pub/H2092/`.
  - Do not post outside the repository; the owner decides that.
- **Time:** this is a long autonomous session. Plan your own phases. Stop when you have a verified result, a proved
  closure worth recording, or when the remaining ideas have failed their cheap tests. In any case stop after about
  72 hours of wall-clock time, unless a calibrated production run is close to its predicted success.

## 5. Deliverable

`report/fable_fresh_answer.md`, written for a mathematician who knows this field:
- the bottom line (matrices, proofs, or neither);
- the ideas you considered, with the reason each was dropped or pursued;
- for each idea pursued: its definition, its calibration, the decisive tests and their results;
- anything built, with its verification and novelty check;
- what you would do next with more compute or time.

Keep it as long as the content needs and no longer.
