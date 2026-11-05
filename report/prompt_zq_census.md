# Prompt (high effort, budget-bounded): census of the ℤ_q all-character bordered family at the 19 composite open orders

*Self-contained, for a fresh Claude Code session on `middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`.
Model usage is near its weekly ceiling: keep reasoning short, let the computers work, commit often, and stop at the
time box (§5).*

**Read first:**
1. this file;
2. `report/b3_push_answer.md` §0, §2, §3 and §4 (the engine, the exact gates, the theorems, and the 669 results under
   the family `z3all`);
3. `README.md` (what is built; the open orders).

Kumar's repository (`https://github.com/vatskishan03/hadamard`, Apache-2.0) has the original construction. Clone it
read-only only if needed: call the `add_repo` tool for `vatskishan03/hadamard` with read access, then run
`GIT_LFS_SKIP_SMUDGE=1 git clone --depth 50 https://github.com/vatskishan03/hadamard /home/user/vatskishan03/hadamard`.
The relevant material there is `search/NOTES.md` §18–19 and `hadamard/2148/README.md`.

## Rules

- **Pushing:** before every push run
  `git fetch -q origin claude/focused-planck-j1645f && git merge -q --no-edit origin/claude/focused-planck-j1645f`,
  then `git push -u origin HEAD:claude/focused-planck-j1645f`.
- **No `main` and no pull requests:** never push to `main` (the owner fast-forwards it) and open no pull requests.
- **Commits:** trailers follow the existing commits. Commit small: one commit per order or per family, with logs.
  Do not commit raw worker logs every few minutes.
- **Processes:** never kill processes you did not start. Never `pkill -f` with a pattern that matches your own shell;
  kill by PID.
- **Restart safety:** runners resume from their output files. Arm one check-in (`send_later`) while workers run.
- **Reuse:** the `theory/bt_*` engine before writing anything new. Credit Kumar for the construction in code headers.
- **Verification:** a construction counts only if verified exactly (integer PAF, H Hᵀ = 4n·I) by a script sharing no
  code with the search.
- **Compute:** 4 cores, about 15 GB of RAM. Install `numpy`, `numba` and `ortools` with pip if they are missing.

## 1. The family (exact)

Kumar `search/NOTES.md` §18, in the form that built H(2148) for n = 537 = 3 · 179.

- **Data:**
  - G abelian of order n − q, for an odd q dividing n;
  - K ≤ G of index q, with Q = G/K cyclic of order q (for q = 9 also Q = ℤ₃²);
  - T = all characters of Q (for odd q the trivial character must be in T), so t = q.
- **Sequences:** four ±1 sequences on G with

      Σᵢ PAFᵢ(g) = −4q  for g ∈ K \ 0,      Σᵢ PAFᵢ(g) = 0  for g ∉ K.

- **Completion:** a border of size 4q, constant on the K-cosets, completes the Goethals–Seidel array over G to an
  H(4n). This needs X, R̂ and Ĉ to be Hadamard of order 4q and X Ĉᵀ + R̂ Mcᵀ = 0, where Mc is the Goethals–Seidel
  array over Q of the compressions.
- **Consequences:** the compressions yᵢ to Q satisfy Σᵢ PAF_{yᵢ} = 4q·δ, and the row sums satisfy Σ rᵢ² = 4q.
  - For q = 3 and even |K| the compressions are three single ±2 spikes in three blocks: 864 patterns, every one
    bordered (Kumar `border3.c`).
  - For q ≥ 5 the patterns and their borders are not yet enumerated (task 2).

**Status so far:**
- Only n = 669 was censused (`z3all`, 239 structures). All of them died at an exact level except the two S9 r = 90
  hosts on ℤ₂ × ℤ₃² × ℤ₃₇ and ℤ₂ × ℤ₉ × ℤ₃₇ (log₂ E_Q = 73.1; they pass G/N = ℤ₆; their ℤ₉ and ℤ₁₈ levels exceed the
  list caps).
- Kumar's scan covered only cyclic G = ℤ_m × GF(r) with m ≤ 30, t ≤ 30 and at most four to six spikes. It found a
  feasible stage 1 only at 537 among the open orders up to n = 1423 (and 1689, later settled), and it did not cover
  about 4,800 structures with more than 22 orbits.
- So non-cyclic G, twisted or coupled symmetry groups S, larger q, and n > 1423 are untried.

**The 19 target orders** (4n open; `data/open_orders_lt10000.json`):

| 4n | n | q |
|---|---|---|
| 2676 | 669 = 3·223 | 3 |
| 5028 | 1257 = 3·419 | 3 |
| 5748 | 1437 = 3·479 | 3 |
| 5892 | 1473 = 3·491 | 3 |
| 6852 | 1713 = 3·571 | 3 |
| 7572 | 1893 = 3·631 | 3 |
| 7716 | 1929 = 3·643 | 3 |
| 4460 | 1115 = 5·223 | 5 |
| 7180 | 1795 = 5·359 | 5 |
| 8380 | 2095 = 5·419 | 5 |
| 8860 | 2215 = 5·443 | 5 |
| 9260 | 2315 = 5·463 | 5 |
| 4676 | 1169 = 7·167 | 7 |
| 7924 | 1981 = 7·283 | 7 |
| 9716 | 2429 = 7·347 | 7 |
| 6876 | 1719 = 9·191 | 3 or 9 |
| 7876 | 1969 = 11·179 | 11 |
| 8684 | 2171 = 13·167 | 13 |
| 9308 | 2327 = 13·179 | 13 |

For each, |G| = n − q, so enumerate the abelian G of that order and the index-q subgroups K. At 669 the run rejected
the cyclic shapes but kept rich hosts on ℤ₂ × ℤ₃² × ℤ₃₇.

## 2. Tasks, in order (each is short; stop after any task if the time box is reached)

1. **Generalize the `z3all` family to ℤ_q** in `theory/bt_census.py` / `bt_engine.py`: targets −4q on K \ 0 and 0 off
   K, t = q, border size 4q.
   - Reproduce the 669 `z3all` numbers exactly as a regression test.
   - Reproduce the 2148 structure (`theory/bt_control2148.py`, or Kumar's stage-1 class at 537) as a positive control.
2. **Q-level patterns and borders for q = 5, 7, 9, 11, 13.**
   - Enumerate the compression patterns, i.e. quadruples on ℤ_q (and on ℤ₃² for q = 9) with Σ PAF = 4q·δ and
     Σ r² = 4q, up to symmetry.
   - For each pattern, decide whether a border exists: X, R̂, Ĉ of order 4q. Use Kumar's `border3.c` or
     `border_probe.py` approach, or a CP-SAT or clique search.
   - A q with no bordered pattern closes the family at every order with that q [M]. Record this first; it may remove
     whole rows of the table at almost no cost.
3. **Census with the cheap exact gates first,** at the orders whose q survives task 2, in order of 4n:
   - **Structures:** G of order n − q (every 2-part and odd part the enumeration limit allows), K of index q, S
     generated by per-axis multipliers, automorphisms of the 2-part and coupled generators, with 12 ≤ r ≤ 120 orbits.
   - **Screens:** coset-Parseval (`theory/b599x4_cosetscreen.py`, adapted to T in `report/b3_push_answer.md` §3.3),
     Theorem α and Lemma R (`theory/bt_thm_alpha.py`, `bt_thm_rigid.py`), the exact Q-level join, then the first
     finer level (`theory/bt_gates.py`, `bt_level1.py`, `bt_ladder.py`).
   - **Then** the stratified count read at the true target (never at the mean).
   - **Commit after each order:** a table row per family with counts dead / alive / best alive (r, log₂ E_Q), and the
     cost of its next exact gate.
4. **Only if a structure is alive with log₂ E_Q ≥ 10, its next exact level is affordable (≤ 10 minutes) and a border
   exists for its pattern:** run its exact levels to the end, as 3964 was found (`theory/bt_chain.py`,
   `bt_gates_levels.py`).
   - Use Kumar's ℤ₂ U/W split (his `search/bgs_uw.py`, the 2148 route) where G has an S-fixed involution, calibrated on
     2148 and its n = 177 analogue (26,095,424 solutions exactly in its stage-1 class).
   - Any core: assemble the border, verify independently, commit at once, and write a bundle `pub/H<4n>/` in the
     format of `pub/H3732/` with the border written out.
5. **Write `report/zq_census_answer.md`:** bottom line, per-q border results, per-order table, any construction, and
   the next computation for anything left alive.

## 3. What counts as done

- Every q in {5, 7, 9, 11, 13}: patterns enumerated and their borders decided (or the q left undecided with the
  reason).
- Every one of the 19 orders: either every structure in range is decided by the gates, or the alive ones are listed
  with their next gate and its cost.

## 4. Calibration warnings (from this repository)

- Estimates are not evidence of solutions. At the open orders every exactly decided structure with a positive
  estimate has so far been empty, up to 2^85 at 669.
- Use estimates only to order structures that pass every affordable exact gate.
- Lattice obstructions (congruences mod 2^k or mod |S|) are common. When S has a kernel of order m ≥ 9 on a coarse
  quotient, test the quotient equations modulo m first (Kumar `search/NOTES.md` §102).

## 5. Time box

- Stop computation after about 24 hours of wall-clock time, or earlier if the remaining work is the long tail of
  list-capped levels. Record the cutoff in the report, commit, and stop.
- Prefer finishing tasks 1–3 for the small orders (4n < 6000) over starting deep runs.
