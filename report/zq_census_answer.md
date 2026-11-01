# The ℤ_q all-character bordered family at the 19 composite open orders: patterns, borders, census, and no new matrix

Date: 2026-11-01/02 (one session; time box 24 h, cutoff recorded in §6). Branch `claude/focused-planck-j1645f`.
Answers `report/prompt_zq_census.md` (the prompt of commit 4064609). Scripts: `theory/zq_*.py` (new) and the
`theory/bt_*` engine. Logs: `logs/zq/`.

**Tags.** [P] proved; [M] exact computation (script named); [E] estimate (stratified Q-level count, read at the true
target, used only to rank).

**Credit.** The construction is Kumar's (`search/NOTES.md` §18–19 and `hadamard/2148` in
https://github.com/vatskishan03/hadamard, Apache-2.0): four ±1 sequences on G of order n − q with
Σᵢ PAFᵢ = −4q on K \ 0 and 0 off K, K of index q, T = all characters of Q = G/K, and a border of size 4q. The border
reduction follows his `border3.c`. His `bgs.stage1` was run read-only as an independent check (§4.2).

---

## 0. Bottom line

<!-- BOTTOM -->
(filled at the cutoff)
<!-- /BOTTOM -->

---

## 1. Tasks 1–3: the generalised family, the regression, the control [M]

- **Family.** `theory/bt_census2.py` now has the families `z5all`, `z7all`, `z9all`, `z11all`, `z13all` (Q = ℤ_q) and
  `z9x3all` (Q = ℤ₃²), next to `z3all`. The census driver for this prompt is `theory/zq_census.py`. Per structure it
  runs, cheapest first:
  - coset-Parseval for T [P];
  - Lemma R [P] (`bt_thm_rigid.rigid`);
  - the Q-level gate [M] (every row pattern with Σ r² = 4q);
  - the Q-level stratified count [E].

  No prefilter is applied, so every structure in 12 ≤ r ≤ 120 is decided by these gates or passes them.
  - The target is Γ = 4q·δ on Q, i.e. −4q on K \ 0 and 0 off K. It is built by the engine from T = all characters
    (`gamma = 'all'`).
  - K is chosen as in z3all: one index-q subgroup per cyclic factor of the odd part divisible by q. For Q = ℤ₃² the
    two ℤ₃ factors are dropped.
  - The odd parts are R.odd_options: cyclic factors, and ℤ_p² with GL(2, p) for p ≤ 5. So ℤ₃ × ℤ₉ (at 1893) is not
    enumerated; this is the enumeration limit.
- **Regression.** `theory/zq_regress669.py`, `logs/zq/census/regress_z3all_669.result`. All 239 records of
  `logs/bt/census2/z3all_669.jsonl` are reproduced exactly: key, verdict, r, e, Q-level solutions, Q-level patterns
  and log₂ E_Q.
  - The generalised enumeration gives 285 structures for 8 ≤ r ≤ 140, the same set as the code of commit 4064609.
  - The old log lacks the 46 structures with r = 104–130. One of them is S6 r = 118, with log₂ E_Q = 111.4 (§3).
  - The exact levels of the S9 r = 90 hosts are reproduced too: G/N = ℤ₆ passes with 130,553,856 weighted
    solutions, the 1.3·10⁸ of `report/b3_push_answer.md`.
- **Positive control.** Two runs on Kumar's 2148 structure:
  - `theory/bt_control2148.py` on his certificate (`logs/zq/control_2148.log`): the core, the border from the
    Paley-H₁₂ frame and the assembled H(2148) are exact.
  - The census at n = 537 (`logs/zq/census/z3all_537.jsonl`) contains his structure: ℤ₂ × ℤ₃ × ℤ₈₉, S of order 11,
    r = 54, e = 27. It passes the Q-level gate with log₂ E_Q = 62.0 and G/N = ℤ₆ with 248,832 weighted solutions. Its
    ℤ₁₇₈ level exceeds the list cap (1.4·10⁹).
  - The CP-SAT level gate of §3 finds its ℤ₆ level feasible. So neither gate kills a structure that has a core.

## 2. Task 4: compression patterns and borders

### 2.1 Patterns [M]

For n = q·p with p prime, |K| = p − 1 is even. So a core's compressions to Q are even, y = 2z, and the core condition
becomes Σᵢ PAF_{zᵢ} = q·δ on Q, with Σ|zᵢ|² = Σ rᵢ² = q (rᵢ the row sums of the z).
- `theory/zq_patterns.py` lists the single blocks by (weight, PAF). It then forms the multisets of four classes with
  weights summing to q and PAFs summing to q·δ, by meet in the middle.
- Each multiset is expanded to block vectors, and the quadruples are reduced up to: per-block translation, negation
  and reversal; block permutation; and a common automorphism of Q.

| q | Q | multisets of (weight, PAF) classes | patterns up to symmetry | log |
|---|---|---|---|---|
| 3 | ℤ₃ | 1 | 1 (three spikes) | `logs/zq/patterns_q3.log` |
| 5 | ℤ₅ | 5 | 3 | `patterns_q5.log` |
| 7 | ℤ₇ | 4 | 3 | `patterns_q7.log` |
| 9 | ℤ₉ | 136 | 67 | `patterns_q9.log` |
| 9 | ℤ₃² | 147 | 32 | `patterns_q9_33.log` |
| 11 | ℤ₁₁ | 756 | 176 | `patterns_q11.log` |
| 13 | ℤ₁₃ | 11,618 | not reduced | `patterns_q13_classes.log` |

For q = 13 the expansion to block vectors and the canonical reduction did not finish in 5.5 h, so the 11,618 class
multisets, each with one representative quadruple, are the q = 13 list. Every q has patterns, so no q is closed by an
empty pattern set.

### 2.2 Borders

With W = Mc/2, the GS array of z over Q (W Wᵀ = q I), the border condition X Ĉᵀ + R̂ Mcᵀ = 0 reads
R̂ = −X Ĉᵀ W/(2q) (Kumar's `border3.c` reduction):
- X and Ĉ must be Hadamard of order 4q;
- R̂ must be ±1;
- R̂ is then Hadamard automatically;
- X may be fixed in its equivalence class, since its left signed permutations go into R̂ and its right ones into Ĉ.

Solvers (all check X, R̂ and Ĉ Hadamard and X Ĉᵀ + R̂ Mcᵀ = 0 exactly):
- `theory/zq_border.py`: random Paley frames. It is too weak even at q = 3.
- `theory/zq_border_sat.py`: X Paley, Ĉ free, one CP-SAT model. It solves q = 3 in 0.6 s and is undecided at q = 5
  after 300 s.
- `theory/zq_border_gs.py`: X and Ĉ both GS arrays over Q (bilinear, XOR products). It solves q = 3 and q = 5
  patterns 0 and 1 when the solver is lucky; its symmetric (Williamson) variant is infeasible for all three q = 5
  patterns.
- `theory/zq_border_fixx.py`: X = the GS array of a perfect quadruple x, sampled uniformly by meet in the middle, and
  Ĉ = GS array of unknown c. For fixed x the condition is linear in c, a CP-SAT model of 4q bits; about 0.5 s per x
  at q = 5.

Results (`logs/zq/border_*.jsonl`):

<!-- ZQ-BORDERS -->
(filled at the cutoff)
<!-- /ZQ-BORDERS -->

No q is closed by the border: q = 3 and q = 5 have bordered patterns. For q ≥ 7 a border has not been found under the
sampled GS ansatz, and an exhaustive decision over all Hadamard X of order 4q ≥ 28 was out of reach. So the border
question for q = 7, 9, 11, 13 is **undecided**. It matters only for a structure that reaches a core, and none did (§3).

## 3. Task 5: census

**Gates, cheapest first.**
1. **Screens and Q-level** (`theory/zq_census.py`): coset-Parseval [P], Lemma R [P], the Q-level gate [M], the
   Q-level count [E].
2. **Exact finer levels** (`theory/bt_gates_levels.py --levels auto`):
   - every S-invariant N ≤ K, coarsest first;
   - 10 s per level, then 60 s for the levels capped at 10 s;
   - lists capped at 2·10⁷;
   - 3.5 GiB of memory;
   - runner `logs/zq/gates/run.sh` and `run_rev.sh`.
3. **CP-SAT level gate** (`theory/zq_xcheck_run.py`, new). It runs `theory/bt_xcheck_level.py` on every invariant
   level and every Q-level row pattern. That checker shares no code with `bt_*`.
   - Time limits: 20 s per level and pattern on the whole census (runner `logs/zq/xcheck/run.sh`), 30–300 s on single
     orders, and 600 s on the positive-estimate hosts (`run_pos.sh`).
   - A level with every pattern INFEASIBLE kills the structure [M].
   - This gate decides in about a second many levels that the list gates cannot reach, because their lists exceed
     the caps. Examples:
     - 1169 S41 at ℤ₁₄;
     - 2171 S41 at ℤ₂₆;
     - 2429 S43 at ℤ₁₄;
     - 1115 S12 r120 at ℤ₁₀.

The table counts each structure once, by its first killing gate. "alive" means it passed every gate run on it. The
last column gives the best alive structure's exact levels passed (N: weighted solutions, or CP-SAT feasible) and its
open levels with what stopped them:
- "list x": the refined list size over the 2·10⁷ cap;
- "10 s", "60 s": the time limit;
- "memory": refused under 3.5 GiB;
- "CP-SAT t s": undecided at that limit.

<!-- ZQ-TABLE -->
| 4n | n | family | structures | [P] coset-Parseval | [P] Lemma R | [M] Q-level empty | [M] finer level empty | undecided | alive | alive, log2 E_Q ≥ 0 | best alive: S, r, e, log2 E_Q | best alive: levels passed / open |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2148 | 537 | z3all | 13 | 0 | 1 | 0 | 0 | 0 | 12 | 4 | S11, r54, e27, 62.04 | N89: 248832, N89: CP-SAT feasible / N2: list 1.36e+09; CP-SAT 60 s |
| 2676 | 669 | z3all | 275 | 80 | 10 | 0 | 176 | 0 | 9 | 6 | S6, r118, e65, 111.41 | N111: 3825664, N74: 1373760 / N37: list 2.69e+08; CP-SAT 20 s, N6: list 8.12e+11; CP-SAT 20 s, N3: 10 s; list > cap; CP-SAT 20 s, N2: memory; CP-SAT 20 s |
| 4460 | 1115 | z5all | 87 | 0 | 6 | 0 | 78 | 0 | 3 | 1 | S18, r90, e47, 1.58 | N111: CP-SAT feasible / N74: list 3.06e+10; CP-SAT 60 s, N37: list 1e+17; CP-SAT 60 s, N6: 10 s; 60 s; CP-SAT 60 s, N3: list 9.5e+10; CP-SAT 60 s, N2: list 6.88e+14; CP-SAT 60 s |
| 4676 | 1169 | z7all | 8 | 0 | 3 | 3 | 2 | 0 | 0 | 0 | — | — |
| 5028 | 1257 | z3all | 47 | 0 | 16 | 0 | 29 | 0 | 2 | 0 | S45, r54, e27, -31.43 | N38: 20487168 / N209: 10 s; 60 s; CP-SAT 20 s, N22: 10 s; 60 s; CP-SAT 20 s, N19: list 4.39e+08; CP-SAT 20 s, N11: list 5.84e+08; CP-SAT 20 s, N2: memory; CP-SAT 20 s |
| 5748 | 1437 | z3all | 11 | 0 | 1 | 0 | 10 | 0 | 0 | 0 | — | — |
| 5892 | 1473 | z3all | 40 | 35 | 2 | 0 | 3 | 0 | 0 | 0 | — | — |
| 6852 | 1713 | z3all | 513 | 178 | 66 | 0 | 260 | 0 | 9 | 0 | S36, r90, e47, -108.76 | N114: 10194336 / N285: list 2.34e+07; CP-SAT 20 s, N190: 10 s; 60 s; CP-SAT 20 s, N95: list 1.53e+13; CP-SAT 20 s, N57: list 1.47e+10; CP-SAT 20 s, N38: list 1.05e+10; CP-SAT 20 s, N30: 10 s; capped (level join p; CP-SAT 20 s, N19: list 3.91e+16; CP-SAT 20 s, N15: list 5.48e+14; CP-SAT 20 s, N10: list 5.64e+10; CP-SAT 20 s, N6: list 6.26e+07; CP-SAT 20 s, N5: memory; CP-SAT 20 s, N3: list 5.48e+14; CP-SAT 20 s, N2: list 1.71e+14; CP-SAT 20 s |
| 6876 | 1719 | z3all | 496 | 30 | 83 | 9 | 10 | 0 | 364 | 0 | S30, r105, e54, -31.96 | N286: 3244032 / N143: memory, N52: 10 s, N44: list 8.09e+11, N26: memory, N22: list > cap, N13: list > cap, N11: list > cap, N4: list > cap, N2: memory |
| 6876 | 1719 | z9all | 109 | 39 | 18 | 0 | 0 | 0 | 52 | 0 | S18, r120, e71, None | — / — |
| 6876 | 1719 | z9x3all | 1910 | 556 | 182 | 0 | 0 | 0 | 1172 | 0 | S36, r84, e47, None | — / — |
| 7180 | 1795 | z5all | 6 | 0 | 2 | 3 | 1 | 0 | 0 | 0 | — | — |
| 7572 | 1893 | z3all | 68 | 64 | 0 | 0 | 0 | 0 | 4 | 0 | S24, r116, e87, None | — / — |
| 7716 | 1929 | z3all | 60 | 20 | 10 | 0 | 0 | 0 | 30 | 0 | S53, r54, e27, -13.37 | — / — |
| 7876 | 1969 | z11all | 23 | 0 | 3 | 17 | 2 | 0 | 1 | 1 | S22, r110, e59, 23.73 | — / N89: CP-SAT 30 s; CP-SAT 300 s; CP-SAT 600 s, N2: CP-SAT 30 s; CP-SAT 300 s; CP-SAT 600 s |
| 7924 | 1981 | z7all | 29 | 0 | 8 | 0 | 0 | 0 | 21 | 0 | S69, r54, e27, -36.95 | — / — |
| 8380 | 2095 | z5all | 67 | 0 | 21 | 8 | 0 | 0 | 38 | 0 | S45, r90, e45, -68.86 | — / — |
| 8684 | 2171 | z13all | 14 | 0 | 5 | 7 | 2 | 0 | 0 | 0 | — | — |
| 8860 | 2215 | z5all | 223 | 0 | 28 | 142 | 0 | 0 | 53 | 0 | S48, r90, e53, -115.09 | — / — |
| 9260 | 2315 | z5all | 85 | 0 | 21 | 0 | 55 | 0 | 9 | 0 | S60, r100, e59, -282.85 | — / N231: list 6e+11, N154: list 8.78e+07, N77: list 1.9e+16, N66: 10 s; 60 s, N42: 10 s; 60 s, N33: list 1.04e+15, N22: list 2.85e+08, N21: list 7.14e+13, N14: list 2.85e+08, N11: list 1.9e+16, N7: list 1.9e+16, N6: list 2.53e+09, N3: list 1.13e+19, N2: list 5.95e+11 |
| 9308 | 2327 | z13all | 48 | 0 | 5 | 38 | 4 | 0 | 1 | 1 | S33, r90, e45, 34.95 | — / N89: 10 s; 60 s; CP-SAT 60 s; CP-SAT 600 s, N2: list 1.66e+15; CP-SAT 60 s; CP-SAT 600 s |
| 9716 | 2429 | z7all | 14 | 0 | 3 | 7 | 4 | 0 | 0 | 0 | — | — |
<!-- /ZQ-TABLE -->

## 4. Notes per order

<!-- ZQ-NOTES -->
(filled at the cutoff)
<!-- /ZQ-NOTES -->

## 5. Not done, and the next computation

<!-- ZQ-NEXT -->
(filled at the cutoff)
<!-- /ZQ-NEXT -->

## 6. Cutoff and reproduction

<!-- ZQ-CUTOFF -->
(filled at the cutoff)
<!-- /ZQ-CUTOFF -->

```
python3 theory/zq_regress669.py logs/bt/census2/z3all_669.jsonl logs/zq/census/regress_z3all_669.jsonl
python3 theory/bt_control2148.py /path/to/vatskishan03/hadamard/hadamard/2148/certificate.json
python3 theory/zq_patterns.py 7 --out logs/zq/patterns_q7.json          # q = 3..11; 13 with --classes-only
python3 theory/zq_border_fixx.py logs/zq/patterns_q5.json --tmax 1200   # borders
logs/zq/census/run.sh 0/1                                               # census, all orders
logs/zq/gates/run.sh; logs/zq/gates/run_rev.sh                          # exact levels
logs/zq/xcheck/run.sh; logs/zq/xcheck/run_pos.sh                        # CP-SAT levels
python3 theory/zq_summary.py --fill report/zq_census_answer.md          # the table
```
