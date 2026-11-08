# New skew-Hadamard orders 4v from the Đoković–Seberry census: coverage, exact screen, calibrated search

Date: 2026-11-07/08. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_skew_census.md`.
Scripts: `theory/skew_*.py`, `theory/skew_*.c`, `scripts/verify_skew.py`. Logs: `logs/skew/`. Tables: `logs/skew/report_tables.md`
(written by `theory/skew_report_tables.py` from the logs).

**Tags.** [P] proved; [M] exact computation (script named); [E] estimate (estimator named). A cap or a timeout is never
a zero.

**Credit.** The low-energy four-list k-tree, the split-value ranking AB·CD/(AB + CD + |A| + |C|) and the reversal-twin
reduction are K. Kumar's (`search/cp_lowe.py`, `search/cp_enum.c`, `search/kt4.c`, `search/listgen`, NOTES §20, §110,
§111; Apache-2.0, github.com/vatskishan03/hadamard at d885022). `kt4.c` was compiled and run unchanged;
`cp_lowe.py` and `cp_hunt.py` were run unchanged from the read-only clone. The list generators `theory/skew_enum_lowe*.c`
and `theory/skew_lg*.c` are re-implementations for arbitrary abelian groups and say so in their headers. The
compression test of the screen is the integer test of his `skew_compress.py` / `census_screen.py` (NOTES §10, §14),
re-implemented; the forced-lag gate L, the quotient gate and the exact stage-1 class weights are new here.

<!-- BOTTOM -->
## 0. Bottom line

1. **No new skew-Hadamard order was found.** No matrix of a census order 4v was built in production, so nothing was
   published and `pub/skew/` was not created.
2. **The open list is 630 values** (227 prime, 82 at most 600), recomputed against Kumar's d885022 and unchanged [M].
   1661 and 1963 were left out as possibly settled by Schumer (6644, 7852); arXiv and Zenodo were unreachable here.
3. **Coverage map** (`logs/skew/coverage.jsonl`): 430 values appear in Kumar's notes or shape tables, 200 never do.
   His leads (1951, 2269, 2437, 661, 883) were not touched.
4. **Screen of all 630 values** (24–130 orbits per block, cyclic and ℤ_p² × ℤ_m groups) [M + E]:
   - exact gates (row sums, a new forced-lag gate L, ℤ_d and ℤ_p² compressions) leave 286 values with a live slice;
   - a calibrated k-tree cost model (pair work per solution in lists of 3–5·10⁷) leaves **exactly two slices within
     reach of 4 cores**: 361 = 19² over ℤ₁₉ × ℤ₁₉ with the split torus ⟨diag(4, 7)⟩ (a group–subgroup pair none of
     Kumar's tools covers), and 325 = 5² · 13 with Singer-3 × ⟨3⟩, both at 10^12.3 pairs per solution, between Kumar's
     solved 1123 (10^11.4) and 761 (10^12.9);
   - the best-ranked slice of all, 549 = 9 · 61 (estimate 10^8.3 solutions), and its non-cyclic variant are exactly
     dead at the ℤ₃ compression [M] — positive estimate, exact zero, again;
   - every other slice costs 10^14.4 or more, including every open prime structure screened (best: 433, 1483).
5. **Calibration** [M]: four of Kumar's orders were rebuilt blind from his structures only — 1204, 1412, 3052 and 4492
   (the comparable-cost control) — and verified by an independent checker (H Hᵀ = 4v I, H + Hᵀ = 2I). The final
   pipeline (general-group low-energy lists, Kumar's kt4 join and split ranking, reversal-twin reduction) found 4492
   in 446 s of join, so it was used in production.
6. **Production** [heuristic; nothing excluded]: 325 — all 14,194 ranked split values, no solution (in-list estimate
   10^0.09); 361 — 10,000 of 33,600 ranked split values in the 12 h budget, no solution (in-list estimate 10^0.84,
   about 0.5 expected twin classes after Kumar's Metropolis correction).
7. **Next:** 361-torus with lists of 3·10⁷ or more per block (or on a larger machine), and the ℤ_p² torus slices above
   v = 999 (§7).
<!-- /BOTTOM -->

---

## 1. The open list, and what was checked for novelty

- **Recomputed** (`theory/skew_census_recompute.py`, `logs/skew/census_recompute.txt`) against Kumar's commit d885022
  (2026-11-06; `git log -1` of the clone, unchanged during the session): census 737; settled 107 (every `skew/<4v>`
  with v odd in the census, `hadamard/3732`, and Schumer's 919); **open 630, 227 prime, 82 at most 600**, identical to
  `data/skew_open_census.json`. The 14 orders of `skew/` that are 8v or 16v (1528, 1800, 2216, 2392, 4344, 4632, 4712,
  7272, 7656, 9424, 9776, 13072, 14352, 14384) are not counted.
- **Sources** (`logs/skew/novelty_sources.md`):
  - arXiv (EKLO arXiv:2610.27175v2) and Zenodo (Schumer's record) are **refused by this session's egress proxy**
    (CONNECT 403); neither could be read. Kumar's NOVELTY.md (read here) records EKLO v1 (19 October 2026; smallest open
    skew order 356) and Schumer's deposit with orders 3676, 5436, 6644 and 7852.
  - **6644 = 4 · 1661 and 7852 = 4 · 1963 are on the open list.** Whether Schumer's matrices of those orders are skew
    could not be checked; both values were excluded from every search here (treated as possibly settled).
  - One web search (2025–2026 skew constructions): Karoui, arXiv:2602.16089, skew order 1252 = 4 · 313 (313 is not a
    census value); Kharaghani–Zaitsev arXiv:2510.06489 (quaternary). Nothing else.

## 2. Coverage map (step 1)

`theory/skew_coverage.py` → `logs/skew/coverage.jsonl`: for each of the 630 values, every paragraph of Kumar's
`search/NOTES.md` (skew sections and every paragraph that mentions skew), `skew/NEGACYCLIC.md`, `skew/PALEY_SKEW.md`, the
`skewplug`, `negaskew` and `listgen` READMEs, and his Paley-product shape table `scan38/data/shapes38.json`, with a
keyword label per paragraph. It is a text map: a match is evidence of contact, not of a decision.

- 430 of the 630 values are mentioned somewhere; 200 never are (the smallest: 229, 459, 485, 527, 565, 569, 595, 709).
- 230 values have Paley-product shapes in his table; for 99 of them every shape is decided there.
- 158 values occur in a paragraph with an exact-negative keyword, 154 with a heuristic-negative one.
- **Leads and running targets (excluded here):** the spin primes 1951, 2269, 2437 (§14), and 661 and 883 (k-trees with
  30 h and 44 h budgets, §115). 1661 and 1963 were also excluded (§1).
- **What his tools do not cover** (from the notes): the plain Goethals–Seidel array over ℤ_p² × ℤ_m with an odd-order
  subgroup of GL(2, p) that is not a field multiplication (his `cp_field` handles only GF(p²)-cyclotomy, i.e. Singer-type
  subgroups; his §12 survey used a capped model and stopped at v ≤ 999); cyclic slices with more than about 46 orbits
  for the exact chain; prime fibres with e > 42 except 251 (`listgen`).

## 3. The screen (step 2)

### 3.1 Slices and gates

`theory/skew_enum.py`, `theory/skew_lib.py`, `theory/skew_screen.py` → `logs/skew/screen.jsonl` (all 630 values).

- **Slices.** G = ℤ_v with every subgroup of ℤ_v^* of odd order generated by at most two elements; and, for p² ∥ v,
  G = ℤ_p × ℤ_p × ℤ_{v/p²} with H generated by (M, u) or by (M, 1) and (1, u), M an odd-order representative of
  GL(2, p) (scalar, split torus, unipotent, Singer) and u an odd-order unit (deduplicated by orbit signature). Kept:
  24 ≤ orbits ≤ 130. Odd |H| gives −1 ∉ H, so a skew block exists unless an orbit is closed under negation [P].
- **Gate R** [P/M]: the skew block has row sum 1; 4v − 1 = b² + c² + d² with each sum attainable by an H-invariant
  block (exact subset-sum DP over orbit sizes).
- **Gate L** [P] (new): if some lag s has f(s) = #{g : g + s ∈ Hg} > v/2, every H-invariant block (the skew one too)
  has PAF(s) ≥ 2f(s) − v > 0, so Σ PAF(s) ≠ 0. Example: v = 153 = 9 · 17, H = ⟨52⟩ (52 ≡ 1 mod 51): f(51) = 102 > 76.5.
- **Gate Z_d** [M]: the compression to ℤ_d along a cyclic factor, for every divisor d of that factor (all primes,
  composite d ≤ 60): each block's compressed vector is a sum over H_d-orbit coordinates, the skew block's is
  (1, y, −y)-antisymmetric; the four lists are joined exactly on the lag-0 energy 4v and every lag class of
  ⟨H_d, −1⟩ on ℤ_d. Validated: Kumar's exact verdicts (201/⟨22⟩ dead at ℤ₃, 265/⟨16⟩ at ℤ₅, 301/⟨4⟩ at ℤ₇; 469/⟨9⟩
  passes with one combination at ℤ₇; 209/⟨26⟩ passes). For d = v/q it is the stage-1 level of the staged slice.
  `quotient_gate` does the same onto a product of axes (ℤ_p²); it agrees with the axis gate on 361 and passes on the
  slice of Kumar's 1100.
- **Caps.** 2·10⁵ compressed vectors per block and a 30 s alarm per gate in the screen; gates were not run for slices
  with more than 65 orbits (beyond any k-tree budget, §3.3). Every cap is recorded as undecided.
- **Estimate** [E]: the lattice Gaussian of the class vector of the sum of the four blocks (uniform samples, exact
  uniform sampling at the row sum), on the hyperplane Σ_c n_c P(c) = 0 with lattice 4ℤ. On Kumar's exact counts:
  209/⟨26⟩ 10^2.26 against 288; 469/⟨9⟩ 10^1.84 against 192 (his model: 3 and 0.2). Replacing the Gaussian's stage-1
  factor by the exact weighted stage-1 count (`stage1_ratio`) changes these by +0.3 and +0.6 decades (toward the truth)
  and leaves 763 unchanged (`logs/skew/` command in §8).

### 3.2 Totals

| values | count |
|---|---:|
| some slice alive | 286 |
| every slice dead (exact gates) | 179 |
| no slice with 24-130 orbits | 158 |
| skipped (Kumar lead / Schumer) | 7 |

| (slice, pattern) records | count |
|---|---:|
| alive | 4405 |
| dead Z_d | 612 |
| dead R | 425 |
| dead L | 257 |
| no slice in range | 158 |
| skipped: Kumar lead | 3 |
| skipped: Kumar running | 2 |
| skipped: Schumer 6644? | 1 |
| skipped: Schumer 7852? | 1 |

"dead Z_d" counts every exact compression kill (ℤ₃ 62, ℤ₅ 60, ℤ₇ 98, ℤ₁₃ 97, ℤ₁₉ 68, ℤ₃₁ 58, …; `logs/skew/screen.jsonl`). The 158 values with no slice in range are mostly primes whose odd |H| give fewer than 24 or more than 130 orbits.

### 3.3 The k-tree cost model and the ranked table

`skew_lib.inlist_estimate` [E]: block k's kept list is modelled as uniform in the whitened ellipsoid |W_k(P − c_k)|² ≤ τ_k
(c_k the conditional mean given the target, τ_k the noncentral χ² quantile of L_k/N_k); the estimate inside the lists is
the lattice Gaussian of the sum. Checked against lists actually built: 353 10^4.44 (lists 10^4.42), 1123 10^3.16 (10^3.24);
against Kumar's printed in-list values it is 0.5–0.6 decades higher (1123, 761, 613), consistently. The cost is
log₁₀(L_A L_B + L_C L_D) − log₁₀ E_in, the pair work per solution of one pass of the k-tree. On this scale the solved
controls are 613: 10.2, 1123: 11.4, 761: 12.9.

One row per (v, classes, row pattern); slices with C ≤ 30 and slice estimate ≥ 10³ (45 slices) plus the prime pass.

| # | v | slice | \|H\| | orbits | classes | rows | log10 E (slice) | log10 E in lists | cost | gates | next step |
|---:|---:|---|---:|---:|---:|---|---:|---:|---:|---|---|
| 1 | 549 | G=Z549, H=[[22]] | 15 | 41 | 20 | [1, 3, 31, 35] | 8.31 | 2.76 | 12.22 | dead at Z3 [M] (all patterns) | none (dead) |
| 2 | 361 | G=Z19xZ19, H=[('torus', [[4, 0], [0, 7]])] | 9 | 45 | 22 | [1, 5, 7, 37] | 11.86 | 2.8 | 12.25 | alive: both Z19 compressions pass [M] | searched (§5) |
| 3 | 325 | G=Z5xZ5xZ13, H=[('singer', [[1, 1], [2, 3]], 1), ('id', None | 9 | 45 | 22 | [1, 7, 17, 31] | 11.1 | 2.72 | 12.33 | alive: Z13 passes; Z5^2 undecided (cap) | searched (§5) |
| 4 | 433 | G=Z433, H=[[27]] | 9 | 49 | 24 | [1, 1, 19, 37] | 14.91 | 0.83 | 14.41 | prime: no Z_d gate | k-tree beyond budget (cost >= 14) |
| 5 | 1483 | G=Z1483, H=[[28]] | 39 | 39 | 19 | [1, 1, 1, 77] | 6.97 | 0.38 | 14.58 | prime: no Z_d gate | k-tree beyond budget (cost >= 14) |
| 6 | 589 | G=Z589, H=[[102]] | 15 | 45 | 22 | [1, 5, 31, 37] | 8.53 | 0.17 | 14.88 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 7 | 549 | G=Z3xZ3xZ61, H=[('torus', [[1, 0], [0, 1]], 12)] | 15 | 45 | 22 | [1, 3, 31, 35] | 5.98 | -1.5 | 16.55 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 8 | 487 | G=Z487, H=[[41]] | 9 | 55 | 27 | [1, 17, 17, 37] | 16.17 | -1.44 | 16.68 | prime: no Z_d gate | k-tree beyond budget (cost >= 14) |
| 9 | 817 | G=Z817, H=[[68]] | 21 | 45 | 22 | [1, 7, 37, 43] | 6.54 | -1.7 | 16.75 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 10 | 725 | G=Z5xZ5xZ29, H=[('singer', [[1, 1], [2, 3]], 7)] | 21 | 45 | 22 | [1, 21, 33, 37] | 4.46 | -1.79 | 16.84 | gates pass | k-tree beyond budget (cost >= 14) |
| 11 | 589 | G=Z589, H=[[64]] | 15 | 49 | 24 | [1, 5, 11, 47] | 8.17 | -1.95 | 17.19 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 12 | 489 | G=Z489, H=[[40]] | 9 | 57 | 28 | [1, 15, 19, 37] | 15.66 | -2.09 | 17.33 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 13 | 481 | G=Z481, H=[[9]] | 9 | 57 | 28 | [1, 1, 31, 31] | 14.83 | -2.35 | 17.59 | gates pass | k-tree beyond budget (cost >= 14) |
| 14 | 723 | G=Z723, H=[[94]] | 15 | 51 | 25 | [1, 29, 31, 33] | 11.82 | -2.45 | 17.69 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 15 | 481 | G=Z481, H=[[16]] | 9 | 57 | 28 | [1, 5, 23, 37] | 14.77 | -2.47 | 17.71 | gates pass | k-tree beyond budget (cost >= 14) |
| 16 | 949 | G=Z949, H=[[16], [74]] | 27 | 45 | 22 | [1, 19, 25, 53] | 3.89 | -3.04 | 18.09 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 17 | 523 | G=Z523, H=[[19]] | 9 | 59 | 29 | [1, 19, 19, 37] | 16.95 | -3.24 | 18.48 | prime: no Z_d gate | k-tree beyond budget (cost >= 14) |
| 18 | 597 | G=Z597, H=[[61]] | 11 | 57 | 28 | [1, 1, 19, 45] | 14.53 | -3.3 | 18.54 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 19 | 817 | G=Z817, H=[[11]] | 21 | 49 | 24 | [1, 7, 37, 43] | 6.15 | -3.65 | 18.89 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 20 | 655 | G=Z655, H=[[176]] | 13 | 55 | 27 | [1, 3, 3, 51] | 11.68 | -4.09 | 19.33 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 21 | 1565 | G=Z1565, H=[[16]] | 39 | 45 | 22 | [1, 3, 3, 79] | 4.58 | -4.41 | 19.45 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 22 | 1201 | G=Z1201, H=[[96]] | 25 | 49 | 24 | [1, 1, 49, 49] | 9.86 | -4.35 | 19.59 | prime: no Z_d gate | k-tree beyond budget (cost >= 14) |
| 23 | 1203 | G=Z1203, H=[[25]] | 25 | 51 | 25 | [1, 1, 47, 51] | 9.15 | -5.09 | 20.33 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 24 | 793 | G=Z793, H=[[22]] | 15 | 57 | 28 | [1, 29, 31, 37] | 11.43 | -5.58 | 20.82 | gates pass | k-tree beyond budget (cost >= 14) |
| 25 | 921 | G=Z921, H=[[64]] | 17 | 57 | 28 | [1, 33, 35, 37] | 11.93 | -5.82 | 21.06 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 26 | 361 | G=Z19xZ19, H=[('torus', [[1, 0], [0, 4]])] | 9 | 57 | 28 | [1, 5, 7, 37] | 4.56 | -6.3 | 21.54 | gates pass | k-tree beyond budget (cost >= 14) |
| 27 | 1301 | G=Z1301, H=[[44]] | 25 | 53 | 26 | [1, 1, 51, 51] | 10.25 | -6.33 | 21.57 | prime: no Z_d gate | k-tree beyond budget (cost >= 14) |
| 28 | 1103 | G=Z1103, H=[[17]] | 19 | 59 | 29 | [1, 37, 39, 39] | 12.72 | -7.37 | 22.61 | prime: no Z_d gate | k-tree beyond budget (cost >= 14) |
| 29 | 899 | G=Z899, H=[[25]] | 21 | 55 | 27 | [1, 15, 39, 43] | 4.58 | -7.59 | 22.83 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 30 | 1255 | G=Z1255, H=[[16]] | 25 | 55 | 27 | [1, 1, 47, 53] | 7.6 | -8.06 | 23.3 | gates undecided/none | k-tree beyond budget (cost >= 14) |

**Reading.**
- After the exact gates **only two slices are within budget: 361 = 19² over ℤ₁₉ × ℤ₁₉ with the split torus
  H = ⟨diag(4, 7)⟩ (order 9; 45 orbits, 22 lag classes) and 325 = 5² · 13 over ℤ₅ × ℤ₅ × ℤ₁₃ with H = ⟨Singer of order 3⟩ ×
  ⟨3⟩ (order 9; 45 orbits, 22 classes)**, both at cost ≈ 10^12.3, between 1123 and 761.
- The best-ranked slice of all, 549 = 9 · 61 with ⟨22⟩ or ⟨16⟩ (|H| = 15, 20 classes, estimate 10^8.3), and the
  non-cyclic 549 slice ℤ₃² × ℤ₆₁ with ⟨(u, 12)⟩, u unipotent, are **exactly dead at the ℤ₃ compression for every row
  pattern** [M] (`logs/skew/gates_candidates.jsonl`): positive estimate, exact zero, at the coarsest quotient — the
  calibration warning of the prompt, again.
- Everything else costs at least 10^14.4 (433, 1483, 589, …), 100 times 1123 or more: beyond 4 cores.
- 361's torus slice is in no tool of Kumar's: `cp_sweep.structures` skips q | m (ℤ₁₉ × GF(19)), and `cp_field` has only
  field multiplications. 325's Singer slice is a GF(25)-fibre structure that his `cp_hunt` can estimate (its sampled
  stage-1 class: 10^−0.03 for the class, 10^−1.97 in lists of 3·10⁷); it was never searched.

### 3.4 Kumar's own estimate pass (cp_hunt, nontrivial K_m)

`cp_hunt.py skew <401 composite values> --estimate-only --L 3e7` (unchanged, `logs/skew/cp_hunt_est.jsonl`): stopped after
30 structures (225 … 695) to free a core. No structure has a positive in-list estimate (best: 623, 10^−0.83; 325
GF(25), 10^−1.97), consistent with his notes; the run is incomplete and decides nothing.

### 3.5 Primes

`theory/skew_prime_pass.py` (`logs/skew/prime_pass.jsonl`): for the open primes, every odd |H| with 30 ≤ e ≤ 80 classes:
20 structures pass gate R. Their k-tree costs (`logs/skew/prime_inlist.jsonl`) are 14.4 (433, |H| = 9) and 14.6 (1483,
|H| = 39) at best, and 19.6 to 33 for the rest (1201, 1301, 1103, 1303, 859, 1741, 2113, 2381). No open prime is
within reach of this method on this machine.

## 4. Calibration (step 3)

All four controls were run blind: only the structure (group, H, row sums) was taken from Kumar; no signs.

| control | structure | tool | lists (A; B, C, D) | E in lists [E] | result | time |
|---|---|---|---|---|---|---|
| 1204 (v = 301) | ℤ₇ × GF(43), e = 6, H = ⟨127⟩ | Kumar's cp_lowe.py, seed 2 (his: seed 1) | 64,000; 1.0·10⁷ ×3 | 10^2.53 (his 10^2.52) | 1 solution at split value 292 | 534 s, 2 threads |
| 1412 (v = 353) | ℤ₃₅₃, H = ⟨22⟩, rows 1, 21, 21, 23 | this pipeline (exact enumeration) | 65,536; 10⁶ ×3 | 10^4.42 | 16 tuples in the first 200 split values; 4 matrices verified | ~1 min lists, 2 s join |
| 3052 (v = 763) | ℤ₇₆₃, H = ⟨9⟩, rows 1, 1, 5, 55 | this pipeline | 32,768; 2·10⁷, 1.5·10⁷, 2·10⁷ | 10^5.38 | 16 tuples at split value 27 (kt4 order); with dedup + ranking: 2 in the first 200 | ~1.5 min lists, 7–10 s join |
| 4492 (v = 1123) | ℤ₁₁₂₃, H = ⟨691⟩, rows 1, 1, 1, 67 | this pipeline, comparable cost (11.4) | 131,072; 4·10⁷, 2·10⁷, 2·10⁷ | 10^3.24 | first version: 16 tuples at split value 718 after 2.4·10¹¹ CD pairs (3354 s join); **final version** (twin reduction, sampled-pair ranking): 2 solutions within 1000 split values, 3.0·10¹⁰ CD pairs | 7 min lists + 446 s join, 2–3 threads |

- Every hit was rebuilt and checked by `scripts/verify_skew.py` (integer PAF on every lag, skew first block,
  invariance, H Hᵀ = 4v I, H + Hᵀ = 2I); hits in `logs/skew/calib/HITS.jsonl`. The verifier also passes Kumar's
  certificate of 1100 (ℤ₅ × ℤ₅ × ℤ₁₁) and fails a copy with one sign flipped.
- Kumar's runs of the same controls: 1204 ≈ 6 min; 1412 16 solutions in 80 s; 3052 found by the exact chain (his
  sampled k-tree had missed it); 4492 2.7 thread-hours with lists 1e8, 3e7, 3e7.
- **Kill rule.** The method reproduced a control of comparable estimated cost (1123, 11.4, against 12.3 for the
  candidates), in its final form, so it was used in production.
- **Bugs found by the controls and fixed before production** (each re-validated): the enumerator's output buffer
  (overflow → τ retry), memory (int16 keys, packed signs, chunked twin reduction; one container OOM kill), the
  Gaussian split ranking (unnormalised; it starved one side on 763) replaced by Kumar's sampled-pair ranking, an int64
  overflow in its split encoding, the Metropolis generator's saturation (stop rules) and speed (whitened incremental
  energy, `skew_lg2.c`).

## 5. Step 4: decisions and runs

### 5.1 Exact stage 1 and gates of the two candidates [M]

For each candidate the stage-1 level (the compressions to the quotients that H acts on) was decided exactly before
any list was built (`logs/skew/gates_candidates.jsonl`):
- **361, ℤ₁₉ × ℤ₁₉, H = ⟨diag(4, 7)⟩** (row patterns (1; 5, 7, 37), (1; 7, 13, 35), (1; 11, 19, 31), (1; 17, 23, 25)):
  the compression along the first axis (H acts by 4, order 9: 3 coordinates) has 43–45 PAF-classes of solutions per
  pattern, by two codes (`zp_gate` and the independent `quotient_gate`, equal counts); along the second axis (7, order
  3: 7 coordinates) 2.1–2.3·10⁶. The exact weighted stage-1 classes (`stage1_classes`): 1356 compression quadruples
  for (5, 7, 37), the heaviest holding 0.4 % of the weight. Staging along the first axis would shrink the blocks only
  3–6×, so the slice was searched unstaged. Gate L passes. Alive.
- **325, ℤ₅ × ℤ₅ × ℤ₁₃, H = ⟨Singer of order 3⟩ × ⟨3⟩** (5 patterns): the ℤ₁₃ compression passes (4.2–4.7·10⁴ classes);
  the ℤ₅² compression exceeds its cap of 2·10⁶ vectors per block (undecided, not a zero); Singer-3 has no invariant
  line, so there is no ℤ₅ compression. Alive as far as decided.
- Planted controls of the gates: the slice of Kumar's 1100 (ℤ₅² × ℤ₁₁, a known solution) passes both the ℤ₅² and the
  ℤ₁₁ gate; the known solutions of 209, 469 and 763 pass theirs.

### 5.2 Searches [heuristic: a finished search without a hit excludes nothing]

Both with `theory/skew_kt.py` (lists by the whitened-energy Metropolis generator `skew_lg2.c`, 2·10¹⁰ moves per block,
reversal-twin reduction, Kumar's split-value ranking, `kt4.c`), budget fixed in advance at 12 h per candidate.

| candidate | rows | lists (A; B, C, D) after twin reduction | E in lists [E] | split values | pairs | result |
|---|---|---|---|---|---|---|
| 325 Singer | (1, 7, 17, 31) | 4,194,304; 6,311,844, 6,222,074, 6,121,921 (from 1.2–1.3·10⁷) | 10^0.09 | all 14,194 seen on both sides | 7.3·10⁹ AB, 2.0·10¹¹ CD | **no solution** (52 min join, 2 threads; lists 1.5 h) |
| 361 torus | (1, 17, 23, 25) | 4,194,304; 9,157,446, 9,940,573, 10,026,076 (from 1.8–2.0·10⁷) | 10^0.84 | 10,000 of 33,600 (the top-ranked; stopped at the 12 h budget) | 1.3·10¹¹ AB, 2.9·10¹² CD | **no solution** (join 03:15–15:05 UTC: 2 threads to 05:00, then 4; lists 2.2 h) |

- The lists are below their targets (3·10⁷ for 361, 2·10⁷ for 325) because each block stopped at the cap of 2·10¹⁰
  Metropolis moves (about 30–45 min per block on 2 threads); exact enumeration of 7·10¹¹-member blocks was measured
  at about 7 h per block and abandoned. Kumar measured, for Metropolis lists at 251, a yield of 0.07 of the Gaussian
  count inside them; with that factor the expectations are about 0.1 (325) and 0.5 (361) twin classes in the whole
  candidate sets.
- A container restart at 02:42 UTC stopped both runs; they resumed from the cached lists. 361's first driver had been
  started with the earlier (Gaussian) split ranking and no twin reduction; it was stopped after its lists, before
  its join, and the join was rerun with the final code. A second restart at 15:16 UTC coincided with the budget stop
  (the last complete chunk ended at split value 10,000, 15:05 UTC); the run was not resumed.
- **Reading** [E]. The share of a solution's split distribution held by the first 10,000 of the 33,600 ranked values is
  not known exactly; on the controls the hits came in the first 1–3 % of the ranked values. With an expectation of
  order 0.5 twin classes in 361's whole lists, a miss is the likeliest single outcome and says nothing about the slice.
  Neither run excludes anything: the lists hold 10⁻⁴·⁵ of each block.


## 6. Matrices

None. No skew-Hadamard matrix of a census order was found in production. The verified matrices of this session are
calibration rebuilds of known orders (1204, 1412, 3052, 4492; `logs/skew/calib/HITS.jsonl`), not published:
`pub/skew/` was not created. `theory/skew_publish.py` is ready for a hit (bundle with certificate, sequences,
standalone `verify.py`, README).

## 7. Next computation

1. **361-torus, the remaining 23,600 split values** of the same lists (`skew_kt.py` resumes from
   `sigmas_done.txt` and the cached lists; about 12 more hours on 4 cores at the measured rate), or, better, **larger
   lists**: the Metropolis lists stopped at 2·10¹⁰ moves with 1.8–2.0·10⁷ members against a target of 3·10⁷; at
   3·10⁷ per block (after twin reduction 1.5·10⁷) the in-list estimate rises by about 10^0.5. On Kumar's 32-core
   machine (lists of 6·10⁸ as for 761) this slice is in the regime where 761 and 841 were solved.
2. **325-Singer with larger lists** (it had 1.2–1.3·10⁷ per block, estimate 10^0.09); and decide its ℤ₅² compression
   exactly (the gate needs more than 2·10⁶ compressed vectors per block: a C version of `quotient_gate`).
3. **The other split-torus and mixed slices of ℤ_p² × ℤ_m above v = 999** (1083, 1175, 1225, 1369, 1681, 1849, 2209, …):
   the screen ran them only with 24–130 orbits and gates for at most 65 orbits; none reached the cost table, but the
   in-list ranking used C ≤ 30 only. A ranking with C ≤ 36 and lists of 10⁸ is the cheap next screen.
4. **Exact gates for the undecided screen entries** (4405 alive records include every 'some gate undecided' and every
   slice above 65 orbits, whose gates were not run): a C join for the Z_d level would decide them; on this session's
   evidence (549, 153, 201, 265, 301) a large share would die there.

## 8. Reproduce

```
python3 theory/skew_census_recompute.py /home/user/vatskishan03/hadamard
python3 theory/skew_coverage.py /home/user/vatskishan03/hadamard
python3 theory/skew_screen.py 2                     # logs/skew/screen.jsonl (resumable; a few hours on 2 cores)
python3 theory/skew_inlist_rank.py 30 3             # logs/skew/inlist_rank.jsonl
python3 theory/skew_report_tables.py > logs/skew/report_tables.md
# k-tree (binaries: cc -O3 -march=native -fopenmp {skew_enum_lowe2,skew_lg2}.c -lm; kt4.c from Kumar's search/)
SKEW_WORK=<dir with bin/> python3 theory/skew_kt.py theory/specs/calib_1123.json --L 4e7,2e7,2e7 --nsig 100000 --threads 2 --K 1 --seed 3 --dedup
SKEW_WORK=<dir> python3 theory/skew_kt.py theory/specs/p361_torus_a.json --L 3e7,3e7,3e7 --enummax 0 --gensecs 5400 --T 0.6 --dedup --nsig 200000 --K 4
python3 scripts/verify_skew.py CERT.json
```
