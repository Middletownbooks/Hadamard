# Prompt: Alpöge-style builders (b = 0 and b ≥ 3) for the remaining open orders, with the staged method (self-contained, maximum reasoning)

*For a maximum-reasoning Claude Code session on `middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`.*

- **Tags:** [M] verified by exact computation in the repo (name the script), [P] proved (give the proof),
  [E] Gaussian estimate (name the estimator and its calibration), [T] theoretical argument, [L] literature recollection
  (check it; arXiv and OEIS are blocked from these containers).
- **Branch:** other sessions push to it. Fetch and merge before every push. Never kill processes you did not start,
  and never `pkill -f` with a pattern that matches your own shell; kill by PID.

## 0. Goal and deliverables

Construct a Hadamard matrix H(4v) for some still-open v, using the array types of Alpöge, Voinov, Reynolds-Haertle and
Claude (2026, every open order below 2000):
- **`g`:** the Goethals–Seidel array, b = 0;
- **`u`, `v`, `w`:** group-developed GS with an inflated border, b ≥ 3.

All three of our recent constructions used the b = 1 bordered GS array. Your job is to find out whether the method that
built them transfers to b = 0 and b ≥ 3.

**Open v:**
- below 3000 (4v < 3000): 537, 571, 599, 643, 659, 669, 719, 739;
- 3000–4000: 751, 789, 823, 839, 859, 863, 883, 907, 917, 919, 933, 947, 955, 971, 991.

That second list was compiled by the user from public sources. Treat any order as open unless the repo has built it.

Deliverables, in order of value:
- **(A)** A verified H(4v): H Hᵀ = 4v·I exactly, checked by an independent script that shares no code with the search
  (pattern: `scripts/verify_h2092.py`), plus a bundle `pub/H{4v}/` in the format of `pub/H2092/` (with the b ≥ 3
  border and corner written out explicitly).
- **(B)** A method costed per solution at ≤ about 2^45 operations (≤ about 100 core-hours on 4 cores) for a specific
  (v, b, G, slice). It must be validated end to end on a smaller order of the same shape, with the exact first
  computation to start it and, if within budget, the run itself.
- **(C)** Exact obstructions (proofs or INFEASIBLE certificates) that remove configurations.
- **(D)** A ranked table of every (v, builder, configuration) you examined, with exact status and next computation.

Write to `report/alpoge_design_answer.md`, with scripts in `theory/ab_*.py` and logs in `logs/ab/`. Commit and push as
results arrive; containers restart. Commit trailers follow the repository's existing commits; no pull requests.

## 1. Where the project stands

**Built here, each verified three ways:**

| order | v | construction | record |
|---|---|---|---|
| H(2092) | 523 | b = 1, Z₁₈ × quartic classes on Z₂₉ | `report/paley_523_answer.md` |
| H(2292) | 573 | b = 1, Z₅₂ × Paley(11) with an order-3 multiplier | `report/RESULT_2292.md` |
| H(2588) | 647 | b = 1, the doubly staged ×353 slice on Z₂ × Z₁₇ × Z₁₉ | `report/paley_647_answer.md` |
| H(2884) | 721 = 7·103 | Cooper–Wallis: Schneider's T-matrices of order 103 with Williamson matrices of order 7 | `cw/RESULT_2884.md` |

**What made the b = 1 cases work:**
1. **A multiplier slice with an exact two-stage split.**
   - Fix a factor Z_F pointwise and act on a prime factor Z_q by an index-s cyclotomic subgroup ("generalized Paley").
   - Compressing along Z_q gives a small, exactly solvable stage 1 on Z_F: the letters are y = α + f·σ.
   - Stage 2 is a family of cyclotomic norm equations in the Gauss-period field.
   - Doubly staged slices act on two primes and need a joint stage 1 plus a cross level.
2. **Exact screens before any estimate.** The whole-slice Gaussian was wrong more than 15 times, by 30–100 bits, on
   exactly empty slices. The screens:
   - lag 0, rows and the Z₂ compression;
   - the α-congruence Σ_k PAF(β_k) ≡ target (mod 4f);
   - the compressed stage 1 by CP-SAT or exact S/D enumeration;
   - the joint test for doubly staged slices.
3. **A hierarchical stage 2.**
   - The Z₂ split gives an S level and a D level.
   - The S level uses a two-level Wagner tree with a Z₃ middle value; the D level a compiled D-check.
   - Gaussians are calibrated against exact counts. For cyclic stage 2 they were within 1–2 bits (523, 573).
4. **Production discipline.** Restart-safe workers with checkpoints, a tally of expected cores, and patience. H(2092)
   arrived at about 3 expected cores (a 5% wait).

## 2. What is known about the Alpöge et al. builders

Read `report/alpoge_methods_forensics.md`, `report/bge3_construction.md`, `report/bge3_answer.md` and
`report/unattempted_2000_3000.md` in full.

### 2.1 b = 0 (`g`)
- Four ±1 circulants of odd length v with Σ_k PAF_k = 4v·δ; rows satisfy Σ r_k² = 4v.
- Their records (v = 223, 283, 311, 487, 491) use multiplier orbits with d = 3–9 on prime v and individually flat
  blocks.
- We reproduced the method at v = 487 (flat-pool k-tree, n′ ≈ 83 bits).
- At 523 / d = 9 (n = 142, flat n′ ≈ 91) about 0.3 expected hits were accumulated with none found. That route is now
  moot, since H(2092) exists.
- **Row-sum lattice lemma** (`report/theory_arm_2000_3000.md`): high-symmetry b = 0 slices die on row sums at
  d ≳ 25–40 for v ≈ 520–740.
- **Twisted multipliers:** 571 with d = 15 and a 3-cycle on (B, C, D) is feasible (n ≈ 95 [E]), but it collapses to a
  lopsided 2-list problem (§2 of the unattempted report).

### 2.2 b ≥ 3 (`u`, `v`, `w`)

Their records: v = 347 (b = 5), 359 (b = 7), 479 (b = 3). The architecture [M] (`report/bge3_construction.md`):
- c = b + 1 is even, c | v + 1, h = (v + 1)/c − 1, G is abelian of order m = v − b = c·h, and H ≤ G has index c.
- **Core:** a GS quadruple over G with Σ_k PAF_k = −4b on H \ {0} and +4 off H.
- **Border:** an H(4c) with the 4 rows H4 ⊗ 1_c removed, inflated by h. The corner is a forced H(4b):
  C = −K·Cbᵀ·Kcᵀ/(4c), and it must come out ±1.
  - c = 4: every compression has a border.
  - c = 6: borders are known only for 347's placement.
  - c ≥ 10: no border is known. Settling this is part of the job.
- **The records sit at about 2^71–2^78** on the list-based cost model [E]. That model gives a barrier
  T ≳ 2^max(n/3, (n − log₂E)/2) for any method that treats the four blocks as list elements (`report/bge3_answer.md`
  §2.1).
- **Exact high-symmetry slices** (182 of them) for 571/c22, 599/c20, 647/c12, 659/c6, 739/c10 and 719/c8, c10 are all
  empty.
- The norm screen `bge3_normscreen.py` and the compression DP `bge3_cost2.py` killed most cheap configurations up to
  300 bits.
- **Forensics:**
  - all three published cores use h = h₁h₂ with coprime h₁, h₂ > 2;
  - 359's core (G = Z₂² × Z₈ × Z₁₁, H = Z₄₄) is almost exactly Legendre(11)-twisted on Z₁₁: 116 of 128 orbit pairs.
    That is the b = 7 analogue of the Paley(11) + multiplier slice that built H(2292).

### 2.3 Other signals
- **Half-shift template (167, 179):** dead for open orders (log₂E ≈ −13 per half at 523) [E].
- **479's GL(2,2) twist:** one parameter per block in the compression.

## 3. The central question

**The b ≥ 3 and b = 0 verdicts above (October 6) predate the staged method.** They used orbit slices, whole-list joins
and whole-slice bit counts. The same style of estimate put the b = 1 routes at 172–213 bits for 523, 573 and 647, and
all three were built. So:

1. **b ≥ 3 has a built-in compression.**
   - Compressing along H gives a stage 1 on Q = G/H (order c), whose targets are fixed by the core condition:
     - at 0: 4m − 4b(h − 1);
     - at nonzero Q-lags: +4h.
   - Compressing along a prime factor Z_q of h, with cyclotomic classes on Z_q and the rest fixed pointwise, gives a
     generalized-Paley-type slice with a two-valued target.
   - Derive the stage-1 and stage-2 systems exactly [P], and port the screens: the α-congruence (with the b ≥ 3
     targets), the rows/Z₂ screen, compressed CP-SAT and the joint test.
   - Check them on the published 479, 347 and 359 cores. Their cores must satisfy every condition you derive [M].
2. **b = 0 at composite v.**
   - The open composite v are 537 = 3·179, 669 = 3·223, 789 = 3·263, 917 = 7·131, 933 = 3·311 and 955 = 5·191.
   - Z_v = Z_F × Z_q with F ∈ {3, 5, 7} carries generalized-Paley slices with an exact stage 1 on Z_F: target 4v at 0
     and 0 elsewhere, with Σ r_k² = 4v.
   - 537 and 933 have **no** b ≥ 3 option (v + 1 = 2·269 and 2·467), and every b = 1 slice of 537 is dead. b = 0 is
     the remaining array-type route for them, apart from T-matrix routes (`report/tmat_routes.md`).
   - Screen all of these exactly. Mind the row-sum lattice lemma and the α-congruence for b = 0.
3. **Stage 2 for any survivor.** Port the hierarchical solver: S/D on a Z₂ factor where one exists, Wagner middle
   levels, a compiled D-check. Cost it with calibrated estimators (`theory/nd_stage2_est.py`, `theory/nd_levels.py`
   adapted to the new targets).

## 4. Configuration table (computed: c even, 4 ≤ c ≤ 30, c | v + 1; h = (v + 1)/c − 1; b = c − 1; m = c·h)

| v | b = 0 over Z_v | b ≥ 3 options c (b): h |
|---|---|---|
| 537 | 3·179 | none |
| 571 | prime | 4 (3): 142 = 2·71; 22 (21): 25; 26 (25): 21 = 3·7 |
| 599 | prime | 4: 149; 6: 99 = 3²·11; 8: 74 = 2·37; 10: 59; 12: 49; 20: 29; 24: 24; 30: 19 |
| 643 | prime | 4: 160 = 2⁵·5; 14: 45; 28: 22 = 2·11 |
| 659 | prime | 4: 164 = 2²·41; 6: 109; 10: 65 = 5·13; 12: 54; 20: 32; 22: 29; 30: 21 |
| 669 | 3·223 | 10 (9): 66 = 2·3·11 |
| 719 | prime | 4: 179; 6: 119 = 7·17; 8: 89; 10: 71; 12: 59; **16 (15): 44 = 2²·11**; 18: 39; 20: 35; 24: 29; 30: 23 |
| 739 | prime | 4: 184 = 2³·23; 10: 73; 20: 36 |
| 751 | prime | **4: 187 = 11·17**; 8: 93 = 3·31; 16: 46 |
| 789 | 3·263 | 10: 78 = 2·3·13 |
| 823 | prime | 4: 205 = 5·41; 8: 102 = 2·3·17 |
| 839 | prime | **4: 209 = 11·19**; 6: 139; 8: 104; 10: 83; 12: 69; 14: 59; 20: 41; 24: 34; 28: 29; 30: 27 |
| 859 | prime | 4: 214 = 2·107; 10: 85 = 5·17; 20: 42 |
| 863 | prime | 4: 215 = 5·43; **6: 143 = 11·13**; 8: 107; 12: 71; 16: 53; 18: 47; 24: 35 |
| 883 | prime | **4: 220 = 2²·5·11**; 26: 33 = 3·11 |
| 907 | prime | 4: 226 = 2·113 |
| 917 | 7·131 | 6: 152 = 2³·19; 18: 50 |
| 919 | prime | 4: 229; 8: 114 = 2·3·19; 10: 91 = 7·13; 20: 45 |
| 933 | 3·311 | none |
| 947 | prime | 4: 236 = 2²·59; 6: 157; 12: 78 |
| 955 | 5·191 | 4: 238 = 2·7·17 |
| 971 | prime | **4: 242 = 2·11²**; 6: 161 = 7·23; 12: 80; 18: 53 |
| 991 | prime | 4: 247 = 13·19; 8: 123 = 3·41; 16: 61 |

- G need not be cyclic: Q can be Z_c or a noncyclic group of order c (479 used Z₂², 359 used Z₂² × Z₈ with a
  non-split H).
- **Bold rows have 11 | h:** a Paley(11) or Legendre(11) slice on H is possible. That was 359's structure and our
  H(2292)'s.
- **c = 4 (b = 3)** always has a border and gives the largest h, so the most room for cyclotomic slices on h's prime
  factors.

## 5. Tools (read before writing new code)

| file | what it does |
|---|---|
| `forensics/*` | the authors' matrices and records; builders `gbuilder.py`, `hbuilder.py`; census, half-shift and twist analyses |
| `scripts/bge3_border.py`, `bge3_skeleton*.py`, `bge3_cost2.py`, `bge3_normscreen.py`, `bge3_screen2.py`; `report/bge3_answer.md` §8 | the b ≥ 3 border, skeleton, cost model and screens; positive controls (the v = 47 rung, `posctrl/`) |
| `data/orig_bge3/` | the 479, 347 and 359 cores in group coordinates |
| `theory/paley_wscreen.py` (`bge3` mode), `scripts/paley_bge3.py` | Paley-coset screens with b ≥ 3 targets (T_J = 4p off H′, −4bp on H′ \ 0) |
| `theory/gp_z2screen.py`, `gp_alpha_screen.py`, `slice_compress_cpsat.py` | generalized-Paley screens (b = 1 targets; generalize them) |
| `theory/nd_s1_sd.py`, `nd_joint_exhaust.py`, `nd_cross.py`, `nd_ejoin.py`, `nd_lag0_screen.py`, `nd_stage2_est.py`, `nd_levels.py` | exact compressed stage 1 (α-congruence + S/D), the joint test, cross level, energy-split joins, calibrated estimators (`report/next_design_answer.md` §9) |
| `theory/q523_*.py`, `theory/p647_*.py` | the H(2092) and H(2588) production pipelines (Wagner tree, D-check, restart-safe production). Read only; copy to `theory/ab_*` to adapt. |
| `scripts/kt_*.py`, `kaggle/` | the flat-pool k-tree for b = 0 (v = 487 control) |
| `pub/H2092/verify.py`, `scripts/verify_h2092_rebuild.py` | verification patterns |

## 6. Suggested order (reorder if you find better)

1. **Derive** the exact staged systems for b ≥ 3 GP-type slices and for b = 0 GP slices [P].
   - Verify them on the published 479, 347 and 359 cores, and on the H(892) b = 0 core (`data/hadamard_892.csv.gz`,
     Z₂₂₃, d = 3), which must satisfy them under their own symmetry.
2. **Screen exhaustively.**
   - Every b = 0 GP slice of the six composite v.
   - Every b ≥ 3 GP-type slice for c = 4 and c = 6, for all open v, over cyclic and noncyclic Q.
   - Then the 11 | h cases with larger c, once border existence is settled.
   - Record exact status for each. Whole-slice Gaussians decide nothing here.
3. **Settle the border question for c ≥ 8–10** (existence, or a construction from known H(4c)), before any large-c
   search.
4. **Positive control.** Reproduce a solution of a small same-shape b ≥ 3 instance (and b = 0 composite, e.g.
   v = 3q or 5q) end to end with the staged pipeline, to a verified matrix. Ideally reproduce 479's or 347's own slice
   in ≲ 2^40, as `report/bge3_answer.md` §7 recommends.
5. **Stage 2 and cost** for the best survivor, calibrated against exact counts. Run it if within budget.
6. If nothing survives, give the cleanest obstructions and the most promising unexplored family, with its first
   computation.

## 7. Rules

- A construction counts only if verified exactly. Separate theorems, proofs and estimates.
- Calibrate every Gaussian against exact counts. Lattice and 2-adic obstructions are real: a slice predicted at 2^9
  had exactly 0 because the target was unreachable mod 32.
- **Compute:** 4 cores and about 15 GB of RAM. Install `ortools`, `numba` and `sympy` with pip if missing.
- Write restart-safe runners. Arm a check-in so the session never idles while workers run (an idle session was
  reclaimed once and killed production workers). Commit results as they arrive.
