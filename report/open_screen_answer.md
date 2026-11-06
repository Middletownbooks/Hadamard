# The 85 open orders below 10000: a ranked screen for the cheapest staged slice, and exact stage-1 decisions

Date: 2026-11-05/06. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_open_screen.md`.
Scripts: `theory/open_screen_*.py` (new). Logs: `logs/open_screen/`. Table: `logs/open_screen/table.jsonl`.

**Tags.** [P] proved; [M] exact computation (script named); [E] estimate (estimator, sampling and calibration
named). A cap or a timeout is never a zero.

**Credit.**
- The bordered k4 family is Kumar's form of [EKLO] Theorem 7; the ℤ_q all-character family is Kumar's (2148); both
  as in `report/b3_push_answer.md` and `report/zq_census_answer.md`.
- Theorem B (the k4 border for odd |K|) and the compression normal form are the OpenAI model's
  (`report/openai_4p3_review.md`).
- The slice machinery (`theory/ab_lib.py`, `ab_gauss.py`, `ab_sweep.py`, `ab_b0_stage1.py`), the lattice Gaussian
  (`theory/q523_gauss.py`), the deficit model (523, 933, 863), `scripts/b0_open_screen.py`, `scripts/gs_twist.py`,
  `scripts/tmat_primev.py` and `scripts/tmat_estimate.py` are the repository's earlier work, reused here.

---

## 0. Bottom line

<!-- BOTTOM -->
1. **No open order has a staged slice that is in the regime with a stage 1 shown non-empty.** No core and no
   Hadamard matrix was found, and none was sought (step 6).
2. **The estimates put 13 slices at 10 orders in the regime** (tier 0: D_def ≤ 26, h_D ≤ 40, E per y ≥ 2^10,
   C₂ ≤ 150, and log₂ E_tot ≥ 0) [E]. The orders are 3452, 3788, 4156, 4748, 5036, 5692, 5708, 5804, 7572 and 9260.
   - Earlier, 6876 (1719 b0) had a tier-0 slice. Its stage 1 is now exactly empty, so it is no longer tier 0.
   - The calibrated cost model reproduces the two production runs (933: 2^24.5 against 2^25.7 measured; 523: 2^44.3
     against 2^44.2) [E, §2.3].
   - The tier-0 costs run from 2^27 to 2^41 per core, all inside the 2^46 budget [E].
3. **Step 5: every candidate decided so far is exactly empty at stage 1 [M]** (§4). In cost order:
   - **6876 (v = 1719 = 9·191), b0, the H(3732) shape:** ℤ₉ × ℤ₁₉₁ (s = 10, f = 19, U = −1), estimated cost 2^42.8 and
     D_def 25.6, and both of its ℤ₃ × ℤ₃ variants. All three are empty by two methods each, with every planted
     control recovered.
   - **3788 (v = 947), b1:** ℤ₁₁ × ℤ₂ × ℤ₄₃, s = 6, U of order 5 (cost 2^26.0). Empty by two methods.
   - **5036 (v = 1259), b1:** ℤ₁₇ × ℤ₂ × ℤ₃₇, s = 4, U of order 4 (cost 2^24.8). Empty by two methods.
   - **5804 (v = 1451), b1:** ℤ₂₅ × ℤ₂ × ℤ₂₉, s = 4, U of order 10 (cost 2^27.0). Empty by one exhaustive method,
     the complete bucketed join over all 402,884 pruned blocks with planted 3/3. CP-SAT is UNKNOWN at 1 h, so this
     is a single-method decision.
   - **3452 (v = 863), k4:** ℤ₄ × ℤ₅ × ℤ₄₃, s = 6, U = −1 on ℤ₄ × ℤ₅ (cost 2^35.1; border by Theorem B):
     <!-- K863 -->see §4<!-- /K863 -->.
   - Earlier exact screens had already emptied the two cheapest-looking slices, the Paley slices of 3788
     (cost 2^17.2) and 2396 (2^24.7). Their estimates were 2^33 and 2^52 expected cores, so this is again "positive
     estimate, exact zero".
4. **Calibration warning, confirmed again.** Every slice decided exactly here had a positive estimate, from
   log₂ E_tot = 3.5 (1259) to 32 (1719), and every decided one is empty.
   - Several tier-0 estimates use the energy-only stage-1 fallback, which is optimistic (6.6 bits at its one check).
   - The remaining tier-0 slices are targets for exact stage 1, not for a solver (§6).
5. **The other families have no route within budget [P, M, E]** (§5):
   - direct b0 at prime v: the cheapest is n = 157.5 bits (1103, d = 19), against n = 111.9 for the slice that built
     H(3532);
   - Cooper–Wallis: the cheapest T-slice is T(463) with k = 11 (n = 101.5), the route to H(9260), at about 2^63
     operations by the T-route scale of `report/tmat_routes.md` §10.5;
   - ℤ_q all-character: no survivor beyond the census.
<!-- /BOTTOM -->

---

## 1. What was screened

<!-- SCOPE -->
- **Orders.** All 85 open orders of `data/open_orders_lt10000.json`, v = 599 … 2459.
- **Staged families, 3242 slices** (`logs/open_screen/screens_s*.jsonl`):

| family | orders with ≥ 1 slice | slices | A dead | B dead | ℤ₂ dead | w-screen dead | stage 1 empty [M] | alive | costed (sGeff ≤ 101) |
|---|---|---|---|---|---|---|---|---|---|
| b0 (composite v) | 14 of 19 | 145 | 95 | 2 | – | – | 3 | 45 | 24 |
| b1 | 68 of 85 | 1475 | 614 | 121 | 11 | 4 | 2 | 723 | 509 |
| k4 (v ≡ 3 mod 4) | 63 of 72 | 1622 | 1104 | 23 | – | – | – | 495 | 342 |

  - "alive" counts slices that pass every exact gate run on them: 1263 in all, of which 699 have B undecided (more
    than 4000 α-keys). A cap is never a zero.
  - Every slice has a decided screen A.
  - Key fix (2026-11-06): on 2-power axes of order ≥ 8 the units have several cyclic subgroups of the same order. The
    first keys recorded only (axis, order), so 114 keys stood for 290 slices. Every slice was screened; the screen
    rows were re-keyed from their recorded generators and the 65 affected cost rows recomputed.
- **The ℤ₂ and w-screens** (`theory/open_screen_z2gate.py`, `logs/open_screen/z2gate.jsonl`) apply to b1 slices with
  U = 1 (ℤ₂), and to b1 and k4 slices with s = 2 and U = 1 (w-screen). They reproduce the earlier kills of 947 s = 2
  (`logs/gp_screen/gp_z2screen_positive.jsonl`) and 599 s = 2 (Theorem 5, `report/paley_design_answer.md` §4).
- **Costs** were estimated for every alive slice with sGeff ≤ 101.
  - Rows: 873 alive slices costed (952 cost rows in `logs/open_screen/deficit_s*.jsonl`). One is a memory error and
    26 are 900 s timeouts, all undecided.
  - Why sGeff ≤ 101: the measured ratio C₂/sGeff is at least 1.48 for every costed slice with sGeff ≥ 60, so sGeff >
    101 means C₂ > 150, outside the regime [E].
  - Check: no costed slice with sGeff > 80 has C₂ ≤ 150 [E].
- **Other families:**
  - direct (unstaged) b0 multiplier slices at the 64 prime v (`scripts/b0_open_screen.py`);
  - Cooper–Wallis T(t) at the 15 prime cofactors t ≥ 167 of the composite v (`scripts/tmat_primev.py`);
  - the ℤ_q all-character family, from its census (`report/zq_census_answer.md` §3).
<!-- /SCOPE -->

## 2. Method

### 2.1 Slices (step 1) [P]

`theory/open_screen_enum.py` enumerates, at every open order 4v, the staged slices of three families as
`theory/ab_lib.Slice` objects:

| family | group G | target | applies |
|---|---|---|---|
| b0 | abelian, order v, G = G′ × ℤ_q | Σ_k PAF_k = 4v·δ | v odd composite |
| b1 | abelian, order v − 1, G = G′ × ℤ_q | 4m at 0, −4 elsewhere | always |
| k4 | order v − 3 = 4h, G = Q × H, Q = ℤ₂² or ℤ₄, K = H | 4m at 0, −12 on K \ 0, +4 off K | v ≡ 3 (mod 4) |

- The staging prime q exactly divides |G| (so ℤ_q is a direct factor); G′ runs over every abelian group of order
  |G|/q (for k4, H′ over every abelian group of order h/q).
- The symmetry is S = U × M, with M the index-s subgroup of F_q^* acting on ℤ_q and U trivial, one cyclic subgroup of
  units on one axis of G′, or −1 on several axes (the menu of `theory/ab_sweep.py`).
- Orbits per block: r = n_B (1 + s), with n_B the number of U-orbits on G′ [P: M acts freely on ℤ_q \ 0 and U fixes
  the ℤ_q coordinate]. Kept: 24 ≤ r ≤ 160. The prompt's "about 40" was widened down to 24 because the 933 success
  has r = 33.
- sGeff = n_B·s is the number of stage-2 class equations up to a factor of about 2; it orders the work.
- **Not enumerated:** diagonal multipliers acting on several axes with one generator (the 647 type), and the
  families k6 and k8. Each is recorded as not covered (§7).

### 2.2 Exact gates (step 2) [P, M]

Cheapest first, per slice (`theory/open_screen_enum.py`, logs `logs/open_screen/screens_s*.jsonl`):
- **A**: stage-1 lag-0 energy and row sums, and for k4 the Q-level, by the DP and 2+2 join of
  `ab_sweep.screen_A`. For b1 the row condition is Σ r_k² = 4 (passed to screen_A as v = 1).
- **B**: α-congruence (`ab_sweep.screen_B`, re-implemented with a vectorised join in `screen_B_fast`; B is left
  undecided when the α-keys exceed 4000).
- Controls [M]: 933 (ℤ₃ × ℤ₃₁₁, s = 10) passes A and B; 523 and 235 (ℤ₂ × ℤ₉ × ℤ_q, s = 4) pass A; 739 ×127
  (ℤ₂ × ℤ₉ × ℤ₄₁, s = 4) fails A, as proved in `report/paley_523_answer.md` §7.

### 2.3 Costs (step 3) [E]

`theory/open_screen_deficit.py` generalises the 523 / 933 / 863 deficit model to any slice and any split:
- **Stage 1.** `ab_gauss.stage1`: the exact count of U-invariant y with the stage-1 lag-0 energy and the row sums,
  times the lattice Gaussian of the other stage-1 lag classes, read at the true target. Where its (E, r) tables
  exceed memory, an energy-only fallback is used (flagged `stage1_method: energy-only`). It counts y exactly by
  energy alone and leaves the row condition to the Gaussian; on the one slice where both ran (599 b1,
  ℤ₂₃ × ℤ₂ × ℤ₁₃, s = 4) it is 6.6 bits optimistic (25.2 against 18.6).
- **Stage 2, given a y.** A synthetic y is sampled uniformly among the U-invariant y with the stage-1 energy and
  rows; real y are used in step 5. The estimator then gives:
  - H₂ = log₂ #fibre (exact);
  - C₂ = the Gaussian bits of every class equation, at the sample mean for a synthetic y, at the target for a real
    one;
  - E per y = H₂ − C₂.
- **Splits.** For every S-invariant subgroup L of G′ (the 10 smallest, and G′ itself):
  - the S-level is the compression along L; C_S is its Gaussian bits;
  - H_S = H(S | y) is exact. It is computed by convolving each cell's pattern distribution, mapped linearly to its
    contribution at one coset per S-orbit of G/L, over each connected component of cells. When that map is
    certified injective (distinct signed subset sums per column), H_S is the sum of the cell entropies [P]. Past a
    size cap the lattice-Gaussian entropy is used, clipped at the cell-entropy sum and flagged `HS_gauss`.
  - h_D = H₂ − H_S (the mean D-fibre); C(D|S) = C₂ − C_S; **D_def = C(D|S) − h_D**.
- **Cost per core:** log₂ cost = D_def + log₂(2^S-step + 2^(h_D/2)). S-step = max(0, C_S − H_S/2) is a plain 2+2 join
  of the S-level with full lists; h_D/2 is the 2+2 D-check.
- **Expected cores in the slice:** log₂ E_tot = stage-1 log₂E + E per y.
- **Regime** (`report/h863_lift_answer.md` §5): D_def ≤ 26, h_D ≤ 40, E per y ≥ 10, C₂ ≤ 150.

**Calibration against the known successes [M sizes, E bits; `--controls`]:**

| control | split | quantity | this estimator | recorded |
|---|---|---|---|---|
| 933, ℤ₃ × ℤ₃₁₁, s = 10, real y | L = G′ (dual first) | C₂ / C_S / C(D\|S) / D_def / cost | 67.4 / 25.6 / 41.7 / 7.3 / 24.5 | 67.5 / 25.6 / ≈ 42 / 7.4 / 25.7 measured (`report/alpoge_design_answer.md` §4) |
| 523, ℤ₂ × ℤ₉ × ℤ₂₉, s = 4, production y | L = ℤ₂ | C₂ / h_D / C(D\|S) / D_def / cost | 146.25 / 36.43 / 62.53 / 26.10 / 44.3 | 146.2 / 36.44 / 62.5 / 26.14 / 44.2 (`report/paley_523_answer.md` §3.6) |
| 863 k4 S7, ℤ₂² × ℤ₅ × ℤ₄₃, synthetic y | L = ℤ₂ / ℤ₅ | C₂ / E per y / D_def | 275.5 / 64.5 / 51.0 / 34.6 | 275.5 / 63.5 / 51.5 / 35.4 (`report/h863_lift_answer.md` §0, §4) |
| 235, the 523 shape at q = 13, synthetic y | L = ℤ₂ | D_def | 24.3 | 23.3 (`report/paley_523_answer.md` §4) |

The cost formula reproduces both production costs (933 and 523) to about 1 bit.

## 3. Ranked table (step 4)

<!-- TABLE -->
`logs/open_screen/table.jsonl` has one line per order and family. Each line gives the best slice; the gate counts of
every slice; stage-1 log₂E; H₂, C₂ and E per y; the best split with C_S, H_S, D_def, h_D, S-step and cost; log₂ E_tot;
the regime tests; the exact stage-1 status; and the next gate with its log₂ cost (the stage-1 join, ≈ log₂ N_y / 2).
It also holds the b0-direct, Cooper–Wallis and ℤ_q lines of §5. The table is rebuilt by
`python3 theory/open_screen_table.py --md FILE`.

**Tiers.**
- 0: in the regime (D_def ≤ 26, h_D ≤ 40, E per y ≥ 10, C₂ ≤ 150) and log₂ E_tot ≥ 0;
- 1: in the regime, but log₂ E_tot < 0 (stage 1 is expected empty);
- 2: log₂ E_tot ≥ 0, but outside the regime;
- 3: the rest.

The best staged slice of each of the 85 orders is in tier 0 at 10 orders, tier 1 at 5, tier 2 at 39 and tier 3 at 10.
The other 21 orders have no costed alive slice. Most of those have no staged slice with 24 ≤ r ≤ 160, or every slice
died at an exact gate; at the rest, every alive slice has sGeff > 101 (so C₂ > 150 [E]).

**The tier-0 slices, the best per order** (cost = log₂ operations per core; † = stage-1 estimate by the energy-only
fallback, which is optimistic):

| 4v | v | fam | slice | stage-1 log₂E | E/y | C₂ | D_def | h_D | cost | log₂ E_tot | exact stage 1 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 5804 | 1451 | b1 | ℤ₂₅ × ℤ₂ × ℤ₂₉, s = 4, U order 10 | 0.9† | 11.8 | 83.4 | 15.6 | 18.5 | 27.0 | 12.7 | **empty, method B only** (§4) |
| 4156 | 1039 | k4 | ℤ₄ × ℤ₇ × ℤ₃₇, s = 4, U = −1 | −22.0† | 31.3 | 119.7 | 17.5 | 28.9 | 32.2 | 9.3 | not run |
| 5708 | 1427 | b1 | ℤ₂₃ × ℤ₂ × ℤ₃₁, s = 10, U order 22 | 7.9† | 13.3 | 99.3 | 20.4 | 21.2 | 32.9 | 21.2 | not run |
| 3452 | 863 | k4 | ℤ₄ × ℤ₅ × ℤ₄₃, s = 6, U = −1 | −3.9 | 47.9 | 138.6 | 17.8 | 34.5 | 35.1 | 43.9 | <!-- K863T -->see §4<!-- /K863T --> |
| 5036 | 1259 | b1 | ℤ₁₇ × ℤ₂ × ℤ₃₇, s = 12, U order 16 | 15.6† | 23.0 | 120.4 | 24.2 | 27.4 | 38.1 | 38.5 | not run (its s = 4 slice is empty, §4) |
| 7572 | 1893 | b1 | ℤ₁₁ × ℤ₂² × ℤ₄₃, s = 6, U order 10 | 5.5† | 13.7 | 117.0 | 24.2 | 25.9 | 39.1 | 19.1 | not run |
| 5692 | 1423 | b1 | ℤ₃² × ℤ₂ × ℤ₇₉, s = 6, U = −1 | −19.6† | 38.2 | 132.6 | 21.6 | 35.5 | 39.4 | 18.5 | not run |
| 9260 | 2315 | b1 | ℤ₁₃ × ℤ₂ × ℤ₈₉, s = 8, U order 6 | −6.9† | 17.2 | 123.9 | 25.1 | 28.7 | 40.2 | 10.3 | not run |
| 3788 | 947 | b1 | ℤ₁₁ × ℤ₂ × ℤ₄₃, s = 14, U order 10 | 18.2† | 42.5 | 138.6 | 22.0 | 38.3 | 41.1 | 60.7 | not run (its s = 6 slice is empty, §4) |
| 4748 | 1187 | k4 | ℤ₄ × ℤ₈ × ℤ₃₇, s = 4, U = −1 | −22.2† | 32.7 | 148.2 | 21.5 | 33.0 | 41.4 | 10.4 | not run |

- All the splits are a ℤ₂ (or ℤ₂ in Q) S/D split, as at 523. The b0 dual-first (L = G′) split wins only at 1719.
- All tier-0 slices are b1 or k4 with a nontrivial U. The pure GP slices with U = 1, the 523 shape, die at the ℤ₂
  screen, the w-screen or screen A, or have stage-1 estimates far below zero.

**Full table (best slice per order and family; the costed families):**

| 4v | v | fam | slices: dead by gate / alive | best slice (tier) | s1 log2E | E/y | C2 | split | D_def | h_D | cost | log2 E_tot | exact stage 1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2396 | 599 | b1 | 16: alive* 5, alive 10, w 1 | G'=Z13xZ2, q=23, s=2, f=11 U=<12> (order 2) on axis 0, r=42 (1) | -37.51 | 12.88 | 40.12 | L2<(0, 1)> | 2.76 | 12.5 | 9.41 | -24.63 | - |
| 2396 | 599 | k4 | 3: A 3 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 2676 | 669 | b0 | 3: A 1, alive 2 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 2676 | 669 | b1 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 2876 | 719 | b1 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 2876 | 719 | k4 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 2956 | 739 | b1 | 36: A 23, B 1, alive 12 | G'=Z9xZ2, q=41, s=8, f=5 U=<8> (order 2) on axis 0, r=90 (2) | 11.31 | 65.66 | 174.0 | L2<(0, 1)> | 25.5 | 49.94 | 50.47 | 76.97 | - |
| 2956 | 739 | k4 | 19: A 19 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 3004 | 751 | b1 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 3004 | 751 | k4 | 37: A 10, alive* 15, alive 12 | G'=Z4xZ17, q=11, s=2, f=5 U=<2> (order 8) on axis 1, r=36 (2) | 12.52† | -8.03 | 45.03 | L34<(0, 1),(0, 2)> | 12.97 | 13.91 | 20.76 | 4.49 | - |
| 3452 | 863 | b1 | 1: A 1 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 3452 | 863 | k4 | 44: A 26, alive* 8, alive 10 | G'=Z4xZ5, q=43, s=6, f=7 U=-1 on axes [0, 1], r=77 (0) | -3.94 | 47.85 | 138.58 | L2<(2, 0)> | 17.81 | 34.51 | 35.08 | 43.91 | - |
| 3628 | 907 | b1 | 12: A 7, B 2, Z2 1, alive 2 | G'=Z3xZ2, q=151, s=25, f=6 U=<2> (order 2) on axis 0, r=104 (3) | 8.15 | -153.19 | 508.82 | G' | 217.56 | 102.65 | 279.69 | -145.04 | - |
| 3628 | 907 | k4 | 17: A 17 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 3788 | 947 | b1 | 26: A 3, B 3, Z2 1, alive* 4, alive 14, S1 1 | G'=Z11xZ2, q=43, s=14, f=3 U=<2> (order 10) on axis 0, r=60 (0) | 18.22† | 42.48 | 138.59 | L2<(0, 1)> | 22.0 | 38.28 | 41.14 | 60.7 | - |
| 3788 | 947 | k4 | 9: A 9 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 3884 | 971 | b1 | 26: A 8, B 3, alive* 3, alive 12 | G'=Z5xZ2, q=97, s=8, f=12, r=90 (2) | 22.34 | 13.96 | 229.28 | L5<(1, 0),(2, 0)> | 7.31 | 126.66 | 70.64 | 36.3 | - |
| 3884 | 971 | k4 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 4076 | 1019 | b1 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 4076 | 1019 | k4 | 21: A 20, alive 1 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 4124 | 1031 | b1 | 12: A 6, alive* 1, alive 5 | G'=Z103xZ2, q=5, s=2, f=2 U=<8> (order 17) on axis 0, r=42 (2) | 61.72† | -39.62 | 77.62 | G' | 30.6 | 22.84 | 47.24 | 22.1 | - |
| 4124 | 1031 | k4 | 9: A 9 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 4156 | 1039 | b1 | 1: A 1 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 4156 | 1039 | k4 | 77: A 36, alive* 22, alive 19 | G'=Z4xZ7, q=37, s=4, f=9 U=-1 on axes [0, 1], r=75 (0) | -21.97† | 31.3 | 119.7 | L2<(2, 0)> | 17.51 | 28.9 | 32.24 | 9.33 | - |
| 4348 | 1087 | b1 | 20: A 10, B 6, Z2 2, alive 2 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 4348 | 1087 | k4 | 23: A 23 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 4412 | 1103 | b1 | 31: A 9, B 1, alive* 4, alive 17 | G'=Z19xZ2, q=29, s=4, f=7 U=<4> (order 9) on axis 0, r=30 (2) | 6.77 | 3.81 | 54.13 | L2<(0, 1)> | 11.08 | 11.06 | 19.8 | 10.58 | - |
| 4412 | 1103 | k4 | 22: A 2, alive* 16, alive 4 | G'=Z4xZ25, q=11, s=2, f=5 U=<4> (order 10) on axis 1, r=60 (2) | 36.17† | -17.44 | 75.44 | L10<(0, 5),(0, 10)> | 22.5 | 15.77 | 38.56 | 18.73 | - |
| 4460 | 1115 | b0 | 5: A 2, alive* 1, alive 2 | G'=Z223, q=5, s=2, f=2 U=<40> (order 6) on axis 0, r=114 (3) | 66.51† | -239.25 | 332.25 | G' | 238.77 | 79.74 | 278.64 | -172.74 | - |
| 4460 | 1115 | b1 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 4460 | 1115 | k4 | 9: A 8, alive 1 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 4676 | 1169 | b0 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 4676 | 1169 | b1 | 115: A 50, B 26, alive* 2, alive 37 | G'=Z8xZ2, q=73, s=6, f=12, r=112 (2) | 22.86† | -9.98 | 283.26 | L4<(2, 0),(4, 0)> | 26.78 | 122.1 | 88.1 | 12.88 | - |
| 4748 | 1187 | b1 | 3: A 2, alive 1 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 4748 | 1187 | k4 | 63: A 38, alive* 25 | G'=Z4xZ8, q=37, s=4, f=9 U=-1 on axes [0, 1], r=90 (0) | -22.23† | 32.65 | 148.2 | L2<(2, 0)> | 21.51 | 32.95 | 41.43 | 10.42 | - |
| 4792 | 1198 | b1 | 45: A 14, alive* 27, alive 4 | G'=Z9xZ7, q=19, s=2, f=9 U=<3> (order 6) on axis 1, r=54 (2) | 46.15† | -8.96 | 66.96 | L3<(3, 0),(6, 0)> | 16.96 | 21.02 | 28.46 | 37.19 | - |
| 4892 | 1223 | b1 | 11: A 4, alive* 1, alive 6 | G'=Z47xZ2, q=13, s=4, f=3 U=<2> (order 23) on axis 0, r=30 (2) | 14.79 | -12.46 | 58.31 | G' | 17.3 | 16.84 | 27.47 | 2.33 | - |
| 4892 | 1223 | k4 | 68: A 44, alive* 10, alive 14 | G'=Z2xZ2xZ5, q=61, s=12, f=5 U=<2> (order 4) on axis 2, r=104 (2) | 14.58† | 67.98 | 241.87 | L2<(1, 1, 0)> | 40.18 | 66.76 | 73.56 | 82.56 | - |
| 5028 | 1257 | b0 | 8: A 6, alive 2 | G'=Z3, q=419, s=38, f=11 U=<2> (order 2) on axis 0, r=78 (2) | 11.17 | 82.08 | 197.65 | G' | 71.24 | 0.0 | 72.24 | 93.25 | - |
| 5028 | 1257 | b1 | 38: A 22, B 10, alive 6 | G'=Z2xZ2xZ2, q=157, s=12, f=13, r=104 (2) | -8.29 | 84.19 | 229.73 | L2<(1, 0, 1)> | 32.63 | 68.98 | 67.12 | 75.9 | - |
| 5036 | 1259 | b1 | 39: A 9, B 2, alive* 9, alive 18, S1 1 | G'=Z17xZ2, q=37, s=12, f=3 U=<3> (order 16) on axis 0, r=52 (0) | 15.59† | 22.95 | 120.35 | L2<(0, 1)> | 24.21 | 27.38 | 38.08 | 38.54 | - |
| 5036 | 1259 | k4 | 17: A 17 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 5132 | 1283 | b1 | 5: A 4, alive 1 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 5132 | 1283 | k4 | 4: alive* 4 | G'=Z4xZ64, q=5, s=2, f=2 U=<5> (order 16) on axis 1, r=144 (3) | -559.72† | -148.99 | 248.99 | L4<(1, 48),(2, 32)> | 128.32 | 33.46 | 182.26 | -708.71 | - |
| 5164 | 1291 | b1 | 23: A 10, alive* 9, alive 4 | G'=Z5xZ3xZ2, q=43, s=6, f=7 U=<2> (order 4) on axis 0, r=84 (2) | 34.52† | 22.24 | 176.83 | L2<(0, 0, 1)> | 37.73 | 39.76 | 58.51 | 56.76 | - |
| 5164 | 1291 | k4 | 16: A 2, alive* 12, alive 2 | G'=Z4xZ7xZ2, q=23, s=2, f=11 U=<3> (order 6) on axis 1, r=48 (2) | 6.67† | -4.11 | 58.11 | L2<(2, 0, 0)> | 12.12 | 12.0 | 25.12 | 2.56 | - |
| 5276 | 1319 | b1 | 2: A 2 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 5276 | 1319 | k4 | 18: A 18 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 5308 | 1327 | b1 | 24: B 1, alive* 15, alive 8 | G'=Z13xZ3xZ2, q=17, s=2, f=8 U=<3> (order 3) on axis 0, r=90 (2) | 51.2† | -24.34 | 132.34 | L13<(1, 0, 0),(2, 0, 0)> | 31.53 | 44.14 | 56.48 | 26.86 | - |
| 5308 | 1327 | k4 | 22: A 22 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 5468 | 1367 | b1 | 4: A 3, alive 1 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 5468 | 1367 | k4 | 55: A 32, alive* 13, alive 10 | G'=Z4xZ31, q=11, s=2, f=5 U=<6> (order 6) on axis 1, r=72 (2) | 22.2† | -20.52 | 99.52 | L31<(0, 1),(0, 2)> | 27.87 | 39.05 | 47.41 | 1.68 | - |
| 5692 | 1423 | b1 | 27: A 15, B 2, alive* 2, alive 8 | G'=Z3xZ3xZ2, q=79, s=6, f=13 U=-1 on axes [0, 1], r=70 (0) | -19.63† | 38.16 | 132.64 | L2<(0, 0, 1)> | 21.64 | 35.52 | 39.4 | 18.53 | - |
| 5692 | 1423 | k4 | 41: A 32, alive* 4, alive 5 | G'=Z4xZ71, q=5, s=2, f=2 U=<20> (order 7) on axis 1, r=132 (2) | 224.57† | -117.18 | 219.18 | L71<(0, 1),(0, 2)> | 103.53 | 61.63 | 137.54 | 107.39 | - |
| 5708 | 1427 | b1 | 22: A 6, B 1, alive* 4, alive 10, w 1 | G'=Z23xZ2, q=31, s=10, f=3 U=<5> (order 22) on axis 0, r=44 (0) | 7.89† | 13.33 | 99.28 | L2<(0, 1)> | 20.4 | 21.17 | 32.85 | 21.22 | - |
| 5708 | 1427 | k4 | 32: A 18, B 8, alive* 1, alive 5 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 5748 | 1437 | b0 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 5748 | 1437 | b1 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 5752 | 1438 | b1 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 5756 | 1439 | b1 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 5756 | 1439 | k4 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 5804 | 1451 | b1 | 29: A 11, B 1, alive* 8, alive 9 | G'=Z25xZ2, q=29, s=4, f=7 U=<4> (order 10) on axis 0, r=50 (0) | 0.89† | 11.78 | 83.43 | L2<(0, 1)> | 15.6 | 18.47 | 26.97 | 12.67 | stage 1: one method only |
| 5804 | 1451 | k4 | 30: A 30 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 5884 | 1471 | b1 | 3: alive* 2, alive 1 | G'=Z49xZ3xZ2, q=5, s=2, f=2 U=<3> (order 42) on axis 0, r=54 (3) | -1280.92† | -67.31 | 101.31 | L6<(0, 0, 1),(0, 1, 0)> | 51.54 | 14.16 | 77.23 | -1348.23 | - |
| 5884 | 1471 | k4 | 2: A 2 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 5892 | 1473 | b0 | 9: A 6, alive 3 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 5892 | 1473 | b1 | 56: B 2, alive* 52, alive 2 | G'=Z32xZ2, q=23, s=2, f=11 U=<3> (order 8) on axis 0, r=54 (2) | 11.19† | -1.39 | 66.39 | L4<(0, 1),(16, 0)> | 11.02 | 20.75 | 23.81 | 9.8 | - |
| 5948 | 1487 | b1 | 2: A 2 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 5948 | 1487 | k4 | 34: A 24, alive* 6, alive 4 | G'=Z4xZ7, q=53, s=4, f=13 U=-1 on axes [0, 1], r=75 (1) | -37.27† | 30.37 | 121.8 | L2<(2, 0)> | 17.67 | 29.51 | 32.87 | -6.9 | - |
| 5996 | 1499 | b1 | 3: A 2, alive 1 | G'=Z107xZ2, q=7, s=3, f=2 U=<3> (order 53) on axis 0, r=24 (3) | 17.44 | -48.28 | 83.15 | G' | 29.46 | 10.67 | 60.38 | -30.84 | - |
| 5996 | 1499 | k4 | 23: alive* 16, alive 7 | G'=Z4xZ17xZ2, q=11, s=2, f=5 U=<4> (order 4) on axis 1, r=120 (2) | 41.26† | -39.54 | 172.54 | L17<(0, 1, 0),(0, 2, 0)> | 54.05 | 53.74 | 81.29 | 1.72 | - |
| 6172 | 1543 | b1 | 6: A 5, B 1 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 6172 | 1543 | k4 | 6: alive* 6 | G'=Z2xZ2xZ7xZ5, q=11, s=2, f=5 U=<3> (order 6) on axis 2, r=120 (2) | 103.79† | -65.97 | 175.97 | L4<(0, 1, 0, 0),(1, 0, 0, 0)> | 74.05 | 45.77 | 98.62 | 37.82 | - |
| 6236 | 1559 | b1 | 43: A 15, B 4, alive* 6, alive 18 | G'=Z41xZ2, q=19, s=2, f=9 U=<10> (order 5) on axis 0, r=54 (2) | 20.92† | -7.24 | 71.24 | L2<(0, 1)> | 13.45 | 15.0 | 31.74 | 13.68 | - |
| 6236 | 1559 | k4 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 6268 | 1567 | b1 | 35: B 3, alive* 26, alive 6 | G'=Z27xZ2, q=29, s=2, f=14 U=<10> (order 3) on axis 0, r=90 (2) | 35.65† | -11.18 | 123.18 | L6<(0, 1),(9, 0)> | 20.06 | 49.18 | 44.96 | 24.47 | - |
| 6268 | 1567 | k4 | 29: A 8, B 2, alive* 9, alive 10 | G'=Z4xZ17, q=23, s=2, f=11 U=<3> (order 16) on axis 1, r=24 (2) | 5.77 | 0.1 | 19.9 | G' | 0.91 | 8.38 | 6.44 | 5.87 | - |
| 6284 | 1571 | b1 | 21: A 5, B 4, Z2 1, alive* 2, alive 9 | G'=Z157xZ2, q=5, s=2, f=2 U=<9> (order 39) on axis 0, r=30 (2) | 35.56 | -30.62 | 60.62 | G' | 23.66 | 16.17 | 37.56 | 4.94 | - |
| 6284 | 1571 | k4 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 6316 | 1579 | b1 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 6316 | 1579 | k4 | 11: A 11 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 6332 | 1583 | b1 | 31: A 12, B 4, alive* 5, alive 10 | G'=Z7xZ2, q=113, s=16, f=7 U=<3> (order 6) on axis 0, r=68 (2) | 2.09 | 47.19 | 169.02 | L2<(0, 1)> | 32.53 | 42.14 | 53.6 | 49.28 | - |
| 6332 | 1583 | k4 | 35: A 34, alive 1 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 6476 | 1619 | b1 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 6476 | 1619 | k4 | 32: A 27, B 5 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 6652 | 1663 | b1 | 8: A 6, B 1, alive 1 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 6652 | 1663 | k4 | 12: A 12 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 6668 | 1667 | b1 | 21: alive* 15, alive 6 | G'=Z49xZ2, q=17, s=2, f=8 U=<8> (order 7) on axis 0, r=78 (2) | 73.61† | -36.06 | 116.06 | L49<(1, 0),(2, 0)> | 32.55 | 45.78 | 55.71 | 37.55 | - |
| 6668 | 1667 | k4 | 10: alive* 10 | G'=Z4xZ32, q=13, s=2, f=6 U=<5> (order 8) on axis 1, r=120 (2) | 102.12† | -77.43 | 201.43 | L4<(1, 8),(2, 16)> | 76.92 | 52.93 | 112.97 | 24.69 | - |
| 6852 | 1713 | b0 | 10: A 9, alive 1 | G'=Z3, q=571, s=30, f=19, r=93 (2) | 16.92 | 82.94 | 243.1 | G' | 30.81 | 120.05 | 90.83 | 99.86 | - |
| 6852 | 1713 | b1 | 16: A 16 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 6876 | 1719 | b0 | 21: A 9, B 2, alive 7, S1 3 | G'=Z3xZ3, q=191, s=10, f=19, r=99 (2) | 17.17† | 53.14 | 232.73 | L3<(1, 1),(2, 2)> | 40.76 | 102.75 | 92.13 | 70.31 | - |
| 6876 | 1719 | b1 | 8: A 8 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 6876 | 1719 | k4 | 18: alive* 18 | G'=Z4xZ13xZ3, q=11, s=2, f=5 U=<4> (order 6) on axis 1, r=108 (2) | 77.39† | -59.65 | 168.65 | L26<(0, 1, 0),(0, 2, 0)> | 69.11 | 43.8 | 92.76 | 17.74 | - |
| 6988 | 1747 | b1 | 45: A 23, B 5, alive* 2, alive 15 | G'=Z3xZ3xZ2, q=97, s=6, f=16 U=<2> (order 2) on axis 0, r=84 (3) | 2.15† | -37.73 | 243.52 | L6<(0, 0, 1),(0, 1, 0)> | 43.62 | 117.87 | 102.56 | -35.58 | - |
| 6988 | 1747 | k4 | 50: A 45, alive 5 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 7132 | 1783 | b1 | 12: alive* 9, alive 3 | G'=Z81xZ2, q=11, s=2, f=5 U=<10> (order 9) on axis 0, r=126 (2) | 119.81† | -65.48 | 184.48 | L18<(0, 1),(9, 0)> | 64.79 | 55.62 | 97.23 | 54.33 | - |
| 7132 | 1783 | k4 | 37: A 26, B 3, alive* 4, alive 4 | G'=Z4xZ5, q=89, s=8, f=11 U=-1 on axes [0, 1], r=99 (2) | -17.27† | 54.31 | 213.12 | L2<(2, 0)> | 31.65 | 51.35 | 57.44 | 37.04 | - |
| 7148 | 1787 | b1 | 12: A 5, alive* 1, alive 6 | G'=Z47xZ2, q=19, s=6, f=3 U=<2> (order 23) on axis 0, r=42 (2) | 27.07† | -19.13 | 97.45 | G' | 29.09 | 30.4 | 44.81 | 7.94 | - |
| 7148 | 1787 | k4 | 8: A 8 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 7180 | 1795 | b0 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 7180 | 1795 | b1 | 16: alive* 8, alive 8 | G'=Z13xZ3xZ2, q=23, s=2, f=11 U=<2> (order 12) on axis 0, r=36 (2) | 21.5† | -9.02 | 43.02 | L39<(0, 1, 0),(0, 2, 0)> | 11.92 | 11.2 | 20.6 | 12.48 | - |
| 7180 | 1795 | k4 | 4: alive* 4 | G'=Z4xZ64, q=7, s=2, f=3 U=<3> (order 16) on axis 1, r=132 (3) | 59.42† | -103.93 | 207.93 | L4<(1, 0),(2, 0)> | 105.63 | 38.14 | 136.86 | -44.51 | - |
| 7324 | 1831 | b1 | 31: A 9, B 1, Z2 1, alive* 14, alive 6 | G'=Z5xZ3xZ2, q=61, s=4, f=15 U=-1 on axes [0, 1], r=80 (1) | -45.59† | 27.79 | 134.72 | L2<(0, 0, 1)> | 24.18 | 33.59 | 41.04 | -17.8 | - |
| 7324 | 1831 | k4 | 17: A 17 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 7388 | 1847 | b1 | 44: A 16, B 5, alive* 8, alive 15 | G'=Z71xZ2, q=13, s=2, f=6 U=<20> (order 7) on axis 0, r=66 (2) | 36.56† | -23.02 | 97.02 | L71<(1, 0),(2, 0)> | 19.84 | 52.81 | 46.25 | 13.54 | - |
| 7388 | 1847 | k4 | 12: A 12 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 7484 | 1871 | b1 | 10: alive* 9, alive 1 | G'=Z11xZ5xZ2, q=17, s=2, f=8 U=<2> (order 10) on axis 0, r=60 (2) | 44.58† | -39.78 | 102.78 | L5<(0, 1, 0),(0, 2, 0)> | 27.0 | 31.18 | 55.69 | 4.8 | - |
| 7484 | 1871 | k4 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 7516 | 1879 | b1 | 14: A 10, B 3, alive 1 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 7516 | 1879 | k4 | 50: A 40, alive* 8, alive 2 | G'=Z4xZ7, q=67, s=6, f=11 U=-1 on axes [0, 1], r=105 (2) | -29.2† | 44.09 | 211.49 | L2<(2, 0)> | 35.79 | 49.44 | 60.94 | 14.89 | - |
| 7572 | 1893 | b0 | 21: A 18, alive 3 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 7572 | 1893 | b1 | 44: A 10, alive* 10, alive 22, w 2 | G'=Z11xZ2xZ2, q=43, s=6, f=7 U=<2> (order 10) on axis 0, r=56 (0) | 5.48† | 13.66 | 117.02 | L2<(0, 0, 1)> | 24.17 | 25.9 | 39.14 | 19.14 | - |
| 7628 | 1907 | b1 | 6: A 6 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 7628 | 1907 | k4 | 24: alive* 24 | G'=Z4xZ7xZ4, q=17, s=2, f=8 U=<2> (order 3) on axis 1, r=144 (2) | 79.61† | -61.73 | 231.73 | L4<(0, 0, 2),(2, 0, 0)> | 63.78 | 77.91 | 107.82 | 17.88 | - |
| 7716 | 1929 | b0 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 7716 | 1929 | b1 | 74: A 50, B 18, Z2 3, alive 3 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 7876 | 1969 | b0 | 1: A 1 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 7876 | 1969 | b1 | 51: A 33, alive* 18 | G'=Z16xZ3, q=41, s=4, f=10 U=<5> (order 4) on axis 0, r=120 (2) | 52.65† | -27.53 | 265.74 | L6<(0, 1),(0, 2)> | 49.81 | 106.73 | 103.18 | 25.12 | - |
| 7916 | 1979 | b1 | 20: A 6, alive* 3, alive 11 | G'=Z23xZ2, q=43, s=6, f=7 U=<2> (order 11) on axis 0, r=42 (2) | 15.25† | 0.35 | 95.32 | L2<(0, 1)> | 22.06 | 18.67 | 38.16 | 15.6 | - |
| 7916 | 1979 | k4 | 40: A 8, alive* 23, alive 9 | G'=Z4xZ13xZ2, q=19, s=2, f=9 U=<2> (order 12) on axis 1, r=48 (2) | 19.12† | -15.21 | 56.21 | L4<(1, 0, 1),(2, 0, 0)> | 18.72 | 17.94 | 28.29 | 3.91 | - |
| 7924 | 1981 | b0 | 7: A 3, alive* 1, alive 3 | G'=Z283, q=7, s=3, f=2 U=<4> (order 47) on axis 0, r=28 (2) | 34.98 | -30.5 | 68.54 | G' | 11.75 | 17.01 | 41.01 | 4.48 | - |
| 7924 | 1981 | b1 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 7948 | 1987 | b1 | 15: A 13, B 2 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 7948 | 1987 | k4 | 55: A 20, alive* 35 | G'=Z4xZ16, q=31, s=2, f=15 U=-1 on axes [0, 1], r=102 (1) | -99.37† | 15.39 | 115.61 | L2<(0, 8)> | 17.56 | 30.5 | 35.17 | -83.98 | - |
| 8108 | 2027 | b1 | 5: A 3, alive 2 | G'=Z2, q=1013, s=44, f=23, r=90 (2) | 7.58 | 85.04 | 242.32 | G' | 39.3 | 69.85 | 74.23 | 92.62 | - |
| 8108 | 2027 | k4 | 15: alive* 13, alive 2 | G'=Z4xZ11xZ2, q=23, s=2, f=11 U=<3> (order 5) on axis 1, r=72 (2) | 29.71† | -19.97 | 100.97 | L4<(1, 0, 1),(2, 0, 0)> | 30.86 | 37.36 | 49.55 | 9.74 | - |
| 8156 | 2039 | b1 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 8156 | 2039 | k4 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 8252 | 2063 | b1 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 8252 | 2063 | k4 | 34: A 34 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 8332 | 2083 | b1 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 8332 | 2083 | k4 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 8348 | 2087 | b1 | 9: A 6, alive 3 | G'=Z149xZ2, q=7, s=3, f=2 U=<5> (order 37) on axis 0, r=40 (3) | 42.06† | -55.29 | 113.93 | G' | 38.85 | 30.99 | 69.11 | -13.23 | - |
| 8348 | 2087 | k4 | 19: A 19 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 8380 | 2095 | b0 | 15: A 10, alive* 2, alive 3 | G'=Z419, q=5, s=2, f=2 U=<7> (order 19) on axis 0, r=69 (2) | 116.1† | -63.38 | 121.38 | G' | 54.42 | 44.22 | 76.55 | 52.72 | - |
| 8380 | 2095 | b1 | 7: A 4, B 2, Z2 1 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 8380 | 2095 | k4 | 11: A 11 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 8396 | 2099 | b1 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 8396 | 2099 | k4 | 24: A 18, B 5, alive 1 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 8684 | 2171 | b0 | 1: A 1 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 8684 | 2171 | b1 | 21: A 5, alive* 15, alive 1 | G'=Z7xZ5xZ2, q=31, s=2, f=15 U=<3> (order 6) on axis 0, r=60 (2) | 22.41† | -7.62 | 78.62 | L2<(0, 0, 1)> | 16.43 | 16.5 | 34.87 | 14.79 | - |
| 8684 | 2171 | k4 | 23: A 23 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 8812 | 2203 | b1 | 3: A 3 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 8812 | 2203 | k4 | 6: alive* 6 | G'=Z2xZ2xZ25xZ2, q=11, s=2, f=5 U=<2> (order 20) on axis 2, r=72 (3) | -49.33† | -33.74 | 93.74 | L4<(0, 0, 0, 1),(1, 0, 0, 0)> | 37.55 | 24.2 | 51.96 | -83.07 | - |
| 8828 | 2207 | b1 | 4: A 4 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 8828 | 2207 | k4 | 53: A 28, alive* 15, alive 10 | G'=Z4xZ29, q=19, s=2, f=9 U=<7> (order 7) on axis 1, r=60 (2) | 23.38† | -16.75 | 88.75 | L29<(0, 1),(0, 2)> | 24.7 | 29.17 | 39.83 | 6.63 | - |
| 8860 | 2215 | b0 | 15: A 11, alive* 3, alive 1 | G'=Z443, q=5, s=2, f=2 U=<13> (order 17) on axis 0, r=81 (2) | 136.93† | -72.24 | 142.24 | G' | 63.4 | 55.65 | 91.23 | 64.69 | - |
| 8860 | 2215 | b1 | 39: A 21, alive* 12, alive 6 | G'=Z27xZ2, q=41, s=8, f=5 U=<2> (order 18) on axis 0, r=72 (2) | 8.79† | 20.72 | 167.27 | L2<(0, 1)> | 34.55 | 38.46 | 54.92 | 29.51 | - |
| 8860 | 2215 | k4 | 46: A 38, alive* 8 | G'=Z4xZ7, q=79, s=6, f=13 U=-1 on axes [0, 1], r=105 (2) | -37.24† | 41.28 | 214.72 | L2<(2, 0)> | 36.75 | 49.89 | 62.74 | 4.04 | - |
| 9004 | 2251 | b1 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 9004 | 2251 | k4 | 21: A 21 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 9148 | 2287 | b1 | 43: A 37, alive 6 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 9148 | 2287 | k4 | 19: A 19 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 9260 | 2315 | b0 | 28: A 17, alive* 4, alive 7 | G'=Z463, q=5, s=2, f=2 U=<36> (order 33) on axis 0, r=45 (2) | 73.07† | -48.67 | 82.67 | G' | 36.42 | 21.86 | 54.75 | 24.4 | - |
| 9260 | 2315 | b1 | 36: A 10, alive* 7, alive 19 | G'=Z13xZ2, q=89, s=8, f=11 U=<4> (order 6) on axis 0, r=54 (0) | -6.88† | 17.19 | 123.85 | L2<(0, 1)> | 25.11 | 28.72 | 40.24 | 10.31 | - |
| 9260 | 2315 | k4 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 9308 | 2327 | b0 | 1: A 1 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 9308 | 2327 | b1 | 1: A 1 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 9308 | 2327 | k4 | 18: A 18 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 9356 | 2339 | b1 | 3: A 2, alive 1 | G'=Z167xZ2, q=7, s=3, f=2 U=<2> (order 83) on axis 0, r=24 (3) | 16.77 | -55.41 | 90.28 | G' | 32.71 | 10.67 | 67.51 | -38.64 | - |
| 9356 | 2339 | k4 | 63: A 57, alive* 6 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 9388 | 2347 | b1 | 10: A 1, B 1, alive* 5, alive 3 | G'=Z17xZ3xZ2, q=23, s=2, f=11 U=<2> (order 8) on axis 0, r=54 (3) | 8.8† | -15.65 | 75.65 | L34<(0, 0, 1),(1, 0, 0)> | 21.54 | 24.85 | 34.64 | -6.85 | - |
| 9388 | 2347 | k4 | 5: A 5 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 9484 | 2371 | b1 | 18: A 10, alive* 5, alive 3 | G'=Z5xZ3xZ2, q=79, s=6, f=13 U=<2> (order 4) on axis 0, r=84 (2) | 15.71† | 13.14 | 189.75 | L2<(0, 0, 1)> | 42.49 | 41.64 | 67.56 | 28.85 | - |
| 9484 | 2371 | k4 | 59: A 55, alive* 4 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 9532 | 2383 | b1 | 16: A 10, B 4, Z2 1, alive 1 | alive slices not estimated (sGeff > 101 [E: C2 > 150], or an estimator cap) | | | | | | | | | |
| 9532 | 2383 | k4 | 6: alive* 6 | G'=Z4xZ7xZ5, q=17, s=2, f=8 U=<3> (order 6) on axis 1, r=120 (2) | 93.84† | -70.19 | 194.19 | L5<(0, 0, 1),(0, 0, 2)> | 70.8 | 61.7 | 102.5 | 23.65 | - |
| 9596 | 2399 | b1 | 35: A 17, B 2, alive* 5, alive 11 | G'=Z11xZ2, q=109, s=12, f=9 U=<2> (order 10) on axis 0, r=52 (1) | -277.21 | 25.9 | 129.49 | L2<(0, 1)> | 23.03 | 33.45 | 39.81 | -251.31 | - |
| 9596 | 2399 | k4 | 10: A 10 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 9692 | 2423 | b1 | 8: A 6, alive 2 | G'=Z173xZ2, q=7, s=3, f=2 U=<6> (order 43) on axis 0, r=40 (3) | 40.58† | -58.13 | 116.77 | G' | 40.39 | 30.99 | 71.95 | -17.55 | - |
| 9692 | 2423 | k4 | 8: A 2, alive* 4, alive 2 | G'=Z4xZ121, q=5, s=2, f=2 U=<4> (order 55) on axis 1, r=60 (3) | -532.93† | -82.27 | 120.27 | G' | 70.56 | 18.26 | 92.14 | -615.2 | - |
| 9716 | 2429 | b0 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 9716 | 2429 | b1 | 2: A 2 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |
| 9836 | 2459 | b1 | 0: - | no GP slice with 24 <= r <= 160 | | | | | | | | | |
| 9836 | 2459 | k4 | 17: A 17 | **every slice dead at an exact gate [P, M]** | | | | | | | | | |

- In the gate column, A, B, Z2, w and S1 count the slices killed by each gate; "alive" counts slices passing every
  gate run, and "alive*" those with B capped.
- b0 lines at prime v are omitted: a prime v has no staged b0 slice, see b0direct (§5).
<!-- /TABLE -->

## 4. Step 5: exact stage-1 decisions [M]

`theory/open_screen_stage1_exact.py` decides stage 1 with up to three methods:
- **A**, the repository's enumerator `ab_b0_stage1.stage1`;
- **B**, a DFS plus pair join that shares no code with `ab_lib` and derives G′, the U-orbits, the letters
  α + f·σ and T₁ from the definitions; its set-up is checked against `ab_lib`;
- **C**, CP-SAT on the integer model; INFEASIBLE decides, UNKNOWN does not.

Every run includes planted controls: a random U-invariant y-quadruple whose energy is at most the real T₁(0) gives
the target, and every method must recover it.

<!-- STEP5 -->
The candidates were taken in order of estimated cost among the tier-0 slices, after the ℤ₂ and w-screens. They are
listed in the order decided. Logs: `logs/open_screen/stage1_exact.jsonl` (1719) and
`logs/open_screen/stage1_exact_step5.jsonl` (the rest).

| # | 4v (v) | slice | cost [E] | method A | method B | method C (CP-SAT) | planted | verdict |
|---|---|---|---|---|---|---|---|---|
| 1 | 6876 (1719) | b0 ℤ₉ × ℤ₁₉₁, s = 10, f = 19, U = −1 | 42.8 | 0 solutions (944 blocks) | 0 (1888 blocks) | INFEASIBLE (4.8 s) | A 4/4, B 4/4, C 3/3 | **empty [M]** |
| 1′ | 6876 (1719) | b0 ℤ₃ × ℤ₃ × ℤ₁₉₁, s = 10, U = −1 on either ℤ₃ | – | – | 0 (10,944 blocks) | INFEASIBLE (89 s, 8 s) | B 3/3, C 3/3 | **empty [M]** (both) |
| 2 | 3788 (947) | b1 ℤ₁₁ × ℤ₂ × ℤ₄₃, s = 6, U order 5 | 26.0 | – | 0 (832 blocks) | INFEASIBLE (1.9 s, strong model) | B 3/3 (and 3/5), C 3/3 | **empty [M]** |
| 3 | 5036 (1259) | b1 ℤ₁₇ × ℤ₂ × ℤ₃₇, s = 4, U order 4 | 24.8 | – | 0 (45,318 blocks; join 1.4·10⁷ × 3.2·10⁷) | INFEASIBLE (121 s, strong model) | B 2/3, C 2/3 (rest undecided) | **empty [M]** |
| 4 | 5804 (1451) | b1 ℤ₂₅ × ℤ₂ × ℤ₂₉, s = 4, U order 10 | 27.0 | – | 0 (402,884 blocks; join 1.9·10⁷ × 3.7·10⁷) | UNKNOWN at 3600 s | B 3/3 | **empty by one exhaustive method [M]**; second method undecided |
| 5 | 3452 (863) | k4 ℤ₄ × ℤ₅ × ℤ₄₃, s = 6, U = −1 on ℤ₄ × ℤ₅ | 35.1 | – | <!-- K863B -->see below<!-- /K863B --> | UNKNOWN at 3600 s (planted also UNKNOWN) | <!-- K863P -->–<!-- /K863P --> | <!-- K863V -->see below<!-- /K863V --> |

**Pruning by exact consequences [P]** (`--strong`). Each is necessary for a stage-1 solution, so it cannot lose
solutions:
- **Row sums.** Σ_k r_k² = Σ_h T₁(h) (the trivial character), so each block has r_k² ≤ Σ T₁.
- **Real characters.** For every real character χ of G′ constant on the U-orbits, Σ_k A_k(χ)² = Σ_h χ(h)T₁(h), so
  each block has A_k(χ)² at most that.
- **Symmetry breaking (C).** Negating or permuting blocks maps solutions to solutions, so C imposes r_k ≥ 0 and
  r₀ ≥ … ≥ r₃.
- **Bucketing (B).** B joins pairs only within matching buckets of the pair-sum invariants. These are the lag-0
  energy and the A(χ)², plus, for k4, the Q-level PAF; all are linear in the block PAF and fixed by the target.
- **Compression normal form (k4 with Q = ℤ₄ and |K| odd).** By the OpenAI model's theorem
  (`report/openai_4p3_review.md` §1), every block's Q-compression is ±(a, a, a, a) with a odd, or ±(t·1 + 2δ_s)
  with t odd, under the energy bound. That is 134 vectors at |K| = 215.
- **k4 planted controls.** Here the planted quadruples are four real enumerated blocks (one constant compression,
  three singletons), so that the normal-form pruning keeps them.

**Re-reading the deficit at real stage-1 solutions** (prompt step 5) was not possible: no candidate has a stage-1
solution.
<!-- /STEP5 -->

## 5. The other families

<!-- OTHER -->
- **Direct (unstaged) b0 multiplier slices at the 64 prime v** (`scripts/b0_open_screen.py`, untwisted and 3-cycle
  twist, 15–130 orbits; `logs/open_screen/b0_prime_direct.jsonl`).
  - 93 slices at 33 primes. The other 31 primes have no odd d | v − 1 in range.
  - Only 8 primes have a row-sum-feasible slice. The cheapest are:
    - 1103 (d = 19, 59 orbits): n = 157.5;
    - 1783 (d = 27, 3-cycle): n = 190.8;
    - 739 (d = 9, 3-cycle): n = 209.5;
    - 1291 (d = 15): n = 238.1.
  - The slice that built H(3532) has n = 111.9, and the twisted 919 slice n = 91.1. Every open prime v is 45 bits or
    more beyond them [E].
- **Cooper–Wallis: T(t) × Williamson-type of order w = v/t** (`scripts/tmat_primev.py`, `logs/open_screen/tmat_cw.jsonl`,
  t = the prime cofactors ≥ 167).
  - Row-sum feasible:
    - **T(463), k = 11 (43 orbits): n = 101.5, log₂E = 27.5**, the route to H(9260) with Williamson of order 5.
    - T(223), k = 3 (156.8), for 669 and 1115.
    - T(491), k = 7 (169.1), for 1473.
    - T(283), k = 3 (207.1), for 1981.
    - T(631), k = 7 (225.2), for 1893.
    - T(571), k = 5 (280.5), for 1713.
    - T(643), k = 3 (534.2), for 1929.
  - T(167), T(179), T(191), T(347), T(359), T(419), T(443) and T(479) have no row-sum-feasible odd-k slice [P, M].
  - T(599) and T(719) (for 1198 and 1438) have none, by `report/tmat_routes.md` §10.3 [P].
  - The T-route scale (T(131) at n = 51 took about 2^27; 883 / k = 21 at n ≈ 110 is about 2^69) puts T(463) / k = 11
    at about 2^63 [E]. That is far over the budget.
  - The other split, Williamson-type matrices of the large order, is not screened.
- **ℤ_q all-character family** (Kumar's 2148 family; from `report/zq_census_answer.md` §3, not rerun).
  - Positive-estimate survivors: 669 ℤ₃ (S6 r118, 2^111.4), 1969 ℤ₁₁ (S22 r110, 2^23.7), 2327 ℤ₁₃ (S33 r90, 2^35.0)
    and 1115 ℤ₅ (S18 r90, 2^1.6).
  - All are rich hosts (r ≥ 90) whose next levels exceed the list caps or 600 s of CP-SAT; they are not staged slices
    in the sense of this screen.
  - 1169, 1437, 1473, 1795, 1893, 1929, 2171 and 2429 are closed there [P, M].
<!-- /OTHER -->

## 6. Verdict on search

<!-- VERDICT -->
- **No named target for a staged stage-2 solver comes out of this screen.**
  - The five cheapest tier-0 candidates were taken to an exact stage-1 decision. Four are empty, one of them (1451)
    by a single exhaustive method. The fifth (863 k4) is <!-- K863W -->as recorded in §4<!-- /K863W -->.
  - Together with the earlier kills of 2396 and 3788 s = 2, every candidate decided so far has a positive estimate
    and an empty stage 1. This matches the repository's calibration warning.
- **What remains alive [E only]:**
  - the tier-0 slices not run: 4156 k4, 5708, 5036 s = 12, 7572, 5692, 9260, 3788 s = 14 and 4748 k4;
  - the 2396, 5948, 7324, 7948 and 9596 slices of tier 1;
  - the tier-2 slices, whose expected count is positive but whose split misses the regime.

  Most of these tier-0 estimates rest on the optimistic energy-only stage-1 count. All of them have U ≠ 1 and r = 44–90.
- **Next step, if the route is kept.**
  - Run `theory/open_screen_stage1_exact.py KEY --methods B,C --strong` on the remaining tier-0 slices, about
    1–60 minutes each with method B. Treat any that is non-empty as a real candidate: re-read its deficit at the real
    stage-1 solutions with `analyse(..., real=True)`.
  - Only then write a solver, validated first on a same-shape known order.
  - The k4 slices (4156, 4748, 3452) have the border by Theorem B, so a core there would complete at once.
- **What the evidence says about the route.** The cheap staged slices sit at the orders whose v − 1 or v − 3 has a
  mid-size prime (q ≈ 29–89) with a fitting s. Stage 1 of every such slice decided exactly here is empty.
<!-- /VERDICT -->

## 7. Cut-off and what was not covered

<!-- CUTOFF -->
- **Run time.** 2026-11-05 22:20 to 2026-11-06 (see the last commit); about 13 hours of the 30-hour box. It stopped
  because the prompt's step 5 has a cap of 5 candidates, and those were taken to a decision or to an exact undecided
  state.
- **Not covered:**
  - diagonal multipliers acting on several axes with one generator (the 647 type);
  - the families k6 and k8 (border-of-width 24 / 32 Theorem 7 variants);
  - k4 groups that do not split as Q × K (for example ℤ₈ × odd with K ⊃ 2ℤ₈): `ab_sweep.configs_bge3` builds
    G = Q × H only;
  - slices with r < 24 or r > 160;
  - sGeff > 101 (C₂ > 150 by the measured ratio, so outside the regime) [E];
  - B undecided for 699 alive slices (the α-key cap);
  - 26 cost timeouts and 1 cost memory error;
  - the Williamson-type side of Cooper–Wallis;
  - a rerun of the ℤ_q census;
  - exact stage 1 for the tier-0 slices beyond the five.
- **Caps that bind, by order:** see the gate column of the full table ("alive*" = B capped). For example, the 2396
  b1 slices are mostly B-capped.
<!-- /CUTOFF -->

## 8. Reproduce

```
python3 theory/open_screen_enum.py --list                      # 3242 slices
logs/open_screen/run_screens.sh                                # screens A, B (4 shards, restart-safe)
OS_MEMGB=2.8 logs/open_screen/run_deficit.sh 4                 # costs of the survivors with sGeff <= 80
python3 theory/open_screen_deficit.py --controls               # calibration (933, 523, 863, 235)
python3 scripts/b0_open_screen.py $(cat logs/open_screen/prime_v.txt) 20000      # direct b0 at prime v
python3 scripts/tmat_primev.py 167,179,191,223,283,347,359,419,443,463,479,491,571,631,643 20000   # T(t)
python3 theory/open_screen_stage1_exact.py 'KEY' --methods A,B,C --planted 4   # step 5
python3 theory/open_screen_table.py --md TABLE.md               # table.jsonl and the markdown table
```
