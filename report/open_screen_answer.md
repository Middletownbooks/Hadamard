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
(filled at the cutoff)
<!-- /BOTTOM -->

---

## 1. What was screened

<!-- SCOPE -->
(filled at the cutoff)
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
(filled at the cutoff)
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
(filled at the cutoff)
<!-- /STEP5 -->

## 5. The other families

<!-- OTHER -->
(filled at the cutoff)
<!-- /OTHER -->

## 6. Verdict on search

<!-- VERDICT -->
(filled at the cutoff)
<!-- /VERDICT -->

## 7. Cut-off and what was not covered

<!-- CUTOFF -->
(filled at the cutoff)
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
