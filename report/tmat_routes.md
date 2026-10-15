# Open orders 2000–3000 after the 2060 construction: what carries over

Date: 2026-10-06. Scripts: `cw_build.py`, `tmat_estimate.py`, `tmat_core_check.py`, `tmat_search.py` (draft, not run).
Tags: [M] measured, [E] estimate.

## 1. What 2060 is

Schneider's order-2060 matrix is a Cooper–Wallis construction: 4·t·w with t = 103 and w = 5.
- **T-matrices of order t = 103:** four disjoint ternary circulants with ΣPAF = t·δ, invariant under the multiplier
  subgroup K = {1, 46, 56} of order 3. That gives 35 orbits, one label in 0..7 per orbit.
- **A Williamson quadruple of order w = 5.**

**Direct consequence [M]:** 721 = 7·103 gives a Hadamard matrix of order **2884**, verified (`cw/RESULT_2884.md`).
No other open v in 2000–3000 is a multiple of 103.

## 2. Key equivalence (makes the repo's machinery reusable)

T-matrices of order t ⇔ four ±1 circulants A, B, C, D of length t with
- ΣPAF = 0 off zero, i.e. a **b = 0 Goethals–Seidel core**, and
- **a_j b_j c_j d_j = +1 for every j** (the Williamson product condition with ε = +1).

The map is (A, B, C, D) = H₄·(T₁, …, T₄). The 8 sign patterns with product +1 are exactly the 8 orbit labels.

**Consequences:**
- **For prime v, the T-matrix route (w = 1) is a special case of the b = 0 route and buys nothing.**
  - [E] At v = 523, d = 9 the T-matrix model needs n = 141.9 constraint bits; the repo's b = 0 estimate is n ≈ 142.
    The product condition removes 2^m options and about as many solutions.
  - So 523, 571, 599, 643, 647, 659, 719 and 739 are untouched by this method.
- **Only composite v = t·w with w ∈ {3, 5, 7, …} gain**, because Cooper–Wallis multiplies the order. In 2000–3000
  those are **573 = 3·191, 669 = 3·223 and 537 = 3·179**.
- **Symmetry must have odd order [M].** If −1 ∈ K, every row is symmetric and each lag-class sum is forced odd, so
  zero is off the lattice.

## 3. Feasibility screen [E] (Gaussian, `tmat_estimate.py`, 2·10⁵ samples)

n = 3m − log2 E is the number of constraint bits (comparable to the repo's `nbits` n). "Per class" is log2 of the
expected solutions modulo the symmetry group 4!·2⁴·(t−1)/k.

| target v (order) | t | k | orbits m | n (bits): tmat_estimate / nbits.py | log2 E per class | status |
|---|---|---|---|---|---|---|
| (control) 515 (2060) | 103 | 3 | 35 | **62.9 / 62.4** | 28.5 | solved by Schneider; blind control running |
| **573 (2292)** | 191 | 5 | 39 | **79.2 / 80.2** | 24.0 | **best new target: about 17 bits above the solved control** |
| 573 | 191 | 19 | 11 | 22.3 | −1.2 | **decided [M]: no solutions** (SAT, all 20 row-sum cases UNSAT, 4 min) |
| 669 (2676) | 223 | 3 | 75 | 156.8 / 156.6 | 53.4 | plentiful but hard (the 523-b=0 level) |
| 669 | 223 | 37 | 7 | 15.1 | −5.2 | **decided [M]: no solutions** (SAT, all 36 row-sum cases UNSAT, 7.5 min) |
| 537 (2148) | 179 | – | – | – | – | (Z/179)* has no useful odd subgroup (178 = 2·89); unrestricted T(179) is harder than 716 was |

**Prime v (w = 1), reference only.** These are identical to the b = 0 route:
- 523/9: 141.9.
- 571/5: 281.
- 659/7: 237.
- 739/9: 210.

**Correction (same day).** The first version of the estimator did not check that row sums can be realised. Each block
sum is ±1 + k·j with j ≡ m − 1 (mod 2). Without that check it reported 599/23 (66 bits), 571/19 (75), 571/15 (95),
647/19 (87) and 647/17 as cheap. **All five are infeasible**: no row-sum quadruple with a² + b² + c² + d² = 4t is
reachable, which `nbits.py` confirms (n = ∞). The check is now in `tmat_estimate.py` (`rowsum_feasible`).

**The six known b = 0 cores (223, 283, 311, 463, 487, 491) are not in product form [M].** This was checked under all
per-block shifts, reversals and negations (`tmat_core_check.py`). The checker's positive control, a scrambled T(103),
gives 2 hits. So the Alpöge 223 core does not give 669 for free.

## 4. Repo components that transfer (shortcuts)

1. **`wil_sat.py` (SAT + CAS).** It already encodes the product condition as 4-literal parity clauses, with an
   IPASIR-UP PSD propagator.
   - For T-matrices: 4 variables per K-orbit with a XOR-zero constraint (ε = +1, no exceptional position), ΣPAF = 0
     per lag class as cardinality constraints, and the PSD bound |X̂(ω)|² ≤ 4t per block and per pair.
   - Unlike tabu, it is complete for the ansatz: it can also prove a (t, k) empty.
   - The k = 19 and k = 37 slices (m = 11 and 7) can be decided outright in seconds.
2. **`redlib.py`, `engine.py` and `kt_gen.py`.** The reduced U-orbit b = 0 setting at prime length with a subgroup of
   order d is exactly this setting. The incremental orbit-level PAF tabu applies with moves restricted to
   product-preserving relabels (pair-flips). `tmat_search.py` is a standalone draft of that.
3. **`nbits.py` / `b1_multiplier_enum.py`.** These use the same Gaussian costing. Their row-sum conditioning
   (a² + b² + c² + d² = 4t) can sharpen the estimate above.
4. **`kt_join.py` / `kt_generic.py` (k-tree join).** These do not transfer cleanly: the product condition couples the
   four lists pointwise (D = A∘B∘C). A per-pattern 2-list split exists, but it costs about 2^78 at 191/k5, which is
   worse than direct search.
5. **Process.** Calibrate blind on a known-positive instance before the target: T(103, k = 3) here, as v = 487 was
   for the k-tree and the W0/W2 controls were for Williamson.

## 5. Recommended order of work

1. **Decide the tiny slices exactly:** 191/k = 19 and 223/k = 37 (and 523/k = 29 for completeness), using SAT or plain
   enumeration of 8^m.
2. **Positive control:** find T(103, k = 3) blind with (a) tabu and (b) SAT + PSD, and measure time-to-solution.
3. **Target 191/k = 5 (→ 2292)** with the faster method. Extrapolate from the control: 16 extra bits means up to
   about 2^16 times the work for brute-force-like methods, less if pruning scales better.
4. **669 (223/k = 3)** only if step 3 shows good scaling. At 157 bits it is in the range that has resisted the repo's
   523 work.

## 6. Results so far (2026-10-06)

- **Step 1, small slices decided [M]:**
  - 191/k = 19: no solutions (20 cases UNSAT).
  - 223/k = 37: no solutions (36 cases UNSAT).
  - 571/k = 57: 44 of 58 cases UNSAT, then stopped. It needed 9.8 GB, apparently a memory leak across cases, so it was
    deferred.
  - 523/k = 29: deferred.
- **T(103, k = 3) positive control: PASS for the k-tree method [M].**
  - Method (`tmat_kt.py`):
    1. Fix the support S = supp T₁ ∪ supp T₂ with |S| = ⌊m/2⌋ orbits.
    2. Write T₁ and T₂ as (P ± R)/2 on S, and likewise on Sᶜ.
    3. This gives **four independent lists** with PAF(P) + PAF(R) + PAF(P′) + PAF(R′) = 2t·δ, solved by the repo's
       k-tree join (`kt_generic.kt_join4`).
  - Sanity run with the known solution's S: found in 62 s.
  - **Blind run:** a random S (seed 7) gave a solution on the **first S tried, in 720 s**. It differs from Schneider's
    array. Built into an order-2060 matrix, it verifies H Hᵀ = 2060·I (`cw/T103_blind.json`).
  - Comparison on the same instance:
    - Tabu (`tmat_search.py tabu2`, about 42k it/s): best E = 2–3 after 141M iterations (56 min), no solution.
    - SAT + PSD (`tmat_sat.py`): one row-sum case UNSAT in 211 s, the next still unfinished after 2 h.
  - So the support split with the k-tree is the method of choice.
- **Now running:** T(191, k = 5), the route to 573 (order 2292), with 4 seeds of the k-tree method.
  - Lists are about 2^19 per side.
  - Expected solutions per random S are about 2^−1.2 [E] (2^37.8 raw / 2^39 supports), against about 2^7 at
    t = 103, so many S will be needed.

## 7. Correction: T(191, k = 5) is empty, so 573 has no route here (2026-10-06, later)

A second, independent analysis by another LLM pointed this out. I verified it before acting [M].

**The row-sum condition.**
- In T-form each row sum is s_i = T_i(0) + k·(integer), and exactly one row contains position 0.
- So **three of the s_i are ≡ 0 (mod k) and one is ≡ ±1 (mod k), with Σ s_i² = t.**
- §3 and §6 checked only the GS-form condition (each block sum ±1 mod k), which is necessary but not sufficient.

**Exact enumeration of this condition over all targets.** Viable (t, k) pairs with k odd and at least 4 orbits:
- (103, 3): the known solution's sums 7, 3, 3, 6 have residues 1, 0, 0, 0;
- (223, 3);
- (523, 3) and (523, 9);
- (571, 3) and (571, 5);
- (643, 3);
- (659, 7);
- (739, 3).

**Dead:**
- 191 (k = 5 and 19), 599, 647, 179 and 719: no k gives a viable pair.
- The remaining k at 103, 223, 523, 571, 643 and 739.

This matches the other analysis exactly.

**Consequences:**
- **573 has no T(191) route.** The T(191, k = 5) k-tree runs were stopped after about 1 h, having finished no S.
- 191/k = 19 and 223/k = 37 are empty for this reason too, consistent with the SAT results in §6.
- `tmat_estimate.py` now applies the exact condition, and its `rowsum_feasible` flag returns False for 191/5.

**Independent confirmations from the other analysis:**
- H(2884) was found independently, with a different order-7 Williamson quadruple. I checked that quadruple: it is
  symmetric with ΣPAF = 28·δ.
- The T-ansatz densities agree with ours.
- The Alpöge cores are not in T-form.
- Plain annealing cannot solve T(103).

**Remaining composite target: 669 via T(223, k = 3).** Raw E is about 2^68, density about 2^-157.
- **k-tree cost [E]:** 75 orbits, so a balanced S has about 37 orbits and each side list is about 2^36 entries
  (2^17 at t = 103). Solutions per random S are about 2^68.2 / 2^75 ≈ 2^-7. A complete join on one S costs about 2^72
  pair operations, and reaching 157 bits with a 4-list join needs lists of about 2^52, against the 2^36 one S can supply.
- **Kaggle verdict:** 5 notebooks × 4 cores is roughly 2^5–2^8 more compute than here. That is nowhere near the gap,
  so Kaggle is not useful for 669 with the current ansatz.
- **What could change that:** extra structure that removes about 70 bits. The untested lead from both analyses is a
  *twisted* symmetry: multipliers that also permute or negate the four rows (for example K of order 6 at t = 223,
  with x ↦ −x acting as a row swap). This changes both the orbit count and the row-sum conditions. It is cheap to
  screen locally with the estimator, plus exact SAT on small slices.
  - If a twisted ansatz brings 669 to about 80 bits or less, the per-S k-tree join is embarrassingly parallel across
    S and would suit Kaggle.

## 8. Twisted-symmetry screen (2026-10-06) [M for feasibility and validation; E for log₂E]

**Ansatz** (`scripts/tmat_twist.py`): U = ⟨u⟩ ≤ (Z/t)* of order d, and g = (π, ε) ∈ B₄, a signed permutation of the four
rows with g^d = 1, satisfying T_{π(i)}(u·x) = ε_i·T_i(x).
- Labels propagate along each U-orbit by g, and position 0 must carry a g-fixed label.
- The ΣPAF condition is invariant under this action, so there is one equation per class of ⟨U, −1⟩ on the nonzero
  lags.
- Twisting allows even d, where −1 ∈ U acts as a row swap or negation.
- One representative of each B₄ conjugacy class (signed cycle type) was screened for every d | t − 1 with at least 3
  orbits. For each: exact row-sum feasibility (an ILP over orbit contribution vectors, against every s with Σs² = t),
  then the Gaussian log₂E.

**Validation against brute force** (`twist_validate_final.txt`):
- **Infeasibility claims are sound.** Every (t, d, g) marked infeasible has brute-force count 0. This covers 59
  combinations at t = 13–43.
- **Estimates where solutions exist:**

  | (t, d, g) | brute force | estimate |
  |---|---|---|
  | (19, 3, untwisted) | 2^14.6 | 2^13.0 |
  | (19, 3, 3-cycle twist) | 2^12.3 | 2^11.5 |

  So the estimator underestimates by about 1–1.5 bits.

**Results.** All 34 (d, g) combinations at t = 223 were screened, plus 38 at t = 191 and t = 179 together.

| t (target) | feasible (d, g) | log₂E | n (bits) | reading |
|---|---|---|---|---|
| 223 (669) | d = 3, g = id | 68.2 | 156.8 | unchanged |
| 223 (669) | d = 3, g = 3-cycle on rows | 66.5 | 156.5 | twist gains nothing |
| 223 (669) | d = 6, g = row swap | −62.0 | 175.0 | empty: −1 ∈ U gives 37 classes for 38 orbits |
| 223 (669) | d = 2, g = row swap | −189.4 | 524.4 | empty |
| 191 (573) | d = 2, g = row swap | −151.3 | 438.3 | empty; every d = 5, 10, 19 and 38 twist is row-sum infeasible |
| 179 (537) | d = 2, two classes of g | about −138 | about 406 | empty |

All other combinations fail row sums or have no g-fixed label for position 0.

**Conclusions:**
- **The twisted-symmetry lead is closed for 669, 573 and 537.**
- When −1 ∈ U, the lag classes no longer pair as {s, −s}, so halving the orbits costs as many constraints as it saves.
- The order-3 row twist leaves the orbit count unchanged.
- **669 stays at about 156 bits.** That is about 94 bits beyond the solved T(103) and out of reach of this method on
  this hardware or on Kaggle.
- **Superseded for −1 ∈ U by §10 [P].** Every slice above with even d is exactly empty by the reflection lemma of §10.
  There the twisted Gaussian is unreliable: it predicts 2^0.9–2^4.4 solutions on small cases whose exact count is 0.

## 9. 3000–4000: T(131, k = 5) gives H(3668) (2026-10-15)

The same route reaches **v = 917 = 7·131**. Full record: `report/b0_917_answer.md`.
- **The ansatz.** K is the order-5 subgroup {1, 53, 58, 61, 89} of F₁₃₁^*: 27 orbits, 2⁸¹ labellings.
- **Row sums [P].** The §7 condition is met by (9, 5, 5, 0), and only by it; k = 13 and k = 65 are infeasible.
- **Estimate [E].** log₂E = 30.0 (16.7 per symmetry class).
- **Search [M].** `scripts/tmat_kt.py 131 5` found T(131) on the second random support, 23 s on one core. A second seed
  found another on the eighth support (118 s).
- **Census [M].** Exact complete joins over 40 random supports gave 0.70 solutions per support against 0.5 predicted.
- **Same-shape ladder [M].** T(41), T(61) and T(101) at k = 5, with exact totals at t = 41 and 61 (the estimator
  0.5–0.9 bits low); each built H(28t).
- **Construction [M].** Cooper–Wallis with the order-7 Williamson quadruple (as for H(2884)) gives four ±1 circulants
  of length 917 with Σ PAF = 3668·δ.
- **Verification [M].** The independent `scripts/verify_h3668.py` checks H Hᵀ = 3668·I in integer arithmetic. The
  bundle is `pub/H3668/`.
- **Other composite open orders 3000–4000 [P].**
  - 955 (= 5·191), 537 and 789 have no row-sum-feasible odd k.
  - 669 is unchanged (k = 3, about 157 bits).

## 10. Every open order below 4000, screened for a new T-matrix route (2026-10-15) [P, M, E]

The question: after T(103) → H(2060) and T(131) → H(3668), can a new T-matrix reach any other open v?

Scripts: `tmat_primev.py`, `tmat_twist3_open.py`, `tmat_reflection_check.py`, `tmat_twist_d2_cal.py`,
`tmat_census_2k2.py`, `tmat_twist3_census.py` and `tmat_menu.py`. Logs are in `logs/tmat/`.

**Bottom line.**
- **No open order below 4000 has a T-matrix route within reach.**
- **537, 789 and 955 are closed exactly [P].** No T(179), T(263) or T(191) has a nontrivial multiplier symmetry,
  twisted or not.
- **The cheapest surviving slices** are 919 (T(919) with a 3-cycle twist, d = 27, n = 91.1) and 883 (T(883),
  k = 21, n ≈ 110). Each costs about 2^67–2^69 operations, against about 2^27 for T(131). Next is 669 (T(223, k = 3),
  n = 156.8).
- For prime v a T(v) slice is the b = 0 slice of the same symmetry with the product condition added (§2). The b = 0
  slice therefore has at least as many solutions.
- **The findable T-matrices miss the open list.** Every T(t) now within reach (n ≤ 81, primes 100–420) gives
  Cooper–Wallis orders that are not open below 4000, so they matter only above 4000.

### 10.1 What a T-route needs

Cooper–Wallis turns T(t) and Williamson-type matrices of order w into H(4tw).
- **Composite open v:** 537 = 3·179, 669 = 3·223, 789 = 3·263 and 955 = 5·191 (917 and 933 are built). Since 3 and 5
  are Williamson orders, these need T(179), T(223), T(263) and T(191).
  - The other split needs Williamson-type matrices of order 179–263. That is not screened here.
- **Prime open v:** only w = 1, i.e. T(v) itself.
- All of these t are ≡ 3 (mod 4).

**The symmetric ansätze are exhaustive [P].**
- For prime t, a symmetry with a fixed point is, after a common shift, x ↦ ux combined with a signed row permutation
  g (`tmat_twist.py`). Here U = ⟨u⟩ is cyclic of order d and g^d = 1.
- **d odd:** g has order 1 (the K-invariant ansatz, k = d) or 3. A 3-cycle needs 3 | d and fixes one row. Its row
  sums are s₁ = s₂ = s₃ = (d/3)·a and s₄ = ±1 + d·b, with 3s₁² + s₄² = t.
- **d even:** −1 ∈ U, and the lemma below settles it.

### 10.2 Reflection lemma [P]

This extends §2's remark on −1 ∈ K to every twist.

Let t be an odd prime and T a T-matrix of order t with T_{π(i)}(−x) = ε_i·T_i(x) for an involution (π, ε) of B₄.
1. **Parity.**
   - If π fixes row i, then PAF_i(c) ≡ [c/2 ∈ supp T_i] (mod 2). The map x ↦ −x − c pairs equal terms of PAF_i(c),
     and its one fixed point x = −c/2 contributes ±T_i(c/2)².
   - If π swaps rows i and j, then PAF_i = PAF_j.
   - So Σ_i PAF_i(c) ≡ [c/2 lies in a π-fixed row] (mod 2). For this to vanish at every c ≠ 0, the π-fixed rows must
     carry no nonzero position.
2. **Shape.** Position 0 is fixed by x ↦ −x, so it lies in a π-fixed row. Hence π is a transposition (i j), the two
   fixed rows are ±δ₀ and 0, and T_j(x) = ε·T_i(−x).
3. **Weighing sequence.** The pairs {x, −x} split between rows i and j. So W = T_i lives on a half-system, has weight
   (t − 1)/2 and PAF(W) = ((t − 1)/2)·δ: W is a circulant weighing sequence. In particular (t − 1)/2 = s².
4. **Self-conjugacy, for t ≡ 3 (mod 4), t > 3.**
   - s is odd and 3 | s (otherwise 3 | 2s² + 1 = t). So t ≡ 1 (mod 3) and (3/t) = −(t/3) = −1.
   - Hence −1 ∈ ⟨3⟩ mod t, and complex conjugation fixes every prime of Z[ζ_t] above 3.
   - From W(ζ)·W(ζ⁻¹) = s², each such prime divides W(ζ) exactly v₃(s) times, so 3 | W(ζ).
   - Then all coefficients of W are congruent mod 3. That is impossible for a sequence with both zero and nonzero
     entries.

**So for prime t ≡ 3 (mod 4), t > 3, no T-matrix of order t is invariant under x ↦ −x, even up to a signed row
permutation [P].**
- This covers every even d, twisted or not, at every t in this section.
- §8's even-d rows become exact zeros.
- For t ≡ 1 (mod 4), steps 1–3 still force t = 2s² + 1 and a half-system CW(t, s²). For each such prime with
  1 < s < 40 (t ≤ 2593), some prime factor of s is self-conjugate mod t, so the class is empty there too [M].

**Checks [M]** (`logs/tmat/reflection_check.jsonl`):
- **Orders 3, 5, 7.** All T-matrices (192, 3,840 and 37,632) were tested against the 76 elements of order ≤ 2 in
  B₄.
  - At t = 3 all 192 are reflection-invariant, with row weights (0, 1, 1, 1) as step 2 predicts.
  - At t = 5 and 7 none are.
- **Four symmetric rows:** 0 solutions at t = 11 and 13 (all 8^{(t+1)/2} labellings).
- **Half-system CW(19, 9):** 0 of 2^18. 19 = 2·3² + 1 is the first case that step 3 leaves open.
- **The twisted Gaussian fails here.**
  - At t = 7 and 11 it predicts 2^3.2–2^4.4 solutions for the transposition classes, where the exact count is 0
    (`logs/tmat/twist_d2_calibration_final.txt`).
  - It gives 2^0.9 and 2^1.4 at (19, d = 6), where §8's own brute force (`twist_validate_final.txt`) counts 0.
  - It misses the joint parity condition of step 1.

### 10.3 Per-order results

- n = 3m − log₂E constraint bits [E].
- Obstructions marked [P] are exact: row sums, or the lemma.
- For scale: T(131, k = 5) has n = 51.0 and T(103, k = 3) has n = 62.9.

| v | route | surviving symmetric slices (orbits) | n (bits) |
|---|---|---|---|
| 537 | T(179) | none: d = 89 leaves 3 orbits, 3 ∤ 178, even d by the lemma [P] | – |
| 789 | T(263) | none: d = 131 leaves 3 orbits, 3 ∤ 262, even d by the lemma [P] | – |
| 955 | T(191) | none: k = 5, 19, 95 fail row sums (§7), 3 ∤ 190, even d by the lemma [P] | – |
| 669 | T(223) | k = 3 (75); 3-cycle d = 3 (75) | 156.8; 156.5 |
| 919 | T(919) | **3-cycle d = 27 (35; row sums (9, 9, 9, ±26))**; 3-cycle d = 9 (103); k = 9 (103); k = 3 (307) | **91.1**; 268.6; 269.1; 803.3 |
| 883 | T(883) | **k = 21 (43; row sums (±1, ±21, ±21, 0), 883 = 2·21² + 1)**; k = 7 (127); k = 3 (295) | **111.9**; 330.0; 767.6 |
| 739 | T(739) | 3-cycle d = 9 (83; untwisted k = 9 fails row sums); k = 3 (247) | 209.7; 626.6 |
| 659 | T(659) | k = 7 (95) | 236.6 |
| 571 | T(571) | k = 5 (115); k = 3 (191) | 280.5; 466.4 |
| 991 | T(991) | k = 9 (111); k = 5 (199); k = 3 (331) | 293.1; 525.8; 875.3 |
| 751 | T(751) | k = 5 (151); k = 3 (251) | 383.6; 638.2 |
| 971 | T(971) | k = 5 (195) | 513.8 |
| 643, 823, 859, 907 | T(v) | k = 3 only | 534.3, 708.5, 743.8, 791.4 |
| 599, 719, 839, 863, 947 | T(v) | none [P] | – |

- 3-cycle twists with d = 3 pass row sums at the prime v. They have the same 190–330 orbits as k = 3 and were not
  estimated.
- With no symmetry at all, T(179), T(191) and T(263) cost 362, 391 and 569 bits [E]
  (`logs/tmat/asym_k1_179_191_263.jsonl`).

### 10.4 Calibration of the two cheapest slices [M]

**883/k = 21** belongs to the family t = 2k² + 1, with row sums (1, k, k, 0) and 2k + 1 orbits.
- **t = 19, k = 3:** the exhaustive count is 24,192 labelled solutions (2^14.56, as in §8's brute force), against
  log₂E = 13.1.
- **t = 163, k = 9:** complete joins on 3,000 random balanced supports found 1,060 solutions, each checked. That is
  0.353 per support against 0.077 predicted, so the estimate is 2.2 bits low (`logs/tmat/census_2k2.jsonl`).
- **So read 883/k = 21 as n ≈ 110.**

**919/d = 27 (3-cycle)** belongs to the family t = (4/3)d² − 2d + 1, with row sums ((d/3)·(1, 1, 1), ±(d − 1)).
- **Control [M]:** the per-support exact join of `tmat_twist3_census.py`, run over all 64 supports at (19, d = 3),
  reproduces §8's brute-force count of 4,968.
- That is the only exact check of the 3-cycle Gaussian: 2^12.28 exact against 2^11.5 estimated, 0.8 bits low.
- Larger members of the family (271/d = 15, 547/d = 21) and the cousin 127/d = 9 (row sums ((d/3)·(1, 1, 1), d + 1),
  as 19/d = 3) were not counted. The per-support join needs Σ_s C(m′, s)·6^{s−1} ≈ 7^{m′}/6 list entries
  (10^11 at 127/d = 9, 10^14 at 271/d = 15). An importance-sampled attempt at 127/d = 9 was stopped after 13 min
  without finishing 10 supports.
- **So 919/d = 27 carries only that calibration.** A few bits either way do not change §10.5.

### 10.5 Cost [E]

- **883/k = 21 (the k-tree applies).**
  - About 2^−28 solutions per balanced support after the correction.
  - One complete join on lists of 2^20–2^21 costs about 2^41.
  - Total: about 2^69.
- **919/d = 27 (3-cycle).**
  - The P/R split of `tmat_kt.py` needs the rows in fixed pairs; the 3-cycle mixes three rows. So only a 2-list split
    per support is available: the 3-row block (6^s labels) against the fourth row (2^{34−s} signs).
  - At the best support size (s ≈ 11): about 2^−38.7 labelled solutions per support, i.e. 2^−42.3 after the 12
    symmetries the lists remove. Each support costs about 2^26 operations.
  - Total: about 2^68, or 2^67 with the 0.8-bit correction.
- **Comparison.** Both are about 2^40–2^42 beyond T(131), which took one support of about 2^27 pair operations
  (23 s). 669 stays at about 2^79 (§7).

### 10.6 Findable T-matrices miss the open list [E]

Menu of primes 100–420, odd k, 5–60 orbits (`logs/tmat/menu_100_420.jsonl`). Row-sum-feasible slices with n ≤ 81:

| t/k | orbits | n (bits) | v = t·w < 1000 (w odd ≥ 3) |
|---|---|---|---|
| 113/7 | 17 | 31.4 | 339, 565, 791 |
| 101/5 | 21 | 37.6 | 303, 505, 707, 909 |
| 163/9 | 19 | 37.7 | 489, 815 (1,060 T(163) found in §10.4) |
| 181/9 | 21 | 42.5 | 543, 905 |
| 131/5 | 27 | 51.0 | 393, 655, 917 (built, §9) |
| 313/13 | 25 | 55.2 | 939 |
| 197/7 | 29 | 59.3 | 591, 985 |
| 151/5 | 31 | 60.2 | 453, 755 |
| 103/3 | 35 | 62.9 | 309, 515 (2060), 721 (built, §1), 927 |
| 211/7 | 31 | 64.1 | 633 |
| 109/3 | 37 | 67.3 | 327, 545, 763, 981 |
| 181/5 | 37 | 74.4 | as 181/9 |
| 307/9 | 35 | 77.4 | 921 |
| 127/3 | 43 | 80.7 | 381, 635, 889 |

- None of these v is open.
- These T-matrices matter only for open orders above 4000 that have a factor t from this table and a Williamson-type
  cofactor.

**Reproduce:**

```
python3 scripts/tmat_primev.py 571,599,643,659,719,739,751,823,839,859,863,883,907,919,947,971,991
python3 scripts/tmat_primev.py --asym 179,191,263
python3 scripts/tmat_twist.py screen 263
python3 scripts/tmat_twist3_open.py 571,643,739,751,823,859,883,907,919,991,223
python3 scripts/tmat_reflection_check.py
python3 scripts/tmat_twist_d2_cal.py 7,11
python3 scripts/tmat_census_2k2.py 3000
python3 scripts/tmat_twist3_census.py 19 3 all
python3 scripts/tmat_menu.py 0/4 100 420     # parts 0/4 .. 3/4
```
