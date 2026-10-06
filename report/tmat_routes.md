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
