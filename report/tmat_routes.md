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
