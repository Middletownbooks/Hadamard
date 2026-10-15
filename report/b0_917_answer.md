# H(3668): b = 0 at v = 917, built by Cooper–Wallis from cyclotomic T-matrices of order 131; the index-26 U = ⟨2⟩ slice costed

Date: 2026-10-15. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_917_b0_design.md`.
Scripts: `theory/b917_*.py`, `scripts/verify_h3668.py`. Logs: `logs/b917/`. Bundle: `pub/H3668/`.

**Tags.** [M] exact computation (script named); [P] proved here (proof in the text); [E] Gaussian estimate (estimator
and calibration named); [L] literature recollection, unchecked (arXiv and OEIS are blocked).

---

## 0. Bottom line

1. **(A) H(3668) is built: v = 917 = 7·131, on the open list; b = 0 (Goethals–Seidel, four ±1 circulants of length
   917 with Σ_k PAF_k = 3668·δ) [M].**
   - **Construction (Cooper–Wallis).**
     - T-matrices of order 131: four disjoint {0, ±1} circulants with Σ PAF = 131·δ, invariant under the order-5
       subgroup K = ⟨53⟩ = {1, 53, 58, 61, 89} of F₁₃₁^* (27 orbits).
     - Williamson matrices of order 7: four symmetric ±1 circulants with Σ PAF = 28·δ.
     - Their Cooper–Wallis product on Z₁₃₁ × Z₇ ≅ Z₉₁₇ gives the four sequences (§1.1).
   - **Search.** The T-matrices were found by the repo's support-split k-tree search (`scripts/tmat_kt.py 131 5`), on
     the second support tried, in **23 s** on one core.
   - **Checks.**
     - Σ PAF = 3668·δ on Z₉₁₇ exactly, and on Z₁₃₁ × Z₇ (`theory/b917_cw.py`).
     - **H Hᵀ = 3668·I exactly** for the Goethals–Seidel array, by the independent `scripts/verify_h3668.py`
       (integer dot products, int64 row blocks, no shared code).
     - The same by the standalone `pub/H3668/verify.py`, together with the published normalized matrix.
     - Separately, the Kronecker-form Cooper–Wallis array built straight from the T-matrices
       (`theory/b917_cw_kron.py`, the code that built H(2884)) is also exactly Hadamard.
   - sha256 of the packed rows (GS array): `49eb87f8ca4b8838b1fa176a8579d95cbe4dab07658872ff59dd37f96a9a8d79`.
   - **A second, different H(3668)** came from re-running the search with another seed (`theory/b917_cw.py --search 3`:
     the eighth support, 118 s). Its T-matrices have row sums (−9, 5, −5, 0); its core has row sums (−27, 47, 17, 21).
     It verifies the same way (`scripts/verify_h3668.py`: Hadamard, same slice, sha256 `c94240f8…5dd4`).
   - **Slice.** The core's full multiplier stabiliser in Z₉₁₇^* is exactly K × {±1}, of order 10 [M]. So it lies in
     the index-26 GP slice of Z₇ × Z₁₃₁ with **U = −1 on Z₇**: s = 26, f = 5, 108 orbits per block.
     - That is not the U = ⟨2⟩ slice of the prompt. It is the neighbouring slice of the same sweep
       (`logs/ab/sweep_b0.jsonl`, passes screens A and B), which the earlier ranking had not costed.
     - Inside it the core lies in the Cooper–Wallis product sub-slice: each column of each block is ± one of the four
       Williamson sequences, which leaves 3 bits per Z₁₃₁-orbit.
2. **Why this route works at 917, and only there among the open composite orders [P, E].**
   - T-matrices of order t are a b = 0 GS core of order t with the product condition a_j b_j c_j d_j = +1
     (`report/tmat_routes.md` §2). In the K-invariant ansatz their row sums obey an exact condition: three ≡ 0
     (mod k), one ≡ ±1 (mod k), Σ s² = t (§7 there).
   - **At t = 131, k = 5 the condition is met: (9, 5, 5, 0), 81 + 25 + 25 = 131, 9 ≡ −1 (mod 5) [P].**
     - The estimate (`scripts/tmat_estimate.py 131 5`) is **log₂E = 30.0** labelled solutions, 16.7 per symmetry
       class, against log₂E = 42.1 for the T(103, k = 3) control that the repo found blind in 720 s.
     - k = 13 and k = 65 fail the row-sum condition [P].
   - Census over 40 random supports (§1.4, exact complete joins): 0.70 solutions per support against 0.5 predicted.
     The estimate is low by 0.5 bits, as on the control ladder (§1.5, 0.5–0.9 bits low).
   - The other composite open orders:
     - 955 = 5·191, 537 = 3·179 and 789 = 3·263 are row-sum infeasible for every odd k [P] (§3).
     - 669 = 3·223 survives only at k = 3: 75 orbits, about 157 constraint bits (`report/tmat_routes.md` §7).
   - **917 was the one open order whose T-matrix factor is small.**
3. **The prompt's slice (index 26, U = ⟨2⟩) is costed, not run [P, M, E] (§2).**
   - **Exact dual-energy lemma [P, M]:** N₅ + 2N₇ = 26 (Σ A² = 28) or 25 (Σ A² = 148) over the 104 class columns. It
     costs 35.2 bits, exactly (`theory/b917_levels.py`).
   - With it imposed exactly, the remaining Gaussian is centred (Mahalanobis 0.06, against 6.1 unconditioned), and
     the slice holds **+55.2 bits per y** [E] (the unconditioned estimate was +52.4). The split is C_dual 71.8 and
     C_cross 118.1.
   - The same estimator runs 1.5 bits high against exact core counts at v = 217, s = 6 (§2.5).
   - Given the dual data w, the cross fibre (anti-column signs) is 71–75 bits, against C(cross | w) ≈ 107–108 bits at
     the fibre mean [E].
     - **Dual-first deficit: 32–37 bits at best, about 42 on average**, against 7.4 at 933.
     - With a ≥ 2^35 fibre join per dual solution, that is ≥ 2^67 per core.
   - A plain 2+2 MITM over blocks costs 2^(C − H/2) = 2^66.3 per core.
   - Every level order the prompt lists was costed against these figures (§2.4). The 150 core-hour budget is about
     2^45.6 operations. The realistic orders need 2^62 or more. Even with a hypothetical representation-type solver
     for the linear (C−) level, the best is about 2^53 plus the cost of generating dual solutions.
   - The same-shape control v = 553 (q = 79, s = 26) costs 2^54 by plain MITM [E], so no s = 26 control is runnable
     either.
   - The s = 6, 10, 14 controls range from 2^5.7 to 2^27 (§2.5).
4. **(B)** The method that built H(3668) is:
   - the cyclotomic T-matrix ansatz, screened exactly for row sums (`scripts/tmat_estimate.py`);
   - the support-split k-tree join (`scripts/tmat_kt.py`);
   - Cooper–Wallis with Williamson matrices of the cofactor order (`theory/b917_cw.py`).

   **It was validated end to end on a ladder of the same shape** (k = 5 T-matrices × W(7); §1.5):
   - t = 41 → H(1148), t = 61 → H(1708) and t = 101 → H(2828), each built and verified by `scripts/verify_h3668.py`;
   - exact total counts at t = 41 and 61 agree with the estimator to within 0.9 and 0.5 bits (it undercounts);
   - earlier, T(103, k = 3) → H(2884) (`cw/RESULT_2884.md`).

   **Cost at 917:**
   - about 3 minutes per T(131) by complete exact joins (one independent hit per about 7 random supports, 27 s each);
   - 23 s and 118 s in the two k-tree runs (§1.3, §1.4).
5. **(C) Exact obstructions** (§3):
   - T(131) for k = 13 and k = 65;
   - cyclotomic T-matrix routes for 537, 789 and 955 (every odd k);
   - Williamson-type ⟨2⟩-invariant quadruples of order 7 are never amicable, so Cooper–Wallis cannot reach the U = ⟨2⟩
     slice.
6. **[L]** The literature status of T-matrices of order 131 was not checkable here. The known infinite families (base
   sequences BS(66, 65), Turyn-type TT(44), Golay/Yang lengths) do not obviously give 131, and 3668 is on the user's
   open list, which is consistent with T(131) being new. If T(131) were already known, H(3668) would follow by the
   same Cooper–Wallis step.

---

## 1. The construction

### 1.1 T-matrices and Cooper–Wallis [P]

- **Inputs.**
  - T₁..T₄: {0, ±1} circulants of order t with disjoint supports covering Z_t, and Σ_i PAF(T_i) = t·δ.
  - W₁..W₄ = A, B, C, D: symmetric ±1 circulants of order w with Σ PAF(W_i) = 4w·δ (Williamson).
- **The four sequences.** On G = Z_t × Z_w (products as functions):

      x₁ =  T₁A + T₂B + T₃C + T₄D,     x₂ = −T₁B + T₂A + T₃D − T₄C,
      x₃ = −T₁C − T₂D + T₃A + T₄B,     x₄ = −T₁D + T₂C − T₃B + T₄A.

  - Each x_k is ±1: exactly one T_i is nonzero at each point of Z_t.
  - **Diagonal terms.** Each T_i meets every W_j exactly once across the four x_k, so they add up to
    (Σ_i PAF(T_i)) ⊗ (Σ_j PAF(W_j)) = t·δ ⊗ 4w·δ = 4tw·δ.
  - **Cross terms.** The terms CC(T_i, T_j) ⊗ CC(W_a, W_b), i ≠ j, cancel in pairs under the quaternion sign pattern.
    The reason is that CC(W_a, W_b) = CC(W_b, W_a) is symmetric in the lag when the W's are symmetric (Williamson
    matrices are symmetric circulants).
- **From G to a circulant core.** gcd(t, w) = 1, so G ≅ Z_tw by the CRT, the x_k are circulants of length tw, and the
  Goethals–Seidel array is an H(4tw).
- **At 917.** t = 131 and w = 7. The order-7 Williamson quadruple is the first found by
  `theory/b917_cw.williamson(7)`: A = B = (+++−−++), C = (++−++−+), D = (+−+−−+−).

### 1.2 T(131) in the k = 5 cyclotomic slice [P, E]

- **The slice.**
  - K = ⟨2^26⟩ = ⟨53⟩ = {1, 53, 58, 61, 89} ≤ F₁₃₁^*, of order 5. This is the same M as in the prompt's slice (2 is a
    primitive root mod 131).
  - The T_i are constant on {0} and on the 26 cosets of K: 27 orbits, and a label in {T₁..T₄} × {±} per orbit, so
    8²⁷ = 2⁸¹ labellings.
  - Equations: Σ_i PAF(T_i) = 0 on the 13 classes of ⟨K, −1⟩ on the nonzero lags (−1 ∉ K since 5 is odd).
- **Row-sum condition [P].**
  - Each s_i = T_i(0) + 5·(integer), and exactly one T_i is nonzero at 0. So three s_i ≡ 0 and one ≡ ±1 (mod 5).
  - Σ s_i² = 131 then forces {|s_i|} = {9, 5, 5, 0}:
    - 131 − s₄² must be 25·(a² + b² + c²), so s₄² ≡ 131 ≡ 6 (mod 25), i.e. s₄ ≡ ±9 (mod 25);
    - s₄ = ±9 leaves a² + b² + c² = 2.
  - The core's T has row sums (0, 5, −9, 5).
  - For k = 13: s₄ ≡ ±1 (mod 13) and s₄² ≤ 131 give s₄ ∈ {±1, ±12}. Then 131 − s₄² ∈ {130, −13}, and neither is
    169·(sum of squares). k = 65 likewise. Both are infeasible.
- **Estimate [E].** `scripts/tmat_estimate.py 131 5` (200,000 samples, lattice-aware Gaussian, row-sum filter):
  - log₂E = 30.0 labelled quadruples, 16.7 per symmetry class (4!·2⁴·130/5);
  - the same estimator gives 42.1 for T(103, k = 3);
  - calibration: within 1–1.5 bits (underestimating) against brute force at t = 19 (`report/tmat_routes.md` §8).

### 1.3 Search [M]

`scripts/tmat_kt.py 131 5 --nS 40 --seed 1 --buckets 200` (log `logs/b917/tmat_kt_131_k5_s1.log`):
- Fix S = supp T₁ ∪ supp T₂, a random set of 13 or 14 orbits.
- On S write T₁, T₂ = (P ± R)/2, and likewise P′, R′ on Sᶜ.
- That gives four independent ±1 lists, of 2¹² to 2¹³ entries each, with PAF(P) + PAF(R) + PAF(P′) + PAF(R′) = 2t·δ.
  They are joined by `scripts/kt_generic.kt_join4`.

The first support gave nothing (13.9 s). The second gave a solution after 2.6 s, re-verified inside the search. The T
matrices are saved in `logs/b917/T131_k5.npy` and `.json`, and were re-checked by a separate script (supports partition
Z₁₃₁, Σ PAF = 131·δ by explicit integer sums, K-invariance).

### 1.4 Census: exact counts per support against the estimate [M vs E]

`theory/b917_tmat_census.py --t 131 --k 5` (log `logs/b917/tmat_census_131_5.jsonl`):
- **Join.** A complete exact 2+2 join per random support: all (P, R) on S against all (P′, R′) on Sᶜ, 64-bit linear
  hash, every match re-verified.
- **What it counts.** The side lists fix the first orbit of S and of Sᶜ to sign +1. Each support therefore counts the
  T-matrices with supp T₁ ∪ supp T₂ = S, T₁ = +1 on the first orbit of S and T₃ = +1 on the first orbit of Sᶜ.
- **Prediction** under the estimator's uniform-label model: 2^30.0 · 2^−27 · (1/4)² = 2^−1 ≈ 0.5 per support.

**Result** (40 random supports, sizes 13 and 14, 26.5 s per complete join on one core):

- **28 solutions on 6 supports**, all distinct and all re-verified; per hit support [4, 4, 4, 4, 4, 8].
- **0.7 per support against 0.5 predicted: the estimate is low by 0.49 bits.**
- Solutions come in groups of 4 (sometimes 8). These are the swaps (P, R) ↔ (R, P) and (P′, R′) ↔ (R′, P′), i.e.
  T₂ ↦ −T₂ and T₄ ↦ −T₄ with the same supports.
- Independent hits per support: 6/40 ≈ 0.15.
- **Consequence for the search:** about one T(131) per 7 random supports, at about 27 s per complete support join on one
  core, i.e. **about 3 minutes per T(131)** [M]. The probabilistic k-tree join of `scripts/tmat_kt.py` is faster per
  support (3–14 s) at a small loss of completeness.

### 1.5 Same-shape control ladder [P, M, E]

`theory/b917_tmat_ladder.py --t 11,31,41,61,71,101 --all 41,61 --build` (log `logs/b917/tmat_ladder.jsonl`).

The same ansatz throughout:
- K = the order-5 subgroup of F_t^* (t ≡ 1 mod 10);
- the exact row-sum condition, the estimator, and an exact count over all supports where it is small;
- the k-tree search;
- Cooper–Wallis with the same W(7);
- the independent verifier with the slice check (K × {±1}).

| t | orbits | row-sum patterns (\|s\|) | estimate log₂E | exact labelled count (2 ≤ \|S\| ≤ m − 2) | first solution | built and verified |
|---|---|---|---|---|---|---|
| 11 | 3 | none | – | – | – | – (empty [P]) |
| 31 | 7 | none | −∞ | – | – | – (empty [P]) |
| 41 | 9 | (4, 5, 0, 0) | 13.8 | **2^14.73** (27,136) | 1st support, 3 s | **H(1148)** |
| 61 | 13 | (6, 5, 0, 0) | 18.2 | **2^18.68** (418,560) | 3rd support, 15 s | **H(1708)** |
| 71 | 15 | none | −∞ | – | – | – (empty [P]) |
| 101 | 21 | (1, 10, 0, 0) | 25.4 | – | 9th support, 72 s | **H(2828)** |
| 131 | 27 | (9, 5, 5, 0) | 30.0 | census §1.4 | 2nd support, 23 s (seed 1); 8th, 118 s (seed 3) | **H(3668)** |

- **How the exact counts are taken.** The complete joins count each labelled solution once per 16, by the symmetry
  T₁ ↔ T₂ and the signs on each side [P]. Supports with |S| ∈ {0, 1, m − 1, m} are excluded, so the exact column is a
  slight undercount, and still above the estimate.
- **The estimator undercounts by 0.5–0.9 bits** here, consistent with its 1–1.5-bit undercount at t = 19
  (`report/tmat_routes.md` §8).
- **The search scales smoothly:** 3 s at 9 orbits, 72 s at 21 and 23–118 s at 27. The 2^30 labelled solutions predicted
  at t = 131 were found at the predicted rate (§1.4).
- **The ladder orders are not open.** They serve only as end-to-end controls of the pipeline: search → Cooper–Wallis →
  circulant core → GS array → independent verification.

### 1.6 Verification and bundle [M]

| check | tool | result |
|---|---|---|
| T-matrices: supports partition Z₁₃₁, Σ PAF = 131·δ, K-invariant | inline integer check (§1.3), `cw_build.check_T` | yes |
| Σ_k PAF = 3668·δ on Z₁₃₁ × Z₇ and on Z₉₁₇ | `theory/b917_cw.py` (exact FFT, rounded) | yes |
| Σ_k PAF by integer dot products; GS array; H Hᵀ = 3668·I in int64 blocks | `scripts/verify_h3668.py` (no shared code) | yes; sha256 `49eb87f8…a8d79` |
| invariance under K (mod 131) and −1 (mod 7) | `scripts/verify_h3668.py --slice 131 26 --umult 6` | yes (multipliers 1, 708, 582, 323, 351, 132) |
| full multiplier stabiliser of the core | inline (§0) | K × {±1}, order 10 |
| the published normalized matrix is Hadamard and is the normalization of the built H | `pub/H3668/verify.py` (numpy only) | yes |
| Kronecker-form Cooper–Wallis array from T and W directly | `theory/b917_cw_kron.py` (`scripts/cw_build.build`) | Hadamard; sha256 `2523fe64…fa757` |

- **Bundle** (written by `scripts/make_bundle_b0.py`; b = 0, so no border):
  - `pub/H3668/sequences_Z917.txt`, whose header states the construction and the slice;
  - `pub/H3668/H3668_normalized.txt.gz`;
  - `pub/H3668/verify.py`.
- **Packed matrix:** `data/H3668_packed.npz`.
- **Row sums** (A, B, C, D) = (−17, 9, −7, 57); 17² + 9² + 7² + 57² = 3668.

---

## 2. The index-26 slice with U = ⟨2⟩ (the prompt's target) [P, M, E]

### 2.1 Structure [P, checked M]

`theory/b917_slice.py`; its self-test checks every identity below on random blocks at q = 131, 79 and 31 (errors
≤ 1.5·10⁻¹¹).
- **Character identities.** Stage 1, dual and cross as in the prompt; the dual is paired, D_{a+13} = conj D_a.
- **Half-degree identity.** |D_a|² = (A + Σ_j s_j Re η_{a+j})² + (Σ_j d_j Im η_{a+j})², with s_j = w_j + w_{j+13} and
  d_j = w_j − w_{j+13} (j < 13). Checked to 1.1·10⁻¹¹.
- **Column alphabet.** u = 2u′ with u′ ∈ {0, ±1, ±η, ±η̄} (η = (−1 + i√7)/2), and w² + 24|u′|² = 49.
  - The only fibre of u is the sign of a |w| = 7 column.
  - The only fibre of w is the sign e of an anti column.
- **Coordinates.** With ρ = Re u = x₀ − p and m = Im u/√7 (p = e on sym columns, m = e on anti columns): w = ρ + 7p, so
  w ≡ ρ (mod 7).
- **Cross split.** The cross equations split into
  - (C+) Σ_k [PAF(ρ_k) + 7·PAF(m_k)] = 3668·δ, the symmetric part;
  - (C−) Σ_k [CC(ρ_k, m_k)(z) − CC(ρ_k, m_k)(−z)] = 0, the antisymmetric part.

  Given the dual, (C−) is **linear** in the anti signs and (C+) is a complementary-sequence condition on them.

### 2.2 Stage 1 and the exact energy count [M]

- **Stage 1.** `python3 theory/ab_b0_stage1.py '<cfg>'` reproduces 38,380 y's (up to block order and sign). The
  stage-2 Gaussian of the log (`ab_gauss.stage2_sampled`) is H₂ 277.6, C₂ 225.2, +52.4, Mahalanobis 6.1.
- **Energy lemma [P].** Given y, block k's dual energy is
  A_k² + 5 Σ_j w_j² = A_k² + 5(26 + 24 n₅ + 48 n₇), where n₅ and n₇ are the numbers of columns with |w| = 5 and 7.
  Its sum over the blocks must be 3668, so N₅ + 2N₇ = (3148 − Σ A²)/120 = 26 or 25.
- **Exact count.** `theory/b917_levels.py` counts the configurations per (n₅, n₇) by multinomials. The energy condition
  costs **35.17 bits exactly** at y #0, #1 and #2 (H₂ 278.22 → H_E 243.05).
- **Gaussian under the energy condition.** Uniform sampling of x given (y, energy), by the exact weights, then the
  lattice Gaussian (`theory/q523_gauss.gauss_bits`):

  | y | H_E | C_all | C_dual (rank 12) | C_cross (rank 25) | Mahalanobis | log₂E per y |
  |---|---|---|---|---|---|---|
  | #0 | 243.05 | 187.83 | 71.78 | 118.10 | 0.06 | **+55.22** |
  | #1 | 243.05 | 187.95 | 71.75 | 118.26 | 0.05 | +55.10 |
  | #2 | 243.05 | 187.93 | 71.78 | 118.18 | 0.06 | +55.12 |

  The unconditioned estimate (+52.4) was dominated by the energy tail (Mahalanobis 6.1). With the energy exact, the
  remaining Gaussian is centred, the regime where b = 0 Gaussians were calibrated to about 1.5 bits at 933.

### 2.3 Conditional budgets [M fibres, E Gaussians]

`theory/b917_cond.py`, y #0, six configurations drawn given the energy:

| anti columns per block | F_e = log₂ fibre of w (exact) | C(cross \| w) at the fibre mean | C(C+ \| w) | C(C− \| w) | deficit |
|---|---|---|---|---|---|
| 20, 22, 20, 24 | 75.1 | 107.9 | 52.4 | 67.6 | 32.8 |
| 18, 22, 22, 24 | 75.1 | 106.5 | 52.3 | 66.2 | 31.5 |
| 22, 18, 22, 22 | 73.2 | 108.0 | 51.8 | 68.3 | 34.8 |
| 22, 16, 22, 22 | 71.3 | 108.1 | 51.5 | 68.7 | 36.8 |
| 20, 18, 24, 22 | 73.2 | 107.5 | 52.0 | 67.7 | 34.3 |
| 22, 20, 18, 22 | 71.3 | 107.9 | 51.3 | 68.7 | 36.7 |

- **Averaged over dual solutions.** The fibre mean of the cross does not sit at the target. So the cross constraint
  per dual solution averages about 116 bits: C_cross = 118.1, less the 2.1 bits that the dual and cross share
  (C_dual + C_cross − C_all). The average deficit is about 42 bits.
- **The other direction.** Given u, the dual's fibre (the signs of the |w| = 7 columns) has 0–3.6 bits against the
  72-bit dual, so cross-first is hopeless.

### 2.4 Costs of the level orders [E]

H = H_E ≈ 243, C ≈ 188 per y; F_e ≈ 74; C(cross | w) ≈ 108 (best) to 116 (mean).

| order | cost per core | why |
|---|---|---|
| 933 order: dual first, then the exact 2+2 cross join over fibres | ≥ 2^(deficit + F/2) ≈ 2^67–2^79 | 2^32–2^42 dual solutions per core, 2^35–2^37 per cross join |
| same, with (C−) solved first by a representation-type linear solver | ≥ 2^(32 + 21) ≈ 2^53, plus dual generation | the yield per dual solution does not change |
| plain 2+2 MITM over blocks, all equations | 2^(C − H/2) = 2^66.3 | generic floor |
| 2+2 with the pair dual-sum as middle value | 2^(C(cross \| dual) − log₂\|L(V)\|) ≈ 2^62–2^71 | the dual is rigid: a pair's dual sum carries about 76 bits, nearly its whole w-entropy |
| Wagner 4-list | ≥ 2^62.7 | needs lists of 2^(C/3) = 2^62.7 > 2^60.75 per block |
| subfield split of the dual (8 lists) | no gain | it enumerates dual solutions (not the bottleneck). Fixing the even/odd shapes costs about 13 bits per block, and the cross deficit is unchanged |
| type pattern first | ≥ 2^66 | the types carry about 75–83 bits. Given types the dual is dense on the x₀ row (72 of about 85 bits), but the cross deficit (32–42) remains |
| cross first | hopeless | the dual's fibre given u is 0–3.6 bits against 72 |
| sub-slices that kill (C−) (per-block symmetry, reflected pairs, pure 3-anti + 1-S5 blocks, κ-constant anti signs) | dead or marginal | each costs 55–80 bits against +55 + 15. Two reflected pairs are impossible: 2(r₀² + r₂²) = 3668 needs 1834 = r₀² + r₂², and 1834 = 2·7·131 is not a sum of two squares [P] |

The decisive figure is the cross level: 25 equations of about 4.3–4.7 bits each (C+: 52, C−: 67), against at most 26
private bits per block (one anti sign per column). The energy lemma forces about 20% of the columns to be sym, which
caps F at about 80. So the cross can never be fibre-balanced in this slice [E]. The growth is linear in s: the deficit
is about 1.3·s + 8 bits (§2.5).

### 2.5 Same-shape controls [M stage 1, E budgets]

`theory/b917_levels.py --q Q --s S` (y #0 of each; stage 1 exact):

| v | q | s | f | stage-1 y's | H_E | C_dual | C_cross | log₂E per y | plain MITM per core |
|---|---|---|---|---|---|---|---|---|---|
| 217 | 31 | 6 | 5 | 113 | 41.6 | 9.0 | 17.8 | +15.1 | 2^5.7 |
| 301 | 43 | 6 | 7 | 24 | 41.0 | 9.1 | 18.9 | +13.2 | 2^7.3 |
| 217 | 31 | 10 | 3 | 42,845 | 77.9 | 18.9 | 32.5 | +27.2 | 2^11.7 |
| 497 | 71 | 10 | 7 | 88 | 80.0 | 21.4 | 37.8 | +21.5 | 2^18.5 |
| 301 | 43 | 14 | 3 | 144,198 | 117.4 | 30.4 | 50.4 | +37.5 | 2^21.2 |
| 497 | 71 | 14 | 5 | 3,499 | 118.4 | 32.5 | 55.0 | +31.9 | 2^27.4 |
| 553 | 79 | 26 | 3 | 1,721,346 | 240.6 | 67.4 | 108.9 | +66.3 | 2^54.0 |
| **917** | 131 | 26 | 5 | 38,380 | 243.1 | 71.8 | 118.1 | +55.2 | **2^66.3** |

- **Calibration against exact counts [M vs E]** (`theory/b917_exact_small.py`). At v = 217 (q = 31, s = 6, f = 5,
  the same f as 917), a complete exact 2+2 join over all block configurations counts every core:
  - per y: 12,096 or 12,012 cores (2^13.56), against 2^15.06 from the energy-exact lattice Gaussian;
  - **the Gaussian is 1.5 bits high**, uniformly over 6 y's (`logs/b917/exact_small.jsonl`).

  This is the same ±1.5 bits as at 933, so the 917 figures (+55.2 per y, 2^66 per core) hold to within a couple of
  bits.
- The per-equation costs are stable across the family: dual 4.5–6, cross 3.6–4.7.
- The s ≤ 10 rows are where the 933 pipeline's deficit is small, about 5–10 bits [E, from F ≈ 3.2s − 8 and
  C(cross | w) ≈ 4.2(s − 1)].
- At s = 26 neither q is in reach. The prompt's "same-shape control first" cannot be run at s = 26 in 150 core-hours
  either.

---

## 3. Exact obstructions [P, M]

1. **T(131) with k = 13 or k = 65 is empty [P]** (row sums, §1.2). Hence k = 5 is the only cyclotomic T(131) ansatz.
2. **Cyclotomic T-matrix routes for the other composite open orders [P]** (the row-sum condition of
   `report/tmat_routes.md` §7):
   - 955 = 5·191: every odd k | 190 (5, 19, 95) is infeasible (also §7 there).
   - 537 = 3·179 and 789 = 3·263: the only odd k are 89 and 131, which leave 3 orbits; three row sums ≡ 0 (mod k)
     force them to 0, and t is not a square.
   - 669 = 3·223: only k = 3 survives (75 orbits, about 157 bits).
   - The prime open orders have no Cooper–Wallis split.
3. **Cooper–Wallis cannot reach the U = ⟨2⟩ slice [P].**
   - A ⟨2⟩-invariant ±1 sequence on Z₇ is (w₀, w₁ on R₁, w₂ on R₂).
   - A symmetric one is constant on Z₇^* and has off-peak PAF 2w₀w₁ + 5 > 0, so no Williamson quadruple is
     ⟨2⟩-invariant.
   - The ⟨2⟩-invariant GS quadruples of order 7 (one sym, three anti, PAF 3 − 1 − 1 − 1 = 0) are not amicable. For an
     anti f = (a, e, −e) and a sym g = (b, −b, −b): CC(f, g)(1) − CC(g, f)(1) = −4be ≠ 0.
   - So the U = ⟨2⟩ slice needs a non-product search, costed in §2.
4. **Two reflected block pairs are impossible in the U = ⟨2⟩ slice [P]:** 1834 = 2·7·131 is not a sum of two squares
   (§2.4).

---

## 4. Notes for the other 917 work and what remains

- **The c = 6 route for 917** (the other session, q = 19, s = 6, D_def 70 [E]) is no longer needed for 917.
- **The U = ⟨2⟩ index-26 slice stays unsolved.** By §2.4 it needs about 2^60 or more per core with known orders. A
  core there would be a second, structurally different H(3668), and no longer a priority.
- **The T-matrix route is worth recording as a screen for any new composite order 7t, 5t, 3t:** an exact row-sum
  check, then the estimate, then the k-tree (seconds to minutes when m ≲ 30 orbits).
- **[L] Literature.** Whether T-matrices of order 131 were known should be checked against Đoković's and Kotsireas's
  T-matrix/T-sequence tables and against recent Hadamard order lists before any announcement.

---

## 5. Reproduce

```
python3 scripts/tmat_estimate.py 131 5,13 200000                      # row-sum feasibility, log2E = 30.0 (k = 5)
python3 scripts/tmat_kt.py 131 5 --nS 40 --seed 1 --buckets 200       # T(131): second support, 23 s
python3 theory/b917_cw.py                                             # Cooper-Wallis -> logs/b917/seq_H3668_Z917.txt
python3 theory/b917_cw.py --search 3                                  # same, re-running the T-matrix search (seed 3)
python3 scripts/verify_h3668.py logs/b917/seq_H3668_Z917.txt data/H3668_packed.npz --slice 131 26 --umult 6
python3 theory/b917_cw_kron.py                                        # second construction (Kronecker CW array)
cd pub/H3668 && python3 verify.py                                     # standalone bundle check
python3 theory/b917_tmat_census.py --t 131 --k 5 --nS 40 --seed 11    # census per support (exact joins)
python3 theory/b917_tmat_ladder.py --t 11,31,41,61,71,101 --all 41,61 --build   # control ladder: H(1148), H(1708), H(2828)
python3 theory/b917_slice.py --q 131 --s 26                           # slice identities (self-test)
python3 theory/b917_levels.py --q 131 --s 26 --ys 0,1,2 --N 12000     # exact energy count + conditional Gaussians
python3 theory/b917_cond.py --q 131 --s 26 --ys 0 --n 6 --N 4000      # cross | w, dual | u budgets
python3 theory/b917_exact_small.py --q 31 --s 6 --ys 0,1,2,3,4,5        # exact core counts at v = 217 vs the Gaussian
```
