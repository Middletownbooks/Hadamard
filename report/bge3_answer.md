# b ≥ 3 bordered Goethals–Seidel cores for the open orders 2000–3000: answers to `prompt_bge3_opus.md`

Date: 2026-10-06. Tags: **[T]** derivation (proof sketch given), **[M]** measured or verified in exact integer arithmetic,
**[E]** estimate. New scripts and logs are listed in §8. Notation follows the prompt: G = Q × H, |Q| = c, |H| = h,
b = c − 1, m = ch, v = m + b; targets Σ_k |X̂_k(ψ, ω)|² = 4v, except 4b at (ψ ≠ 1, ω = 0).

---

## 0. Short answers

1. **Algorithm to about 2^50 (Q1): none exists within list-based methods; the barrier is quantified below.**
   - [T] Any method that treats the four blocks as list elements needs time T ≳ 2^max(n/3, (n − log₂E)/2), with
     T·M² ≈ 2^n once memory M is below the full list size. This covers k-tree, MITM, Schroeppel–Shamir and
     compression-and-lift (§2.1).
   - [M] Calibration against exact solution counts (§2.2):
     - The refined model (`bge3_cost2.py`) is within 1.5 bits whenever solutions exist.
     - The plain Gaussian model is 1–3 bits pessimistic per profile.
     - Neither model is pessimistic enough to change any conclusion.
   - [M/E] On that model, the only open b ≥ 3 configuration still alive with a cost estimate is **571 / c = 22 / d = 10**.
     It needs **≳ 2^80**: n_lift = 168 bits after fixing the compression, log₂E ≈ 23.5. Every cheaper configuration in
     the prompt's table is now dead (§4.2).
   - [E] The published 479 and 347 cores sit at **≈ 2^71 and ≈ 2^78** on the same model. Whoever found them used
     about that many operations or a method that is not list-based. I could not identify which (§5).
   - [T] The Q-character decomposition does not yield more independent lists. Decoupling the per-ψ problems requires
     fixing the column "types" (§2.4):
     - That costs about 2 bits per column orbit, about 2r bits in total.
     - Every configuration has log₂E < 2r: decoupling would need log₂E of about 90–130, and the largest available is 86.
     - The T-matrix support split worked because it cost only 1 bit per orbit, and T(103) had log₂E = 42 > 35.
2. **Algebraic construction (Q2): none known, and the natural places to look are empty.**
   - [M] Exact decision of **182 high-symmetry slices** (r ≤ 22–26 orbits per block) for 571/c22, 599/c20 (both Q),
     647/c12, 659/c6, 739/c10, 719/c8 and 719/c10: **all empty**, with a positive control (§3.3). Any cyclotomic or
     Paley-type construction with large multiplier symmetry would live in these slices.
   - [T] Three natural ansätze are impossible (§3.2):
     - rows constant on QR/QNR (Paley): impossible for h > 7, by counting;
     - rank-1 Legendre cores: always overshoot the target by exactly 4 at every ω ≠ 0;
     - the c = 4 Legendre grid: requires a perfect binary sequence of odd length, which does not exist.
   - I know of no theorem (Miyamoto, Xia–Xia–Seberry, Ito, Paley/Szekeres) that covers any of these (v, c, h).
3. **Screens (Q3):** several configurations that the prompt lists as undecided or queued are now decided.
   - [M] An exact compression DP (`bge3_cost2.py`) kills:
     - 647/c12/d13;
     - 659/c6/d9, the best-looking 659 option (log₂E 86 under the plain model);
     - 599/c20/Z2×Z10, four groups;
     - 599/c20/Z20, group 0.
   - [T/M] **The border is a real extra condition for c ≥ 6.** The corner is forced: C = −K·Cbᵀ·Kcᵀ/(4c) must be ±1.
     - c = 4: every compression has a border (2048/2048).
     - c = 6: borders are known only for 347's exact compression placement (at most 1 of 150 symmetry-equivalent
       placements works with 347's frame).
     - **c ≥ 10: no border is known.** That covers every large-c open configuration (571, 599, 647, 739, 719/c20, 643, 573).
   - [T] A row-sum granularity lemma kills multiplier symmetries with large orbits on fixed rows.
   - [M] Lattice (2-adic) obstructions are real: one calibration slice predicted at 2^9 solutions has exactly 0, because
     the target is unreachable mod 32.
4. **Forensics (Q4)** (§5):
   - All three published cores use h = h₁h₂ with coprime h₁, h₂ > 2. This fits a two-sided compression method, but the
     cost formula in §2.6 shows such a method does not beat the barrier in general.
   - 479 and 359 use monomial compressions, in a placement that makes the border easy to build.
   - 479's twist reduces the compression to one parameter per block.
   - There is no hidden product or low-rank structure: an apparent 15σ signal vanishes against a symmetry-aware null.
5. **What to implement first** (§7):
   1. Use `bge3_cost2.py` plus the repo's joint-residue filter as the standard screen.
   2. Settle border existence for c ≥ 10 before any large-c search.
   3. Do not put Kaggle or local time into b ≥ 3 for open v: 523/b0/d9 remains the best route.
   4. If the goal is to learn the published method, the right positive control is to reproduce a solution in 479's own
      slice. Any method claiming 2^50 on open configurations should do that in ≲ 2^40.

---

## 1. Facts used throughout

- [T] **Character form.** Σ_k |X̂_k(ψ, ω)|² = 4v − 4m·[ψ ≠ 1, ω = 0].
  - In the group ring this reads Σ_k X_k X_k^(−1) = 4v − 4c·H + 4·G.
  - Equivalently, the four ±1 blocks form a divisible difference family relative to H, with λ₁ = λ₂ − c.
- [T] **Q-compression (ω = 0).** Write x̄_k = μ_k·1 + d_k, with x̄_k(q) ≡ h (mod 2). Then:
  - Σ_k |d̂_k(ψ)|² = 4b for ψ ≠ 1;
  - Σ_k ‖d_k‖² = 4b²/c;
  - Σ_{k,q} x̄_{k,q}² = 4(h + b).
- [T] **For each ω ≠ 0** the c-vectors Z_{k,ω}(q) = r̂_{k,q}(ω) are complementary over Q with power 4v. Equivalently, the
  4c × 4c GS array of the row-Fourier coefficients is √(4v) times a unitary matrix.
- [T] **Parseval per ψ.** Σ_k ‖Y_{k,ψ}‖² = 4v − 4c·[ψ ≠ 1], and Σ_k PAF_H(Y_{k,ψ}) = 4vδ − 4c·[ψ ≠ 1].
  - For c = 4 (Z2²) this reads 4N₁ + n_★ = v and 4N_ψ + n_★ = v − 4 for ψ ≠ 1. Here N_ψ counts the columns concentrated on
    ψ and n_★ counts the 3–1 columns.
  - Hence n_★ ≡ v (mod 4) and n_★ ≥ v mod 4.

---

## 2. Q1: algorithms

### 2.1 The barrier for list-based methods [T]

Take four lists of sizes N_k, with n constraint bits and E = ΠN_k·2^−n expected solutions.
- The bucketed k-tree finds a solution in time T ≈ 2^n/M² using lists of size M ≤ N.
- When N³ ≥ 2^n, the k-tree with sublists of size 2^(n/3) takes T = 2^(n/3).
- Hence T_min ≈ 2^max(n/3, (n − log₂E)/2), reached only with memory ≈ 2^((n + log₂E)/4).
- With realistic memory (M ≤ 2^30 entries) the cost is T ≈ 2^(n − 60).

Representations (§2.5) change memory, not T_min.

More independent lists help only in the k-tree regime E ≥ 2^(n/3): 2^j lists give 2^(n/(j+1)) when lists of size
2^(n/(j+1)) exist. Every live b ≥ 3 configuration has E < 2^(n/3) (571: log₂E ≈ 23.5, while n/3 ≈ 56–64). The one
exception was 659/c6/d9, which is dead. Below that threshold, extra lists do not lower T below 2^((n − log₂E)/2).

The rule reproduces the repo's own successes:

| instance | n | log₂E | T_min | observed |
|---|---|---|---|---|
| T(103) support split | 63 | 7 per S | 2^29 (lists 2^17) | 12 min |
| 487, b = 0, flat | n′ ≈ 83 | ≈ 29 | ≈ 2^28 | found |
| 523, b = 0, flat | n′ ≈ 91 | small | ≈ 2^31 per hit | ≈ 0.3 hits so far |

Applied to b ≥ 3 (refined numbers from §2.3 where available):

| configuration | n (compression fixed) | per-quadruple log₂E | T_min = 2^ | T at M = 2^30 | status |
|---|---|---|---|---|---|
| **479 published** (c4, d12, r46) | 145.0 | 2.5 | **71** | 2^85 | solved by others |
| **347 published** (c6, d6, r70) | 187.2 | 32.1 | **78** | 2^127 | solved by others |
| 571 c22 d10 r62 | 168.0 | 7.8 (total 23.5) | **80** | 2^108 | alive; border unknown (c = 22) |
| 719 c8 d22 | 135 (plain) | 18 | 58.5 | 2^75 | **dead** (no admissible compression; confirmed with `bge3_cost2`) |
| 659 c6 d9 | 204 (plain) | 86 | 68 | 2^144 | **dead** (§4.2) |
| 647 c12 d13 | 179 (plain) | 37 | 71 | 2^119 | **dead** (§4.2) |
| 599 c20 (Z2×Z10, all groups; Z20, both groups) | 184–200 | ≤ 12 | ≥ 87 | – | **dead** (§4.2 and the norm-budget screen) |
| 643 c14 d12 | 206.6 (refined) | 2.0 total (plain 29) | ≈ 102 | 2^147 | alive; border unknown (c = 14) |
| 647 c18 d12 | 233 (plain) | 36 | ≥ 98 | – | DP overflow |
| 523 c4 d12 | 219.2 (refined, group 1) | **−32** (plain 3.5); group 0 flagged infeasible | – | – | effectively dead |

**Reading.** b ≥ 3 is intrinsically harder per orbit than b = 0:
- n/L ≈ 5.5–6 bits per lag class, against ≈ 4.9 for b = 0;
- r/L ≈ 1.6–2.

So the exponent n − 2·log₂N ≈ L(β − 2r/L) is about 2L for b ≥ 3, against about L for b = 0.

### 2.2 Is the model right? Exact calibration [M]

`scripts/bge3_exact_small.py` and `scripts/bge3_exact_count.py` enumerate all 2^r S-invariant blocks of a slice and
apply the exact per-block PSD bound |X̂(χ)|² ≤ T̂(χ). They then count ordered 4-tuples by pair-hashing.

c = 4, Q = Z2²; log₂ of solutions:

| v (h) | d | r | L | plain model, per profile | exact, per profile | refined model, all profiles | exact, all profiles |
|---|---|---|---|---|---|---|---|
| 23 (5) | 2 | 11 | 8 | 7.4 | **10.0** | 16.1 | **17.58** |
| 31 (7) | 3 | 10 | 5 | 14.9 | **16.56** | 24.0 | **24.15** |
| 55 (13) | 3 | 18 | 9 | 30.7 | **32.02** | 39.5 | **40.5** |
| 55 (13) | 3 | 20 | 11 | **22.3** | **0** | no admissible compression (0) | 0 |
| 63 (15) | 4 | 19 | 13 | 4.3 | **0** | 9.1 (sampler fixed) | 0: target unreachable mod 32 |
| 63 (15) | 4 | 21/22 | 17 | −1.4 / −0.9 | 0 | – | 0 |
| 6 further slices | | 16–24 | | −30 … −5 | 0 | – | 0 |

- **When solutions exist,** the refined model is within 1.5 bits and the plain one is 1.3–2.6 bits pessimistic.
- **When the model predicts solutions but there are none,** the cause is an exact obstruction:
  - a compression with no admissible quadruple, which `bge3_cost2.py` detects;
  - a 2-adic lattice obstruction, detected by the residue probe or the repo's joint-residue filter.

So the model errs on the optimistic side, never by tens of bits on the pessimistic side.
**The 2^71–2^78 figures for the published cores are not a modelling artifact at these sizes.**

### 2.3 Refined cost model [M/E] (`scripts/bge3_cost2.py`)

How it works:
- The compression is linear in the orbit signs (x̄ = M·y), so its distribution is computed exactly by a pruned DP over
  orbits.
- All admissible compression quadruples are enumerated, together with their block counts.
- The Gaussian model is applied only to what is left, using samples drawn from the exact compression fibre.

**Sampler fix.** The fibre moves must include changes of the "+ count" along small integer relations
Σ_g D_g·M_g = 0 between orbit-signature groups. Without them some directions look deterministic.
- This happened at first for 523 and 643: "null mismatch", now corrected.
- Count vectors are not unique: for 523, two size-1 orbits ↔ one size-2 orbit in the same coset.

Results:
- 479 real: compression costs 10 bits, n_lift 145.0, log₂E 11.0 (plain 7.0).
- 347 real: compression costs 37 bits, n_lift 187.2, log₂E 41.3 (plain 55.5).
- 571 c22 d10: 522 compressions per block on 6 class-orbits; 60,384 admissible quadruples. **The compression costs
  45.8 bits; n_lift = 168.0; log₂E = 23.5** (plain: n = 192, log₂E = 30).
- 523 c4 Z4 d12, rerun with the fixed sampler:
  - group 1: compression 54 bits, n_lift 219, **log₂E = −32** (plain model +3.5);
  - group 0: 6 lag-class combinations stay deterministic given the compression and miss the target, so it is flagged
    infeasible. This agrees with the norm-budget screen (`bge3_normscreen.py`), which fails group 0 at t = 2 and t = 10.
- 643 c14 d12: compression 69 bits, n_lift 207, **log₂E = 2.0** (plain 29).
- Both get worse, not better.

### 2.4 The per-character grid does not decouple [T]

The structure is a 4 × c grid of sequences Y_{k,ψ} on H.
- **Columns** (fixed ψ): one 4-list condition each, Σ_k |Ŷ_{k,ψ}(ω)|² = target_ψ(ω).
- **Rows** (fixed block k): coupled pointwise, since (Y_{k,ψ}(z))_ψ is the Q-Fourier transform of a ±1 column.

Compare T-matrices: there the coupling (disjoint supports) and the condition (one sum over the same four sequences) live
on the same index. Here they live on different indices.

**c = 4, Q = Z2²: exact type structure.**
- Each column x ∈ {±1}^4 is either ±(a row of H₄) or a 3–1 vector.
- Equivalently, x̂ = ±4e_ψ (type e_ψ, 8 columns) or x̂ ∈ {±2}^4 with sign product −1 (type ★, 8 columns).
- Type entropy is 2 bits, and sign entropy given the type is 2 bits.
- **Given all types and the ★ signs, the four ψ-problems are independent 4-list problems on H.**
- The cost of fixing them is about 3.5 bits per column orbit (2 bits for types alone), i.e. 2r–3.5r bits.
- Solutions per type map ≈ E·2^−2r ≪ 1 for every configuration: 479 has 2r = 92 against log₂E 11.
- [M] In the real 479 core, n_★ = 77 + 61 + 29 + 52 = 219 ≡ v (mod 4), as forced by 4N₁ + n_★ = v. About half the
  columns are ★, as for random columns.

**Any Q with an index-2 subgroup (the c = 2 version of the same thing).**
- Split each block into U = X|Q₀×H and V = X|(Q₀+t)×H, and set P = (U + V)/2, R = (U − V)/2.
- These are ternary with complementary supports. The conditions become:
  - Σ_k |P̂_k|² = v − (v − b)·[φ ≠ 1, ω = 0];
  - Σ_k |R̂_k|² = v − (v − b)·[ω = 0].
- That is literally two T-matrix-like problems joined by supports.
- The support pattern carries 1 bit per orbit (2r in total), and the cross condition
  Σ_k [CC(U_k, V_k) + CC(V_k, U_k)] = 2(T_P − T_R) does not vanish because T_P ≠ T_R.
- The T-matrix rotation (P, R) ↦ (P ± R)/2 just returns (U, V). No free decoupling exists.

**What a k-tree would face per list if types were given.** Each ψ-problem would carry about n/4 bits (about 36–56). That
is easy for a k-tree, since lists of signs per block have 2^(r/8)–2^(r/4) elements. But a random type map has
E·2^(−2r) ≈ 2^(−80) expected solutions for 479, and worse elsewhere, so about 2^80 type maps would have to be tried.

For Q = Z4, the column transform takes values in {0, ±2, ±2i, ±2 ± 2i, ±4}. The coupling information per column is again
about 2 bits, so the same conclusion holds. For large c (571: 22 characters per column) the coupling information per
column is larger still.

**Conclusion [T]:** every exact decoupling found fixes type or support information worth ≈ 2r bits. It is affordable only
if log₂E > 2r, which holds for no b ≥ 3 configuration in range.

### 2.5 What structure can buy: memory, not time

- **Compression-first prefilter.**
  - Generate each block list directly from the exact compression fibre, by DP and sampling.
  - Lists shrink by the per-block compression information and n shrinks by n_comp, so T_min is unchanged and memory
    drops.
  - For 571 this is worth about 46 bits of n.
- **Representations through lag separability [T].**
  - For a subgroup K of index 2, the PAF at lags in K is additive over the two K-cosets of each block.
  - So lists of blocks with prescribed ℓ bits of K-lag PAF can be generated by an inner MITM over half-blocks in time
    2^(r/2) + output size.
  - Result: T·M² ≈ 2^(n−ℓ) for ℓ ≤ n_K (about n/2), but still T ≥ 2^(n − 2r).
  - This is Schroeppel–Shamir's memory saving, not a change of exponent.
  - For 479: T ≥ 2^78 at M = 2^30, against 2^85 without it.
- **Filters such as flatness** reduce n by about the number of log₂E bits they spend (487: Δn ≈ 47 for Δlog₂E ≈ 49).
  They help memory, but cannot go below n − log₂E: for 571, 168 − 7.8 ≈ 160 bits, so T ≥ 2^100 at M = 2^30.

### 2.6 Two-sided (CRT) compression [T]

Setup: H = Z_{h₁} × Z_{h₂}. Compression A goes onto Q × Z_{h₁} and compression B onto Q × Z_{h₂}. The lift handles
the mixed characters, i.e. ω₁ ≠ 0 and ω₂ ≠ 0.

- Let H(y | A, B) be the per-block entropy left once both marginals are known.
- Pairs (A, B) to try per solution: ρ = 2^(n_mix − 4H(y|A,B)).
- Each lift is an exhaustive 4-list problem with lists of size 2^H, costing 2^(2H).
- Total: ρ·2^(2H) = 2^(n_mix − 2H).
- Against the direct 2^(n − 2r) this differs by 2^(2I − n_AB), where I = r − H is the per-block information in the
  marginals and n_AB counts the compression constraint bits over all four blocks.

**CRT wins only if n_AB > 2I.** [E] Typically I is 25–40 bits per block while n_AB ≈ 0.2n ≈ 30–45, so it loses. This
explains the repo's experience:
- the lift given the true (A, B) of 347 took 4 s;
- the blind controls at v = 95, 131 and 215 found nothing: stage A produced 10^5.9–10^7.3 compressions, and none of
  the roughly 2.6·10^5 that were lifted exactly had a lift.

### 2.7 571 (c = 22, h = 25): direct answers to §5 Q2

- **H-compression first?** It carries only the characters (1, ω), about 1/22 of the lift bits (≈ 8 bits).
- **Is stage 2 separable per z?** No:
  - lags (q, z ≠ 0) pair column z′ with column z′ + z, so every ψ ≠ 1 condition couples all columns;
  - only the pure-Q lags (q, 0) are column-separable: a 4h-list "column set" condition with ≈ c/2 classes, a few bits.
- **The right first stage is the Q-compression.**
  - It is a strong per-block filter, worth 45.8 bits on the refined model.
  - After it, 168 bits remain in one coupled 4-list problem, so T ≥ 2^80.
- **Parseval number.** Σ_k ‖d_k‖² = 4b²/c = 80.2.
  - With h odd, x̄ is odd. For example, 76 entries at ±1 and 12 at ±3 give Σx̄² = 184 = 4(h + b).
  - Monomial compressions need Σ_k δ_k² = 84 with δ_k even: δ = (8, 4, 2, 0) or (6, 4, 4, 4).

### 2.8 Positive controls with known answers [M]

Use these for any new method before touching an open v:
- **c = 4 exact slices**, from `bge3_exact_count.py`, all with Q = Z2²:

  | v | h | d | r | L | solutions per ordered profile | profile |
  |---|---|---|---|---|---|---|
  | 23 | 5 | 2 | 11 | 8 | 1024 | (2,4,6,6) |
  | 31 | 7 | 3 | 10 | 5 | 96,624 | (2,2,4,10) |
  | 55 | 13 | 3 | 18 | 9 | 3.85·10⁹ | (2,2,4,14) |
  | 55 | 13 | 3 | 18 | 9 | 4.34·10⁹ | (2,4,10,10) |

- **Negative controls**, with exact answer 0:
  - v = 55, d = 3, r = 20: no admissible compression;
  - v = 63, d = 4, r = 19: blocked mod 32.
- **c = 6 rungs** v = 35 and 47 (repo).
- **The published 479 slice** (d = 12, r = 46): log₂E ≈ 11, T_min ≈ 2^71. A method that genuinely beats the barrier
  must find a solution here.

---

## 3. Q2: algebraic and cyclotomic constructions

### 3.1 Literature check (from memory; please verify citations)

- **Miyamoto (1991):** H(4q) for prime powers q ≡ 1 (mod 4) given an H(q − 1). All open prime v here are ≡ 3 (mod 4).
- **Xia–Xia–Seberry-type cyclotomy:** gives v = q². None of the open v is a square.
- **Paley I/II and Szekeres:** orders q + 1 and 2(q + 1). Not applicable.
- **Ito / cocyclic over dicyclic groups:** known only for small t (≲ 50).
- **The b ≥ 3 architecture itself:** I know of no theorem producing it. Structurally it is "an H(4c), inflated by h, with
  four rows replaced by a GS core". Note that 4v + 4 = 4c(h + 1): the border is the point "∞" of H ∪ {∞}, minus the
  four directions H₄ ⊗ 1_c.

### 3.2 Obstructions for the natural ansätze [T]

- **Rank-1 Legendre core.** Take X_k = w_k ⊗ χ + ε_k ⊗ δ₀ with (w_k) and (ε_k) both complementary over Q (power 4c)
  and amicable, which is the natural way to make the sum independent of ω.
  - It gives Σ_k |X̂_k|² = 4c(h + 1) = 4v + 4 at every ω ≠ 0, so it overshoots by exactly 4. This is the
    H(4c(h + 1)) "full border" object.
  - Its compression is ±1-valued, so Σ_k a_k² ≤ 4c² < 4v.
- **Paley rows.** Every row is constant on QR and on QNR (multiplier group of order (h − 1)/2, α = id).
  - Row sums then lie in {±1, ±(h − 2), ±h}.
  - Together with Σ_{k,q} x̄² = 4(h + b) this forces (1 + h)N_α + 4N′ = 4c(h + 1) − 4, which has no solution for h > 7.
  - Twisting by an involution α ∈ Aut(Q) gives a sub-case of this one, so it fails too.
- **General lemma (row-sum granularity).** Suppose S fixes row q and u has orbits of size d on H∖{0}.
  - Then x̄(q) ∈ ±1 + 2d·Z.
  - Near-constancy (Σ‖d_k‖² = 4b²/c) together with Σa² = 4v forces 2d ≲ 2√(4b²/c) + 2. Otherwise the admissible blocks
    have |a_k| ≈ c·x with x ∈ {1, 2d ± 1, …}, which misses 4v.
  - This is the mechanism that killed 659/c6/d9 (d = 9, row sums in {±1, ±17, ±19, …}).
  - Twisted α spreads each orbit over several rows and is the only way around it.
- **Legendre grid, c = 4, h = p ≡ 3 (mod 4).**
  - Columns ψ ≠ 1 are filled by one ±χ each, with a ★ at z = 0 in those three blocks. This meets the ψ ≠ 1 targets
    exactly, since the cross terms vanish because χ(−z) = −χ(z).
  - It then forces the fourth block to be 1_Q ⊗ B with PAF_B = pδ, a perfect binary sequence of odd length p > 1.
    None exists. Spreading column 1 across blocks conflicts with the row partition.
- **One algebraic block plus three searched.** If X₁ = w ⊗ χ + ε ⊗ δ₀, then |X̂₁(ψ, ω)|² takes only two values in ω ≠ 0
  (it depends on ω only through χ(ω)). The other three blocks must then satisfy a target that takes two values in ω.
  - This removes about c bits of freedom and no constraints, a net loss.
  - It also inherits the compression problem: x̄₁ = ±ε is ±1-valued.
- **"Half-core" (two Legendre blocks plus a free pair).**
  - The free pair must have power independent of ω, for every ψ.
  - Counting gives E ≈ 2^(c − L/2): it needs L ≤ 2c, which forces orbit sizes that the granularity lemma forbids.
  - For 719 (c = 10, 4), 599 (c = 10) the row-sum arithmetic already fails: 4v − 2c is not a sum of two squares.

### 3.3 Exhaustive high-symmetry slices [M]

`scripts/bge3_slice_exact.py` covers every symmetry group S ⊂ Aut(Q) × (Z/h)* with at most 2 generators and r within
the stated range. For each slice it applies, in order:
- row-sum profiles;
- the exact per-block PSD bound;
- an exact 4-list join, with matches re-verified.

| configuration | slices decided | r range | no profile | PSD kills every block | join empty |
|---|---|---|---|---|---|
| 571 c22 (Z22) | 8 | 12–20 | 0 | 8 | 0 |
| 599 c20 (Z2×Z10) | 82 | 8–22 | 4 | 55 | 23 |
| 599 c20 (Z20) | 20 | 13–24 | 0 | 15 | 5 |
| 647 c12 | 15 | 13–24 | 10 | 0 | 5 |
| 659 c6 | 8 | 10–24 | 8 | 0 | 0 |
| 719 c8 | 24 | 9–25 | 24 | 0 | 0 |
| 719 c10 | 8 | 8–24 | 7 | 1 | 0 |
| 739 c10 | 17 | 8–26 | 13 | 3 | 1 |

**No solutions in any of the 182 slices.** Positive control: the same code finds solutions in the v = 31 and v = 23 slices
(§2.8).

A cyclotomic construction with multiplier symmetry of index ≤ ~4 in (Z/h)* would have shown up here. Remaining algebraic
hopes would need r ≥ 27 (log₂E ≈ 0 territory) or symmetry outside Aut(Q) × (Z/h)*, such as block-permuting twists.

---

## 4. Q3: screens

### 4.1 The border (corner) condition [T/M]

- [T] **The corner is forced.** Let Cb be the compressed core, the 4c × 4c GS array of x̄ in the core's layout. Then
  C = −K·Cbᵀ·Kcᵀ/(4c). Verified on 347's template with maximum error 0 (`bge3_border_check.py template347`).
  - K and Kc must be ±1, with rows mutually orthogonal and balanced in each block.
  - Given that, C is automatically √(4b) times an orthogonal matrix.
  - **The only extra condition is that every entry of K·Cbᵀ·Kcᵀ equals ±4c.**
  - It depends only on the deviations d_k, and on **how** they sit in the GS layout.
- [M] **c = 4 (Q = Z2²).** Every compression is "3 single ±2 bumps + 1 constant block" (§4.4). All 2048 placements and
  signs have a border, with C checked Hadamard. **No restriction.**
- [M] **c = 6.** There are four compression profile types, by c·‖d_k‖² per block: (0,8,8,9) [347's], (0,5,8,12),
  (0,0,8,17) and (0,0,5,20).
  - With frames derived from 347's K (random within-block permutations and signs, Kc solved exactly over 2^19 sign
    vectors), borders were found **only for 347's own placement**.
  - Of 150 copies of 347's compression transformed by core symmetries (per-block translation, negation, reversal, block
    permutation), at most 1 works with 347's frame: 1/150 and 0/150 in two runs with different seeds.
  - For 6 sampled admissible compressions of each of the four profile types, with 20 frames each, no border was found
    (`logs/bge3_answer/border.log`). That includes other placements of 347's own profile type (0,8,8,9).
  - [E] A random frame satisfies all (4b)² entry conditions with probability ~10^−15. **Borders must be built together
    with the compression.**
- **Consequence.** Before any search with c ≥ 6:
  - fix the compression to a placement with a known border, as the c = 6 rungs did with 347's pattern;
  - for c ≥ 10 (571, 599, 647/c18, 739, 643, 573, 719/c20) **a border has to be constructed first.**
  - For c = 22 the frames cannot be Kronecker products: there are no 21 orthogonal balanced ±1 vectors of length 22.

### 4.2 Exact compression DP and conditioned null directions (`bge3_cost2.py`) [M]

This is equivalent to the t = 1 screen with profile restriction, but faster because it works on class-orbits. It also
reports the conditioned lift. It decided configurations the prompt lists as undecided or queued:

| configuration | status in the prompt as posed | result |
|---|---|---|
| 647 c12 Z12 d13 r60 L30 | "t = 1 undecided (re-running)" | **no admissible compression quadruple** (2698 single-block compressions) |
| 659 c6 Z6 d9 r78 L39 | "queued" | **none**: row sums in {±1, ±17, ±19, …}, see §3.2 |
| 599 c20 Z2×Z10 d14 r55 L30, groups 0–3 | "fails t = 1 (all 6 groups)" | confirmed: none |
| 599 c20 Z20 d14 r60 L32, group 0 | "one fails, one undecided" | group 0: none; group 1: my DP too large (> 40 min), but the norm-budget screen (commit 4888fd1, another session) decides it infeasible |
| 647 c18 d12 r74 L39 | queued | DP overflow (5.4·10⁶ states) |
| 571 c22 d10 r62 L35, groups 0–3 | "t = 1 passes; t = 5 re-running" | **admissible**: 60,384 quadruples, n_lift 168, log₂E 23.5 |
| 719 c8 d22 r45 L23; 659 c6 d18 r52 L27; 739 c10 d18 r54 L29 | "fails t = 1" | confirmed: no admissible compression quadruple (cross-check of the two screens) |

### 4.2b Norm-budget screen (another session, commits 4888fd1–1d0d968)

`bge3_normscreen.py` uses the exact identity Σ_k ‖x_k‖² = (4v(ct − c + 1) + 4b(c − 1))/(ct) for compressions onto
Q × Z_t. It is weaker than the full t = 1 DP but much faster.
- **Where they overlap it agrees with every verdict here:** 599 c20 dead (all groups); 647 c12, 659 c6, 719 c8 and
  c10, 739 c10 dead; 571 c22 survives.
- **Its survivors with n ≤ 300:** 571 c22 (n 192), 719 c20 (227), 523 c4 group 1 (226), 643 c14 (232), 647 c18 (233),
  659 c10 (242), 719 c18 (251), 573 c14 (269), 659 c4 Z4 r90 (288) and 719 c6 (289).
- **Ranking.** 571 c22 is the cheapest survivor. Every other survivor has n ≥ 226, so by §2.1 it needs
  ≳ 2^((n − log₂E)/2): at least 2^93 when log₂E ≤ 40.
  - The known values are: 643 c14 ≈ 2^102 (refined), 647 c18 ≈ 2^98, 573 c14 ≈ 2^130.
  - 523 c4 group 1 has refined log₂E ≈ −32, i.e. no expected solutions.

### 4.3 Row-sum granularity [T]

See §3.2. Cheap check before costing: if a fixed row's orbits all have size ≥ d with 2d > 2√(4b²/c) + 2, the
configuration is dead.

### 4.4 Forced compression shapes [T]

- **c = 4, Q = Z2².**
  - h odd: three blocks have x̄ = μ + 2e_q (one bump each) and one block is constant.
  - h even: the same shape, with a_k ≡ 2 (mod 4) in three blocks and ≡ 0 (mod 4) in one. Hence 4v ≡ 12 (mod 16), so
    v ≡ 3 (mod 4).
  - Proof: x̂(ψ) ≡ a (mod 4) and x̂ ∈ 2Z, and Σ_k x̂_k(ψ)² = 12 must be 4 + 4 + 4 + 0.
- **Monomial compressions (any c)** need Σ_k δ_k² = 4b with δ_k even. For example:
  - c = 22: δ = (8,4,2,0) or (6,4,4,4);
  - c = 6: (4,2,0,0), the profile (0,0,5,20).
  - By per-block translation all bumps can be placed on one coset (359's "same special coset" is this normalization).

### 4.5 Lattice (residue) obstructions [M]

The calibration slice v = 63/d4/r19 has admissible compressions and a refined-model estimate of 2^9.1. It has **0
solutions**:
- the target is reachable mod 8 (52/52 row-sum orderings);
- mod 16 for 16/52;
- **mod 32 for 0/52** (`logs/bge3_answer/residue_probe.log`).

Individual lag classes have PAF residues fixed by the row sum. For example, class 15 has PAF ≡ 0 (mod 16) for a ∈ {2, 14}
and ≡ 12 for a ∈ {4, 12}. The Gaussian model assumes spacing 4 everywhere, so it cannot see this.

The repo's joint forced-residue filter (`bge3_residue.py`) is the scalable form. **Run it per configuration, against each
admissible compression quadruple from `bge3_cost2.py`, before any lift.**

---

## 5. Q4: forensics

- **CRT-composite h in all three.**
  - h values: 479 = 7·17, 347 = 3·19, 359 = 4·11 (H = 2Z₈ × Z₁₁).
  - The cores avoid prime h and 2p, consistent with a two-sided (CRT) compression method.
  - [M] Compressions onto Q × Z_{h₁} and Q × Z_{h₂} look random: entry variance within 0.5–1.9σ of a permuted null.
  - Against random blocks with the same symmetry they are *less* lumpy (z = −1.3 to −3.2). That follows from the core
    condition itself (power exactly 4v at those characters), so it is not evidence of a construction
    (`bge3_forensics3.py`).
  - [T] §2.6: two-sided compression does not beat the barrier by itself.
- **Compression placement.**
  - 479 and 359 are monomial; 347 uses small sets ({Z₂×0}, {(0,1),(0,2)}, {(0,0),(1,1),(1,2)}).
  - Border existence depends on placement (§4.1). The "same special coset" in 359 is the placement in which GS(d)
    restricted to V⊥ becomes G₄^±(δ) ⊗ I on the R-even and R-odd parts, which is the natural one for building K and Kc
    by hand.
  - [E] The authors probably chose compression and border first, then searched cores in that fibre, as the rung search
    does.
- **479's twist.**
  - α of order 3 on Z2² leaves 2 class-orbits, so x̄ = (x₀, x₁, x₁, x₁): **one parameter per block**, and the compression
    costs 10 bits instead of 37 for 347.
  - It also identifies the three ψ ≠ 1 problems up to the multiplier u, i.e. Y_{k,ψ∘α} = Y_{k,ψ}∘u. This roughly halves
    the free grid. It is already inside the cost model (n ≈ 144–151).
- **No product or low-rank structure.**
  - Each (block, coset) slice of the cores, as an h₁ × h₂ matrix, carries more top-singular-value energy than i.i.d.
    signs: 0.42 vs 0.29 ± 0.01 for 479, 0.63 vs 0.46 ± 0.01 for 347, i.e. 15–16σ.
  - **Against the symmetry-aware null** (random blocks on the same affine orbits) this vanishes: z = −0.2 and +0.8. The
    multiplier symmetry is the whole effect.
- **Cost.** Under the calibrated model 479 required ≈ 2^71 and 347 ≈ 2^78 of list-method work (§2.1).
  - Either the authors had a budget of 10^21–10^23 operations, or they used a non-list method.
  - The repo's local search scales as about 2^(0.55n), i.e. 2^80+ here, so it is not that.
  - What their method is remains the open question.
  - **Cheapest experiment to resolve it:** reproduce a solution inside 479's slice with any candidate method.

---

## 6. A generalization of the architecture (b < c − 1): derived, not yet realized [T]

The border derivation does not need Ψ = Q̂∖{1}.
- Take any conjugation-closed Ψ ⊂ Q̂∖{1} with |Ψ| = b odd and v = ch + b.
- Target: Σ_k |X̂_k(ψ, ω)|² = 4v − 4m·[ψ ∈ Ψ, ω = 0], i.e. Σ_k PAF_k(q, z) = 4vδ − 4·Σ_{ψ∈Ψ} ψ(q).
- K and Kc must have every block-part a ±1 vector with Fourier support in Ψ: KᵀK = 4c·(I₄ ⊗ P_Ψ).
- The corner formula is unchanged.

The original architecture is the case Ψ = Q̂∖{1}.

For Q = Z2³ and Ψ a basis (c = 8, b = 3, v = 8h + 3):
- the frame problem is literally the c = 4 one, so borders always exist;
- the family covers 523 (h = 65), 571 (71), 643 (80), 659 (82) and 739 (92);
- other Q and Ψ give options for v = 537, which has no standard border, e.g. c = 16, b = 9, h = 33.

**But [M]:**
- no core exists at h = 1 (parity proof: the block with a ≡ 0 mod 4 would have zero power everywhere);
- none exists in the z ↦ −z slice at h = 3 or h = 5 (exact, all row-sum profiles).

Treat it as an unconfirmed lead. The next step is a Gaussian cost run with this target before any search.

---

## 7. What to implement first

1. **Screen every remaining candidate in this order:** `bge3_normscreen.py` (fastest, exact), then `bge3_cost2.py`
   (exact compression, calibrated lift), then the joint-residue filter. Candidates: the norm-screen survivors 571 c22,
   643 c14, 647 c18, 719 c20, 719 c18, 659 c10, 573 c14 and 523 c4 group 1.
   - `bge3_cost2.py` is exact on the compression and calibrated to about 1.5 bits on the lift.
   - Together the screens have already removed 647/c12, 659/c6, 599/c20, 719/c8, 719/c10 and 739/c10.
2. **Border existence for c ≥ 10** before any core search there. Starting point:
   - monomial compressions (δ = (8,4,2,0) or (6,4,4,4) for c = 22) in the 359 placement;
   - search frames by the ε-enumeration method in `bge3_border_check.py`.
   - If no border exists for any admissible c = 22 compression, 571's b ≥ 3 route is closed regardless of cost.
3. **Do not spend Kaggle or local compute on b ≥ 3 for open v.** The best surviving configuration, 571, needs ≳ 2^80
   (≳ 2^100 at 2^30 memory). 523/b = 0/d = 9 with flat pools (n′ ≈ 91, k-tree ready) remains the best route in range.
4. **To learn the published method,** attack 479's own slice (d = 12, r = 46, log₂E ≈ 11) with candidate methods. A
   method that finds it in ≲ 2^40 is the breakthrough the prompt is asking for; the c = 4 calibration slices (§2.8) are
   the smaller rungs.
5. **Optional, cheap:** extend the exhaustive slice decision to r ≤ 30 with a C join, which would close the algebraic
   question for larger cyclotomic indices.

---

## 8. Files

| file | purpose |
|---|---|
| `scripts/bge3_cost2.py` | refined cost model: exact compression DP plus conditioned Gaussian (fibre sampler with kernel moves) |
| `scripts/bge3_exact_small.py` | exact total solution count, small r |
| `scripts/bge3_exact_count.py` | exact per-profile count with PSD prefilter, r ≤ 24 |
| `scripts/bge3_slice_exact.py` | exact decision of high-symmetry slices (positive control in §2.8) |
| `scripts/bge3_border_check.py` | corner formula check; c = 4 exhaustive border existence; c = 6 frame search |
| `scripts/bge3_residue_probe.py` | 2-adic reachability of the target, small slices (§4.5) |
| `scripts/bge3_forensics3.py` | CRT-compression and rank forensics with symmetry-aware nulls (§5) |
| `scripts/bge3_genpsi_exact.py` | generalized-Ψ architecture, exact small search (§6) |
| `logs/bge3_answer/calib_exact.jsonl`, `calib_gauss.jsonl` | calibration (§2.2) |
| `logs/bge3_answer/cost2_open.jsonl`, `cost2_rerun.jsonl` | refined model on open and published configurations |
| `logs/bge3_answer/slices_exact.jsonl` | 182 slice decisions (§3.3) |
| `logs/bge3_answer/border.log`, `residue_probe.log`, `forensics3.log`, `genpsi_exact.log` | outputs of the above |
