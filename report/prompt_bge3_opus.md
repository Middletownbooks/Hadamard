# Prompt for a higher-reasoning session: b ≥ 3 bordered Goethals–Seidel cores for open Hadamard orders 2000–3000

*Self-contained. If you are running in Claude Code, the repository is `middletownbooks/hadamard`, branch
`claude/focused-planck-j1645f`. The file paths cited below exist there. Everything stated as fact is either a theorem or
was verified in exact integer arithmetic. Estimates are marked [E].*

---

## 0. Goal and what I want back

Hadamard matrices of order 4v are unknown for v ∈ {523, 537, 571, 573, 599, 643, 647, 659, 669, 719, 739}, i.e. orders
2092–2956. In 2026 one group (Alpöge, Voinov, Reynolds-Haertle, with Claude) closed every open order below 2000. Three of
their matrices (v = 347, 359, 479) use a **bordered architecture with b ≥ 3 border rows**, described in §1. For several
of the open v, this architecture's cheapest symmetric configurations are estimated at **n ≈ 180–230 constraint bits**
(§3). That is the same range as the published 347 (n ≈ 201) and 479 (n ≈ 151). Our own search methods stall at about
**90 bits** (§4).

**I want one or more of the following, ranked by how likely they are to produce a matrix:**

1. **An algorithm** that exploits the structure in §2 to bring the effective search cost for at least one open
   configuration down to roughly 2^50 elementary operations, i.e. within reach of about 20 CPU-cores × days. For each
   proposal give:
   - a complexity estimate;
   - the decomposition into independent lists or subproblems;
   - a positive control: a smaller instance where the method must provably or demonstrably find a solution, ideally
     one of the published cores or the small "ladder rungs" in §4.
2. **An algebraic or cyclotomic construction** for the two-valued target in §1. Look especially at configurations where
   h (the order of H) is a prime or prime power:
   - 571: c = 22, h = 25 = 5²
   - 599: c = 20, h = 29
   - 647: c = 12, h = 53
   - 659: c = 6, h = 109
   - 739: c = 10, h = 73
   - 719: c = 8, h = 89; c = 10, h = 71

   Possible frameworks: relative or divisible difference families, Paley/Szekeres-type cyclotomy, Xia–Xia–Seberry,
   Miyamoto, or Ito-type constructions. Say whether any existing theorem already covers one of these (v, c, h).
3. **Stronger cheap necessary conditions** (screens) that can kill configurations. The current exact screen only checks
   the compression onto Q × Z_t (§3), which is weak when h is prime.
4. **Forensics:** what the published cores' extra structure (§1.4) suggests about how they were found, if it implies a
   method that scales.

Prefer concrete, checkable statements: formulas, pseudocode, and explicit small cases. Flag guesses as guesses.

---

## 1. The b ≥ 3 architecture (verified on the three published matrices)

### 1.1 Parameters

- b odd ≥ 3, and c = b + 1 (even).
- v + 1 = c(h + 1), so h = (v+1)/c − 1, and m = v − b = c·h.
- G is an abelian group of order m with a subgroup H of order h and index c. In the split case, used throughout,
  G = Q × H with |Q| = c and H cyclic.

| v | b | c | h | G | H | G/H |
|---|---|---|---|---|---|---|
| 479 | 3 | 4 | 119 = 7·17 | Z2²×Z119 | Z119 | Z2² (split) |
| 347 | 5 | 6 | 57 = 3·19 | Z6×Z57 | Z57 | Z6 (split) |
| 359 | 7 | 8 | 44 = 4·11 | Z2²×Z8×Z11 | Z44 | Z2³ (**non-split**) |

### 1.2 The core

Four ±1 functions X₁..X₄ on G (group-developed blocks of a Goethals–Seidel array) with

    Σ_k PAF_{X_k}(g) = 4m            g = 0
                     = −4b           g ∈ H \ {0}
                     = +4            g ∉ H
    Σ_k a_k² = 4v,  where a_k = Σ_g X_k(g).

The core sits inside a bordered GS array; the border is described next.

### 1.3 The border (understood; not the bottleneck)

- The 4b × 4b corner is a Hadamard matrix of order 4b.
- Each border row, restricted to a core block, is constant on H-cosets. So the border is a 4b × 4c sign matrix K,
  inflated by h.
- K is a Hadamard matrix of order 4c with the 4 rows H₄ ⊗ 1_c removed; equivalently, every border row is balanced
  inside each block.
- The template for c = 6 (from 347) is in `data/template347.npz`.

**Given any core satisfying §1.2 for some v, the full Hadamard matrix of order 4v follows.** This is verified for the
c = 6 rungs v = 35 and v = 47, which we found by annealing and completed with 347's border.

### 1.4 Extra structure observed in the published cores

- **479:** a twisted automorphism stabiliser of order 12 inside GL(2,2) × (Z/119)*. It permutes the three nontrivial
  Z2²-cosets together with a multiplier, giving 46 affine orbits per block, against 128 with scalar multipliers only.
- **347:** scalar multipliers only (order 6, 70 orbits).
- **359:** scalar order 5 on Z11 only (96 orbits). All four blocks deviate from constant on the **same single coset**
  of G/H. There is also an exact twisted slice relation, D(1,0,e,f) = (−1)^e · D(0,1,e,f).
- **Compressions onto Q** (coset sums) are "constant per block plus a small deviation". This is forced (§2).
- **Flatness:** b ≥ 3 blocks are statistically uniform (PAF spread is null-like). That is unlike the group's b = 0
  cores, which are unusually flat.

---

## 2. Character decomposition over Q (derived, and verified exactly on 479 and 347)

For a character ψ of Q define the component on H:

    Y_{k,ψ}(z) = Σ_{q∈Q} ψ(q) X_k(q, z).

Taking the Fourier transform over Q of the target in §1.2 gives T̂_ψ(s) = −4c for s ≠ 0 when ψ ≠ 1, and 0 when ψ = 1.
With Ŷ the Fourier transform over H, the core condition is equivalent to:

- **ψ = 1** (compression onto H):
  Σ_k |Ŷ_{k,1}(ω)|² = 4v for every ω ∈ Ĥ, including ω = 0.
  So the four H-compressions form a **complementary integer quadruple over Z_h**, with entries in {−c, −c+2, …, c}.
- **ψ ≠ 1:**
  Σ_k |Ŷ_{k,ψ}(ω)|² = 4v for ω ≠ 0, and = 4b at ω = 0.

Measured values:

| core | ψ = 1, all ω | ψ ≠ 1, ω = 0 | ψ ≠ 1, ω ≠ 0 |
|---|---|---|---|
| 479 | 1916 (= 4v) | 12 (= 4b) | 1916 (= 4v) |
| 347 | 1388 (= 4v) | 20 (= 4b) | 1388 (= 4v) |

**The coupling:** for each block k and each z ∈ H, the vector (Y_{k,ψ}(z))_ψ is the Q-Fourier transform of a column in
{±1}^c. So the c sub-problems over H, each of length h, are coupled only position by position, through the set
F_Q({±1}^c) of 2^c allowed patterns.

**Example, c = 4 with Q = Z2²:**
- 3–1 columns give all four components equal to ±2, under one parity constraint.
- 2–2 or 4–0 columns give a single ±4 component and zeros elsewhere.
- So each (k, z) is a "type" (±4 at one character, or all ±2) plus signs.

**A precedent that this kind of coupling can be broken.** T-matrices (four disjoint ternary circulants with ΣPAF = t·δ)
are equivalent to GS quadruples with the pointwise product condition a·b·c·d = +1. A support split made them tractable:
1. Fix S = supp T₁ ∪ supp T₂.
2. Write T₁, T₂ = (P ± R)/2 on S, and likewise on Sᶜ.
3. The condition becomes PAF(P) + PAF(R) + PAF(P′) + PAF(R′) = 2t·δ, with **four independent lists**.
4. A k-tree (Wagner) 4-list join then found T(103) with K = {1, 46, 56} **blind in 12 minutes** (lists of 2^17 per
   side, about 63 constraint bits).

**Question:** is there an analogous type split for the per-position coupling above that turns a b ≥ 3 core into
independent lists?

---

## 3. Open configurations: costs and screen status

**Cost model** [E] (`scripts/bge3_cost.py`):
- Blocks are invariant under a cyclic symmetry group S = ⟨(α, u)⟩ ≤ Aut(Q) × (Z/h)*, with r orbits per block and L
  lag classes.
- n = −log₂ P(4 random S-invariant blocks with the row-sum profile hit the target), from a Gaussian model.
- log₂E = Σ log₂(#blocks with row sum a_k) − n.
- **Calibration:** the real 479 configuration comes out at n ≈ 151 and the real 347 at n ≈ 201.

**Exact screen** (`scripts/bge3_screen.py`, profile-restricted version `bge3_screen2.py`):
- Compress onto Q × Z_t for t | h.
- Each block's image x̄ is an integer vector constant on class orbits.
- Require Σ_k |F_k(χ)|² = 4v if χ is nontrivial on Z_t, 4b if χ is trivial on Z_t but nontrivial on Q, and
  Σ a_k² = 4v at χ = 0.
- Enumerate achievable x̄ per block (a DP over orbits, pruned by |F| bounds), then meet in the middle over the 4 blocks.
- t = 1 is just the Q-compression, which is weak.

**Screen status after this session** (exact DP screen with profile restriction, plus a new exact norm-budget screen:
Σ_k ‖x̄_k‖² = (4v(ct − c + 1) + 4b(c − 1))/(ct) for the compression onto Q × Z_t, verified on 479 and 347). Every
configuration with log₂E > 0 and n ≤ 300 was screened.

| v | surviving configurations (n bits [E]) | dead at or below 300 bits |
|---|---|---|
| 571 | **c = 22, Q = Z22, h = 25, d = 10, r = 62, L = 35 (192)**: passes t = 1 and the t = 5 norm test | c = 4 (all) |
| 523 | c = 4, Q = Z4, h = 130, d = 12, r = 64, L = 39 (226): passes t = 1, t = 5 and the norm test at every t | – |
| 719 | c = 20, Q = Z2×Z10, h = 35, d = 12, r = 66, L = 41 (227); c = 6 (289) | c = 8, c = 10, c = 18 |
| 643 | c = 14, Q = Z14, h = 45, d = 12, r = 72 or 74, L = 39 (232–234) | – |
| 659 | c = 10, Q = Z10, h = 65, d = 12, r = 70, L = 41 (242); c = 4, r = 90 (288) | c = 6 (all) |
| 573 | c = 14, Q = Z14, h = 40, d = 12, r = 76 or 80 (269–285) | – |
| 599, 647, 739 | **none** | every b ≥ 3 configuration with n ≤ 300 |

For comparison, the published 479 core sits at n ≈ 151 and 347 at n ≈ 201.

Other notes:
- 537 has no border (v + 1 = 2·269).
- 669's cheapest b ≥ 3 option is at n ≈ 590.
- **Twisted symmetries** (multipliers that also permute or negate blocks or coset classes) were searched only for the
  published cores' forensics, not costed for the open v. The 479 core shows that twisted symmetry is used.

---

## 4. What has been tried, and what we can do computationally

### 4.1 Compression-and-lift for b ≥ 3 (c = 6)

The pipeline: stage A enumerates compressions onto Q × Z_p, a forced-residue filter is applied, then a one- or two-sided
exact lift follows (`scripts/bge3_*`, `posctrl/RESULTS.md`).
- Validated on the real 347 core: the true compression is found and lifts back to the core.
- **Failed its blind positive control:** no solution at v = 95, 131 or 215 in any slice, partly for lack of compute.
- On v = 515 (u = 59 and u = 161 configurations) it found nothing.

### 4.2 Local search

The repo's local-search engines give:
- **c = 6 ladder by annealing:** solved h = 5 (v = 35) and h = 7 (v = 47) quickly. h = 9 and 11 stalled. Symmetric
  attempts at h = 57 (347's own size) and h = 85 never got below energy of about 10⁴.
- **Measured scaling for b = 0 cores:** roughly 2^(0.55·n) flips (tabu). So anything above about 70 bits is out of
  reach this way.

### 4.3 k-tree (Wagner 4-list) joins

These work up to about 83–91 "flat" bits when the four lists are independent:
- v = 487, b = 0, d = 9: found two new cores from flat pools of 2^26.
- v = 523: about 0.3 expected hits accumulated so far, none found (Kaggle runs continuing).
- T(103) via the support split: 63 bits, solved blind in 12 minutes.

### 4.4 Ruled out today

- **T-matrices with Cooper–Wallis:** 573 is dead (T-form row sums), 537 has no subgroup, 669 is at 157 bits.
- **Twisted T-ansatz:** gives nothing.
- **Twisted b = 0 for 571 (d = 15, block 3-cycle):** feasible at 94.6 bits, but it collapses to a lopsided 2-list
  problem (2^35 × 2^110).

### 4.5 Resources

- 4 local cores and 15 GB of RAM.
- About 5 parallel Kaggle notebooks (4 CPUs, about 30 GB each, 12-hour sessions).

---

## 5. Concrete questions

1. **Decomposition algorithm.** Using §2, can a b ≥ 3 core be split into independent lists, for example by fixing a
   per-position "type" (which Q-characters carry the column's mass) as in the T-matrix support split? Work out the
   c = 4 case (Q = Z4 or Z2², 523 at h = 130) and the large-c, small-h cases (571: c = 22, h = 25; 599: c = 20, h = 29).
   - How many constraint bits would a k-tree face per list?
   - What would the positive control be? (The published 479 core is c = 4; the rungs v = 35, 47 are c = 6.)
2. **Large c, small h.** At 571, Q = Z22 and h = 25, so G ≅ Z550 is cyclic. The ψ = 1 subproblem is a complementary
   integer quadruple of length 25 with entries in {−22, …, 22}, and the Q-compressions are nearly flat (Parseval forces
   Σ_k ‖x̄_k − mean‖² ≈ 4b(c−1)/c ≈ 80).
   - Does a two-stage method work here: first the H-compression quadruple (small), then the per-z columns?
   - Is the second stage a set of independent subproblems once the first is fixed?
3. **Algebra.** For h a prime power q (25, 29, 53, 73, 89, 109), is there a cyclotomic or relative-difference-family
   construction of four ±1 functions on Z_c × GF(q) with the two-valued PAF (−4b on {0} × GF(q)*, +4 off it)? Classical
   constructions might give one or more blocks, with the rest found by search, which would cut the bit count.
4. **Screens.** Which exact necessary conditions beyond Q × Z_t compressions are cheap? For example: per-ψ conditions
   at ω = 0 (Σ_k |Ŷ_{k,ψ}(0)|² = 4b for each ψ ≠ 1), integrality and parity of the components, a Parseval budget per
   ψ, sum-of-squares conditions per character orbit, or the border's balance condition.
5. **Forensics.** 479's twisted automorphism (order 12, 46 orbits), 359's "same special coset for all blocks" and its
   twisted slice relation: do these point to a construction method, such as an ansatz that fixes the compression and
   searches the residual?

Please return:
- for each idea, a short derivation;
- the predicted cost in bits or operations;
- the smallest test instance with a known answer;
- what to implement first.
