# A fresh look at the open Hadamard and skew-Hadamard orders: what was tried, what is now closed

Date: 2026-11-09 to 2026-11-12. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_fable_fresh.md`.
Scripts: `theory/fable_*.py`. Logs and working notes: `logs/fable_fresh/` (`NOTES.md` is the session's state file).

**Tags.** [P] proved here (proof in the text); [M] exact computation (script named); [E] estimate. A cap or a timeout
is not a zero.

**Credit.** Egan, Kotsireas, London and Ó Catháin (EKLO) for the survey and Theorem 7; Alpöge, Voinov and
Reynolds-Haertle for the twelve matrices of August 2026; K. Kumar for his constructions, his border formulation
(`search/NOTES.md` §102 of his repository) and `kt4.c`; the OpenAI model for Theorem B and the compression normal
form; this repository's earlier sessions for the skew k-tree pipeline that §5 reuses.

---

## 0. Bottom line

<!-- BOTTOM -->
1. **No new matrix.** No Hadamard matrix of an open order and no skew-Hadamard matrix of an open census value was
   built.
2. **Two natural families of product theorems are closed, with proofs or exhaustive computation** (§3, §4). Both
   were aimed at the open orders 4·t·v with a small odd t, where Goethals–Seidel (GS) quadruples of order v already
   exist (v = 223, 283, 491 from Alpöge et al., 571 from Kumar): 2676, 4460, 5892, 6852 and 7924.
   - *Multiplying one GS quadruple by ℤ_t* (every fibre {τ} × G carries ±X or ±X\* for a block X of the quadruple):
     for t = 3 and t = 5 such a formal identity exists **only** for quadruples with at most four real spectral
     degrees of freedom (Williamson type, "propus-like", periodic-Golay-like), and for **no** type with five or
     more [M]. In the product form T ⊗ W with arbitrary T-matrices it exists for the Williamson type only, for
     every t [M]. Types with four degrees of freedom have 2v sign bits and stop existing near v ≈ 160 by the usual
     count [E], so nothing reaches v ≥ 167.
   - *Multiplying two GS quadruples* (orders t and v, an octonion-type composition, output of order tv): all 2048
     monomial compositions were enumerated [M]. Each forces, on both factors, a family of ternary sequences with
     symmetric, disjoint supports and zero total autocorrelation, which cannot exist at an odd order above 1 [P].
3. **EKLO's Theorem 7 with a large index k and a tiny subgroup S is closed as a route** (§6). For S = ℤ₃ the core is
   equivalent to a four-block conference matrix over the sixth roots of unity with a sign condition [P]; Paley's
   conference matrices, quadratic or sextic, satisfy the conference equation and cannot satisfy the sign
   condition [P]; and the border exists exactly when two ±1 bases are "Hadamard-unbiased" through the compressed
   core [P], which a ±1-valued compression has no reason to allow [E]. This explains why borders are
   known only for spike compressions.
4. **Bordered arrays cannot reach the skew census** (§6.4) [P]. From the shape alone, a coset-constant border of a
   Goethals–Seidel core is one of Kumar's character borders, and when the group has odd order (as a skew-type block
   requires) its width is 4·odd, so the order is ≡ 0 (mod 8). The 630 open values 4v, v odd, can only come from
   unbordered arrays over a group of order v or from a different kind of array.
5. **A family of bordered arrays that the censuses have not used** (§7). The same shape analysis shows that over
   Q = ℤ₂^a *any* set T of characters may carry the border, not only T = {1} (Theorem 5), all nontrivial characters
   (Theorem 7) or all characters; the border then depends only on |T| integer quadruples [P]. One instance outside
   the classical shapes was built and verified: four affinely independent characters of ℤ₂³ over ℤ₂³ × ℤ₃, a
   Hadamard matrix of the (known) order 112 [M]. For 4p the new shapes start at ℤ₂⁴ (|T| = 7 or 11) and exist at
   every open prime. Two first instances (599 and 751) die at the Q-level for every T [P]; no census was run.
6. **The idle cores reran the one calibrated slice this repository had left unfinished** (skew 361, split torus
   over ℤ₁₉²) with fresh lists. Status in §8.
7. **Reading.** Every structure cheap enough for 4 cores is now either searched or closed by a theorem of the kind
   above; Kumar's own ranking puts the nearest live targets at about 10³ core-hours (4792 = 8·599 by a Spence
   plug-in; skew 661, 883, 1483). §9 says what I would run with that much compute.
<!-- /BOTTOM -->

---

## 1. Novelty sources checked

- Kumar's repository `vatskishan03/hadamard` is still at commit d885022 (2026-11-06) [M: `git fetch`, 2026-11-10].
  The open lists of `data/open_orders_lt10000.json` (85 orders) and `data/skew_open_census.json` (630 values) are
  unchanged.
- arXiv and Zenodo are not reachable from this container; nothing newer than the lists was checked there.

## 2. Ideas considered

| # | Idea | Decisive test | Outcome |
|---|---|---|---|
| 1 | Multiply a GS quadruple of order v by ℤ_t (reach 4·t·v from the cores at 223, 283, 491, 571) | exhaustive formal search, t = 3, 5 | closed [M], §3 |
| 2 | Octonion-type product of two GS quadruples | classification of all monomial compositions + parity lemma | closed [M, P], §4 |
| 3 | Williamson pattern with T-matrices that satisfy three "half-amicability" relations (one free block, three symmetric) | exhaustive over cyclic T-matrices, t ≤ 11 | none exist [M], §3.4 |
| 4 | Theorem 7 with large k, S = ℤ₃, F₅ or F₇ (p = 4k − 1, 6k − 1, 8k − 1) | cell-energy identity, sixth-root reformulation, border analysis | closed [P, E], §6 |
| 5 | The torus base point: G = F_p\* (the Theorem 5 group of order p − 1), blocks χ(x + c) | exact autocorrelation, energy count | no nearby solution [P, E], §5.2 |
| 6 | Norm-equation form over the fixed field of the multiplier group | reformulation only | no algorithm gained [P], §5.1 |
| 7 | Uniform cyclotomy (punctured lines of F_q²) in the fibres of Theorem 7 | energy identity | impossible for q ≥ 7 [P], §6.5 |
| 8 | Ito / dicyclic matrices (negaperiodic Golay pairs of length 2p) | sign count | negative for p > about 160 [E]; dropped |
| 9 | Legendre pairs over a group of order 2p − 1 | sign count; Kumar screened the four-block form | dropped |
| 10 | 10012 = 4·2503 over ℤ₂² × F₆₂₅ (the EKLO family at q = 25) | Kumar's notes | already examined by Kumar (F₂₅ × F₂₅ exactly empty; "2 % per run" otherwise); dropped |
| 11 | 4792 = 8·599 by Kumar's Spence plug-in ℤ₇ × ℤ₁₉ × F₉ | his calibrated cost | 1.5·10³ core-hours against 288 available; dropped |
| 12 | Skew T-structures for Paley blocks beyond t = 3 | Kumar's exhaustive t = 9, 19; count at t = 33 | about 10⁻³ expected [E]; dropped |
| 13 | Unstructured small skew values (89, 101, 107, 149) | size | 2³¹¹ sign vectors at 89 with no multiplier; dropped |

Ideas 8–13 were dropped on the stated one-line reason and are not discussed further.

---

## 3. Multiplying a Goethals–Seidel quadruple by ℤ_t

### 3.1 Definition

Let (A, B, C, D) be ±1 sequences on an abelian group G of order v with
PSD_A + PSD_B + PSD_C + PSD_D = 4v at every nontrivial character. Write X\*(g) = X(−g).

- A **formal multiplication by ℤ_t** is an assignment of a *letter* u(j, τ) ∈ {±A, ±B, ±C, ±D, ±A\*, ±B\*, ±C\*, ±D\*}
  to each slot (j, τ), j = 1..4, τ ∈ ℤ_t. It defines four ±1 sequences Z_j(τ, g) = u(j, τ)(g) on ℤ_t × G.
- It is **valid for a type** if Σ_j PSD_{Z_j} = 4tv at every nontrivial character of ℤ_t × G for *every* quadruple
  of that type. A **type** is a set of relations among the letters: each block is free (f), symmetric (s: X\* = X),
  or equal to another block or to its reverse. Its **dimension** is 2·#free + #symmetric, the number of real
  parameters of the spectral vector y = (Â, B̂, Ĉ, D̂)(ψ).
- At a character (ρ, ψ) the transform of Z_j is Σ_τ ρ^τ u(j, τ), with u = ±y_k or ±conj(y_k). Validity for all
  ρ ∈ μ_t is equivalent to: for every shift d ≠ 0 the formal polynomial
  P_d = Σ_j Σ_τ u(j, τ + d)·conj(u(j, τ)) in the commuting symbols y_k, conj(y_k) vanishes on the type, and every
  block index occupies t slots [P: compare coefficients of ρ^d; the values y(ψ) of a type fill a sphere, on which
  a quadratic form is constant only if it is a multiple of the norm].
- Specialisation [P]: an assignment valid for a type is valid for every type obtained by adding relations. So a
  type without a multiplication kills every type that specialises to it.

The classical case is the Williamson type (s, s, s, s): Cooper–Wallis with T-matrices.

### 3.2 Results [M]

| t | types searched | with a multiplication | without |
|---|---|---|---|
| 3 | all 164 letter-relation types (`theory/fable_tmult_types.py`, exhaustive meet in the middle over block classes) | 59, all of dimension ≤ 4 | 105: every type of dimension 5 (34), 6 (18), 7 (4), 8 (1); also 28 of dimension 4 and 20 of dimension 3 |
| 5 | the four minimal classes of dimension 5: (f,s,s,s), (f,f,s,=C), (f,f,s,=A), (f,f,s,=A\*) (`theory/fable_tmult_h.py`, 5.0·10⁷ pairs of block classes each, exact re-check of hash matches) | none | all four; hence every type of dimension ≥ 5 by specialisation |
| 5 | generic type (CP-SAT, `theory/fable_tmult_sat.py`) | – | INFEASIBLE (proof, 37 s), consistent with the line above |

- Controls: the Williamson type has 384 solutions at t = 3 and 18,432 at t = 5; the two codes and CP-SAT agree at
  t = 3 on every type run in more than one of them (`logs/fable_fresh/tmult_types_t3.log`, `tmult_h_t5.log`).
- The types of dimension 4 that do multiply at t = 3: Williamson (s,s,s,s), 384; (A, A, B, C) and (A, A\*, B, C) with
  B, C symmetric and A free, 576 each ("propus-like"); (A, B, A\*, B) and its relatives, 2112 (these force a periodic
  Golay pair, so even v only).
- Product form, any t (`theory/fable_cod4.py`): Z_j = Σ_i T_i ⊗ (±W or ±W\*) with unrestricted T-matrices needs a
  4 × 4 array in the letters whose Gram matrix is |y|²·I as a formal identity. Exhaustive over all such arrays:
  288 for the Williamson type, **0** for every type with a free block that was run (one, two or three symmetric
  blocks, D = C\*, generic). This is the order-4 case of the known fact that square complex orthogonal designs of
  full rate stop at order 2.

### 3.3 Consequence

A type of dimension 4 has 2v sign bits against (v − 1)/2 equations. With log₂ of about 4.0–4.3 bits per equation at
v = 167–223, the count 2^{2v − 4.1·v/2} turns negative near v ≈ 160 [E] (at 167: about 2¹ labelled solutions before
dividing by the v translations; at 191: 2⁻⁸; at 223: 2⁻²²). The open orders 4·t·v have v ≥ 167. Types of
dimension 5 or 6, which the same count allows up to v ≈ 600 and v ≈ 2500, do not multiply. **The route "known GS
core × small T-matrices" is closed for t = 3 and 5, and in product form for every t.**

Side remark, not checked against the literature: the propus-like type gives a GS quadruple of order 3v from one of
order v. It reaches only known orders, and none of its 576 assignments returns a propus-like quadruple
(`theory/fable_propus3.py`) [M], so it does not multiply symmetric Hadamard matrices either.

### 3.4 A relation-restricted variant [M]

The Williamson pattern Z₁ = T₁A + T₂B + T₃C + T₄D, Z₂ = T₂A − T₁B + T₄C − T₃D, Z₃ = T₃A − T₄B − T₁C + T₂D,
Z₄ = T₄A + T₃B − T₂C − T₁D fails for a free A and symmetric B, C, D by exactly
−4·Im(Â)·[B̂·(I₁₂ − I₃₄) + Ĉ·(I₁₃ + I₂₄) + D̂·(I₁₄ − I₂₃)], I_{ab} = Im(T̂_a·conj T̂_b) [M: numerical identity, 200
random points, residual 7·10⁻¹⁴]. So it would work for T-matrices with I₁₂ = I₃₄, I₁₃ = −I₂₄, I₁₄ = I₂₃ at every
character. `theory/fable_tamic.py` enumerates all cyclic T-matrices (every assignment of block and sign to the
positions, the block of position 0 fixed): 24, 480, 4704, 64,800 and 644,160 of them for t = 3, 5, 7, 9, 11, and
**none** has the three relations. Not run: t = 13; non-cyclic groups; assignments outside this pattern for t ≥ 7.

---

## 4. Multiplying two Goethals–Seidel quadruples

### 4.1 Definition

Let V = (V₁..V₄) on G_t and W = (W₁..W₄) on G_v be ±1 quadruples with PSD sums 4t and 4v. A **monomial
composition** is z_j = Σ_i ε_{ji}·v_i^{(e_{ji})}·w_{π_j(i)}^{(f_{ji})} (j = 1..4, x^{(+)} = x, x^{(−)} = conj x, π_j
permutations, every v_i and w_k once per row) with Σ_j |z_j|² = |v|²·|w|² identically on ℂ⁴ × ℂ⁴. The octonion
product in Cayley–Dickson coordinates is one. Put

  Z_j(τ, g) = ½ Σ_i ε_{ji}·V_i(e_{ji}τ)·W_{π_j(i)}(f_{ji}g).

Then Σ_j PSD_{Z_j} = (4t)(4v)/4 = 4tv at every character where both factors are nontrivial, and Z_j is ±1-valued
exactly when three of its four terms agree at every point, that is, when Π_i V_i(e_{ji}τ) and Π_k W_k(f_{jk}g) are
both constant, for every row j [P].

### 4.2 Classification [M]

`theory/fable_octo.py` enumerates all compositions (four Latin squares of row permutations; 1024 sign and
conjugation patterns per row; meet in the middle on 64-bit hashes with exact re-check): **2048 compositions, all on
the Klein-four Latin square**. For each factor the set of conjugation patterns of the four rows is, up to
g ↦ −g, one of

  E = {++++, ++−−, +−+−, +−−+}  or  O = {−+++, +−++, ++−+, +++−},

and all four combinations (E or O for V, E or O for W) occur (`logs/fable_fresh/octo.log`).

### 4.3 Both pattern sets are impossible at odd order [P]

Let X = (X₁..X₄) be ±1 sequences on an abelian group of odd order N > 1 with PSD sum 4N.

- **Parity lemma.** Let T₁, …, T_r be {0, ±1} sequences on the group whose supports are disjoint, symmetric
  (supp = −supp) and cover the group. Then Σ_i PAF_{T_i}(s) is odd for every s ≠ 0, so it is not 0.
  *Proof.* A pair (g, g + s) contributes ±1 when both points lie in one support. The map g ↦ −g − s sends such a
  pair to (−g − s, −g), again inside one support by symmetry, and is an involution on these pairs. Its only fixed
  point is g = −s/2, whose pair {−s/2, s/2} lies in one support. So the number of contributing pairs is odd. ∎
- **E.** Π_k X_k ≡ c and, from the other three patterns, each product X_iX_j is an even function. So
  X(−g) = ±X(g) as vectors. Negate X₄ if needed so that c = −1: at each g one coordinate i(g) differs from the
  other three, and i(−g) = i(g). With T_i(g) = X_{i(g)}(g)·[i(g) = i] one has X_j = 2T_j − Σ_i T_i and
  Σ_j |X̂_j|² = 4 Σ_i |T̂_i|², so the T_i have zero total PAF off 0: T-matrices with symmetric supports. The lemma
  forbids them.
- **O.** The four conditions say X_i(g)X_i(−g) = σ(g) for all i (one symmetric/skew splitting of the group for all
  four blocks) and Π_k X_k = c·σ. On the symmetric part the pattern is "one differs" (sequences T_i as above,
  symmetric); on the skew part the pattern is "all equal or two against two", so X = H₄U there with U_m(g) = ±1 in
  exactly one coordinate m(g) (H₄ the 4 × 4 Sylvester matrix), U skew. T̂ is real, Û is imaginary, the cross term
  has real coefficients and vanishes, and Σ_j |X̂_j|² = 4(Σ_i T̂_i² + Σ_m |Û_m|²). The eight sequences T_i, U_m have
  symmetric disjoint supports covering the group and zero total PAF off 0. The lemma forbids them.

**So no monomial composition of this kind turns two GS quadruples of odd orders t, v > 1 into one of order tv.**
Scope: rows in which every v_i and w_k occurs once; four output blocks.

---

## 5. Two reformulations that gave no search

### 5.1 The norm equation [P]

GS quadruples over ℤ_p correspond to (s_f; β_f), f = 1..4, with s_f odd integers, β_f ∈ ℤ[ζ_p],

  Σ_f s_f² = 4p,  Σ_f β_f·conj(β_f) = p,  β_f ≡ s_f/2 (mod 1 − ζ_p),

through β_f = ½ Σ_i a_f(i) ζ^i. *Proof.* Odd coefficients give Σ_i a(i)ζ^i ≡ Σ_i ζ^i = 0 (mod 2), and
Tr(γ) ≡ −(γ mod (1 − ζ)) (mod p) turns "a(i) is an odd integer" into the congruence. Conversely the two norm
conditions and Parseval give Σ a(i)² = 4p over 4p odd integers, so every a(i) = ±1. ∎ With a multiplier group U the
same holds over the fixed field. So a core is a vector of Hermitian norm p in the standard lattice ℤ[ζ]⁴ in a fixed
residue class; the ±1 condition is automatic. The quaternion algebra behind the form is totally definite, so its
norm-one units are the trivial symmetries only, and bucketing by the residues modulo powers of (1 − ζ) costs as
many guesses as it removes candidates. No algorithm came out of it.

### 5.2 The torus base point [P, E]

On G = F_p\* (cyclic of order p − 1, the group of Theorem 5 for 4p) let X_c(x) = χ(x + c), zero at x = −c.

- [P] Σ_x X_c(x)X_c(xy) = −1 − χ(y) for y ≠ 1 (a quadratic character sum with discriminant c²(1 − y)²). The
  transform is a Jacobi sum: |X̂_c(ψ)|² = p for ψ ∉ {1, χ}, and 1 at ψ = 1, χ.
- [P] Four of them with the zeros filled so that Σ_k η_kχ(c_k) = 0 have total autocorrelation −8 on the squares and
  0 on the non-squares. This is H(4) ⊗ Paley(p + 1) seen through the split torus (border 8).
- Theorem 5 needs −4 everywhere, so the difference is 4χ(y): energy 4p has to move into the single character χ. In
  the ℤ₂ × ℤ_{(p−1)/2} splitting, the difference parts must reach row sums with Σ = p from Σ ≤ 4, which takes at
  least about √p sign changes (25–50 at p = 599), after which 299 autocorrelation equations of about 5 bits each
  must be restored: about 2²⁰⁰ candidate change sets against 2¹⁴⁰⁰ of constraint [E]. Generic solutions sit at
  Hamming distance of order 600 from the base point. Dropped.

---

## 6. EKLO's Theorem 7 with a large index and a tiny subgroup

Setting: G abelian, S ≤ G of index k, Q = G/S, four ±1 sequences with Σ PAF = −4(k − 1) on S ∖ 0 and +4 off S,
n = 4|G| + 4(k − 1) = 4k(|S| + 1) − 4, a border of width b = 4(k − 1). For n = 4p this needs k(|S| + 1) = p + 1.
Earlier censuses (here and Kumar's) took k ≤ 16. The question was whether a *tiny* S with a large k, where the core
lives on a group of order about (p + 1)/4·|S| with many automorphisms, opens new slices.

### 6.1 Cell-energy identity [P]

Call the restriction of a block to a coset of S a *cell*. Summing the character equations over Q̂ gives, for each
nontrivial ψ of S, Σ_cells |f̂_cell(ψ)|² = n.

- S = ℤ₃: a constant cell has f̂(ω) = 0, any other cell 2εω^j. So exactly **4k − 1 cells are non-constant and one is
  constant**; p = 4k − 1.
- Paley-type cells sδ₀ + tχ on F_q (energy q + 1), constants (0) and constants with the zero flipped (4):
  (N_const + N_flip)(q + 1) = 4(N_flip + 1), so q = 7 with one flipped constant (p = 8k − 1), q = 5 with two
  (p = 6k − 1), q = 3 as above, and **no solution with these cell types for q ≥ 9**.

### 6.2 S = ℤ₃: a conference matrix over the sixth roots of unity [P]

Let G = Q × ℤ₃. For a non-constant cell with pattern −ε(1,1,1) + 2εe_j put z = εω^j ∈ μ₆ and y = −ε; for the
constant cell z = 0, y = ±3. Then the core equations are equivalent to

  (i) Σ_b |ẑ_b(ρ)|² = 4k − 1 for every ρ ∈ Q̂ (four blocks on Q over {0} ∪ μ₆, one zero, Hermitian PAF sum 0),
  (ii) Σ_b |ŷ_b(ρ)|² = 4k − 4 for ρ ≠ 1 and Σ_b ŷ_b(1)² = n, where y = −z³ off the zero.

*Proof.* The characters of G are (ρ, 1), (ρ, ω), (ρ, ω̄). The PAF pattern is equivalent to PSD sum n at characters
nontrivial on S, b at nontrivial characters trivial on S, and n at the trivial one. X̂_b(ρ, ω) = 2ẑ_b(ρ) and
X̂_b(ρ, 1) = ŷ_b(ρ). ∎ More generally an automorphism of order 3 of an H(4p) with f fixed points gives a Butson
weighing matrix over μ₆ of order (4p − f)/3 and weight p.

**Paley's conference matrices do not fit [P].** In the general form (an automorphism of order 3 with p − 3 fixed
points, no block structure assumed) Z is a (p + 1) × (p + 1) matrix over {0} ∪ μ₆ with one zero in each row and
Z Z\* = pI, Y = −Z³ entrywise with ±3 at the zeros, and the fixed part needs 3FFᵀ = 4pI − YYᵀ for a ±1 matrix F
with p − 3 columns. So YYᵀ must have the eigenvalue 4p with multiplicity at least 4.
- Candidates for Z: Paley's skew conference matrix C₂ = χ₂(det(x, y)) on PG(1, p), and for p ≡ 1 (mod 6) the sextic
  one C₆ = χ₆(det(x, y)), which also satisfies C₆C₆\* = pI; and every monomial equivalent D₁ZD₂ over μ₆.
- Cubing kills μ₃: (D₁C₆D₂)³ = Δ₁C₂Δ₂ entrywise with Δ = ±1. So in every case Y = Δ₁(3D − C₂)Δ₂ with D a ±1
  diagonal, and YYᵀ is similar to (p + 9)I − 3[C₂, D]. The commutator has norm at most 2√p, so the eigenvalue 4p
  would need p − 3 ≤ 2√p. **For p ≥ 11 no matrix monomially equivalent to a Paley conference matrix, quadratic or
  sextic, can be the Z of such an H(4p).** With D = I the miss is uniform: YYᵀ = (p + 9)I, a PAF sum 0 where (ii)
  needs +12.

### 6.3 The border [P for the equivalence, E for the conclusion]

Let M be the GS array of the compressions y (order 4k) and W ⊂ ℝ^{4k} the vectors that sum to 0 on each of the four
blocks. For borders whose E has rows in W (every known one):

  a border (K, E, F) exists ⟺ there are ±1 matrices E, Z of size b × 4k with rows in W, EEᵀ = ZZᵀ = 4kI, and every
  entry of E Mᵀ Zᵀ equal to ±4k.

Then K = E Mᵀ Zᵀ / 4k and F = −Zᵀ. *Proof.* The off-diagonal block of H Hᵀ is the lift of K Fᵀ + E Mᵀ, so
F = −M Eᵀ K / b. With M Mᵀ = I₄ ⊗ (bI + 4|S|J), M(I₄ ⊗ J)Mᵀ = n(I₄ ⊗ J) and EᵀE = 4kI − 4(I₄ ⊗ J) one gets
F Fᵀ = 4kI − 4(I₄ ⊗ J), the required Gram matrix, automatically; Z = −Fᵀ has rows in W and Z Zᵀ = 4kI; and
E Mᵀ Zᵀ = 4kK. The converse is the same computation read backwards. ∎ This is Kumar's formulation with M in place
of his Γ ⊗ I. On W the matrix M/√b is orthogonal, and the condition says that the orthonormal bases M·E and Z of W
are related by a real Hadamard matrix.

- For a spike compression (constant plus a multiple of δ₀, Kumar's case, and the k = 4 normal form behind
  Theorem B) M acts on W as a 4 × 4 matrix tensor the identity, and E, Z can be chosen once for all cores.
- For S = ℤ₃ the compressions are ±1-valued with a single ±3. A spike needs Σ d² = 4(k − 1) with |d| ≤ |S| + 1, so
  k ≤ 17 at |S| = 3, 37 at 5, 65 at 7, and a group acting transitively on Q ∖ 0.
- For a compression without that structure the condition is (4k − 4)² exact equalities, each met with probability
  of order 1/k by an integer of standard deviation 4k [E]. At k = 6 already a random pair has about 10⁻¹¹ expected
  admissible columns. Nothing suggests a border for a generic compression at k ≥ 6.

**So large k with |S| ≤ 7 is closed as a route**: the core would be cheap to state, but its border is a
Hadamard-type object of order 2p − 2 that nobody knows how to build.

### 6.4 Which borders are possible at all, and why none serves the skew census [P]

Let G be abelian of order ks, S ≤ G of index k and order s > 1, Q = G/S, and suppose
H = [[K, E′], [F′, C]] is a Hadamard matrix of order n = 4ks + b with C the Goethals–Seidel array of four ±1
sequences on G, K of order b > 0, and E′, F′ constant on the cosets of S inside each of the four blocks (E, F their
quotients, of sizes b × 4k and 4k × b). Write P = Σ PAF.

1. P = −b on S ∖ 0, and P is constant on every other coset: P = β(q), β(−q) = β(q).
   *Proof.* C Cᵀ = I₄ ⊗ Dev(P), so F′F′ᵀ = I₄ ⊗ (nI − Dev(P)) must depend on cosets only. ∎
2. F Fᵀ = EᵀE = I₄ ⊗ B with B = bI − Circ_Q(β); its eigenvalues are b − β̂(ρ), ρ ∈ Q̂. (CᵀC = I₄ ⊗ Dev(P) as well,
   because Σ XᵢXᵢᵀ commutes with every block and with the inversion matrix.)
3. Let T = {ρ : β̂(ρ) ≠ b}. Then 4|T| = rank E ≤ b, and T is closed under complex conjugation. If n > b² then
   b = 4|T|. *Proof.* KᵀK = nI − sFᵀF has the eigenvalue n on the kernel of F, and ‖K‖² ≤ b². ∎
4. If k is odd, the trivial character lies in T. *Proof.* Otherwise B·1 = 0, so every row of E sums to 0 over a
   block of odd length k. ∎

These are Kumar's (Q, T) character borders, recovered from the shape alone: T = Q̂ ∖ {1} with β = 4 is EKLO's
Theorem 7 (so k must be even there), T = Q̂ with β = 0 is the all-character family.

**Consequence for skew matrices.** A skew-Hadamard H of this shape has a skew-type first block, so |G| is odd, so k
is odd. Then 1 ∈ T, the other characters of T pair off under conjugation, and |T| is odd. With n > b² this gives
n = 4(|G| + |T|) ≡ 0 (mod 8). **No bordered Goethals–Seidel array with a coset-constant border of width
0 < b < √n is a skew-Hadamard matrix of an order 4v with v odd.** The census values can only come from unbordered
arrays over a group of order exactly v, or from arrays of another kind (Kumar's negacyclic arrays give 8v for the
same reason).

### 6.5 Uniform cyclotomy in the fibres [P]

In EKLO's family over Q × F_{q²} (S = F_{q²}), let every cell be constant on the q + 1 punctured lines through 0:
f = f₀δ₀ + Σ_i c_i·1_{L_i}, so f̂(ψ) = (f₀ − σ) + q·c_{i(ψ)}, σ = Σ c_i. Summing over Q̂ for one line i gives
Σ d² + 2q Σ d·c_i = 4k − 4 with d = f₀ − σ odd. Writing Σ d·c_i = −u (then u ≡ 2 mod 4) and summing over the lines
gives Σ f₀σ = 4 − u(q − 1) and Σ σ² = 4 + 2u; since σ is even, |Σ f₀σ| ≤ ½ Σ σ², hence u(q − 2) ≤ 6 and **q ≤ 5**.
The cores EKLO found over F₁₉² and F₂₃² necessarily use finer classes.

---

## 7. Partial character borders: a family of bordered arrays that no census has used

§6.4 shows that a coset-constant border forces a core with

  Σ_k |x̂_k(χ)|² = n at every character χ of G outside T,  = 4|T| at χ ∈ T,  n = 4|G| + 4|T|,

for a set T of characters of Q = G/S. Equivalently Σ PAF(g) = −4 Σ_{ρ∈T} ρ(g) for g ≠ 0. The known cases are
T = {1} (Theorem 5), T = Q̂ ∖ {1} (EKLO's Theorem 7) and T = Q̂ (Kumar's all-character family); Kumar's censuses also
take T defined by character orders.

**Observation.** When Q = ℤ₂^a every character of Q is itself a ±1 vector, so **every subset T ⊂ Q̂ is a candidate
shape**, and:
- [P] the compressed core acts on the border space as a direct sum of 4 × 4 integer blocks: M(e_a ⊗ ρ) = Γ_ρ(e_a) ⊗ ρ
  for ρ ∈ T, with Γ_ρ the Goethals–Seidel matrix of the four Walsh values ŷ_k(ρ) of the compressions (the characters
  are real and inversion is trivial on Q), Γ_ρΓ_ρᵀ = 4|T|·I. So the existence of a border depends only on T and on
  these |T| integer quadruples, not on S or on the rest of the core. Theorem B is the case T = Q̂ ∖ {1}, Q = ℤ₂²,
  where the normal form leaves one quadruple type.
- Translating T by a character (twisting every block by it) and applying GL(a, 2) give equivalent shapes, so shapes
  are affine classes of subsets of F₂^a. A point is Theorem 5; a flat is the all-character family on a quotient; a
  flat minus a point is EKLO's family. **Every other class is new.** There are, for instance, 2 classes of 4-subsets
  of F₂³ (one new), 3 classes of 7-subsets of F₂⁴ (two new) and 2 classes of 11-subsets of F₂⁴ (both new) [M:
  `theory/fable_tclasses.py`, which also lists the open primes per residue class].

**Verified instances [M]** (cores by simulated annealing on Q × ℤ₃, borders by CP-SAT, the full matrix assembled and
H Hᵀ = nI checked; `theory/fable_border2.py`, `theory/fable_border4.py`, logs `logs/fable_fresh/border*`):

| Q | T | shape | n | result |
|---|---|---|---|---|
| ℤ₂² | {1} | point (Theorem 5) | 52 | border found, H(52) verified (control) |
| ℤ₂² | three nontrivial | flat minus point (EKLO k = 4) | 60 | border found, H(60) verified (control; several cores) |
| ℤ₂³ | three independent | flat minus point after a twist (k = 4 with even \|K\|) | 108 | border found for 2 of 2 random cores, H(108) verified |
| ℤ₂³ | {1, χ₁, χ₂, χ₄}: four affinely independent | **new** | 112 | single-character border found (323 s), **H(112) verified** |
| ℤ₂³ | five characters | new | 116 | general border: CP-SAT UNKNOWN at 600 s (not a zero) |
| ℤ₂⁴ | seven characters not in a 3-flat | new, the shape of 599 | 220 | <!-- T7 -->running at the time of writing<!-- /T7 --> |

The orders 112 and 116 are of course known; the point is that a border exists for a character set outside the three
classical shapes, for a core with unstructured compressions.

**Where the new shapes sit.** For n = 4p one needs p = 2^a·s + |T|, so |T| ≡ p (mod 2^a) is determined by a. For
a ≤ 3 the shape is classical. New shapes start at a = 4: |T| = 7 for p ≡ 7 (mod 16) (599, 1031, 1223, 1319, 1367,
1543, 1559, 1783, 1831, 1847, 1879, …) and |T| = 11 for p ≡ 11 (mod 16) (907, 971, 1019, 1259, 1291, 1451, 1499,
1579, 1787, 1979, 2027, …); at a = 5, |T| = 15 for p ≡ 15 (mod 32) (719, 751, 1039, 1103, 1327, 1423, 1487, 1583, 1871, …),
and so on. Every open prime gets new shapes on groups ℤ₂^a × S that the censuses ran only with T = Q̂ ∖ {1}.

**Two first gates, both negative [P].** The Q-level rigidity of the classical shapes is not lifted automatically:
- 599 = 16·37 + 7 over ℤ₂⁴ × ℤ₃₇ with the multiplier group of order 9 on ℤ₃₇ (5 orbits per coset, 80 per block).
  Cell sums are ≡ ±1 (mod 18); Parseval at the Q-level gives Σ y² = 1360 over 64 cells, so exactly two cells of
  size 17, two of size 19 and sixty of size 1. A block containing one big cell B would need its sixteen-point Walsh
  transform to be at most 5 in absolute value at the seven characters of T, so its ±1 part needs |Walsh| ≥ 12 there:
  7·12² > 16·15. Two big cells in one block cancel only where ρ(q₁ + q₂) is constant, which puts T inside an affine
  hyperplane, and then the eighth character of that hyperplane cannot reach 2396. All four big cells in one block:
  their signed Walsh sum ±17 ± 17 ± 19 ± 19 either takes the value ±72 at some character (72² > 2396) or is never
  below 34, so it cannot be small on T. **Dead for every T.**
- 751 = 32·23 + 15 over ℤ₂⁵ × F₂₃ with Paley-type cells. The cell-energy identity forces 125 Paley cells, 2 constant
  cells and 1 constant with a flipped zero; the S-nontrivial equations become Σ_k (ŝ_k(ρ)² + 23 t̂_k(ρ)²) = 3004 for
  the 32 characters of Q (no cross term, because the Walsh values are real). But the three special cells have sums
  ±23, ±23, ±21 and cannot cancel at 15 characters by the same Parseval count. **Dead for every T.**

So the family is real, it is unexplored, and its cheapest instances at the two smallest orders tried die at the
same kind of gate as the classical shapes. A census (shapes × groups × multiplier groups with the exact Q-level
gate, then a border for the surviving Walsh quadruples) is the natural next computation; it was not done here.

## 8. The production run: skew 361 over ℤ₁₉ × ℤ₁₉ [E; a finished run without a hit excludes nothing]

`report/skew_census_answer.md` left one calibrated slice unfinished: v = 361, H = ⟨diag(4, 7)⟩, pattern
(1; 17, 23, 25), searched to 10,000 of 33,600 ranked split values. Its lists were not kept, so the run was restarted
with new Metropolis lists (seed 2, target 3·10⁷ per block) and the same pipeline (`theory/skew_kt.py`, Kumar's
`kt4.c`); launcher `scripts/fable_run361.sh`, log `logs/fable_fresh/p361_s2.log`. Expectation: 0.5 to 7 twin
classes in the whole lists, depending on whether Kumar's measured Metropolis yield (0.07) or the Gaussian count
applies.

<!-- RUN361 -->
Status at the time of writing: lists being generated; container restarts every few hours cost the unfinished block
each time. Result: see the final lines of `logs/fable_fresh/NOTES.md`.
<!-- /RUN361 -->

---

## 9. What I would do next

1. **With about 10³ core-hours:** Kumar's Spence plug-in for 4792 = 8·599 over ℤ₇ × ℤ₁₉ × F₉ (his `search/eightp/gs`,
   calibrated on 5912 within a factor of two), then the skew slices he ranks at 10³ (661, 883, 1483). These are the
   only live targets with a measured rate.
2. **A census of partial character borders** (§7): affine classes of T in F₂^a for a = 4, 5, 6, every group
   ℤ₂^a × S (and the non-split 2-parts) of order p − |T|, every multiplier group, the exact Q-level gate, and a
   border library per Walsh-quadruple type. This is the one structural family found here that is both new and
   cheap to screen with the repository's existing gates.
3. **Borders.** §6.3 turns "does a border exist for this compression" into a concrete question about pairs of ±1
   bases. A theorem for a class of non-spike compressions at k = 6, 8 or 10 would open the mid-range indices
   (k = 10–60 with |S| = 11–59 at 599 and 719), which no census has touched because no border is known there.
4. **Formal multiplications at t ≥ 7** for the dimension-5 types outside the Williamson pattern. The count does not
   exclude a sporadic identity; t = 7 needs a three-level join (7·10⁵ block classes).
5. **Not worth repeating:** product theorems of the two kinds closed in §3 and §4; Theorem 7 with |S| ≤ 7.

## 10. Reproduce

```
python3 theory/fable_tmult_types.py 3            # all 164 types at t = 3
python3 theory/fable_tmult_h.py 5 "f s s s" "f f s =2" "f f s =0" "f f s *0"
python3 theory/fable_tmult_sat.py 5 generic 600  # CP-SAT, INFEASIBLE
python3 theory/fable_cod4.py all_sym AB_sym A_sym D_eq_Cstar generic
python3 theory/fable_octo.py                     # 2048 compositions, pattern sets E / O
python3 theory/fable_tamic.py 3 5 7 9 11
python3 theory/fable_border4.py 3 3 0,1,2,4 1 380 1   # the new border shape, H(112)
python3 theory/fable_border2.py 3 3 1,2,4 2 300 1       # twisted k4 with even |K|, H(108)
scripts/fable_run361.sh 4                        # needs /home/user/skew_work/bin (see logs/fable_fresh/NOTES.md)
```
