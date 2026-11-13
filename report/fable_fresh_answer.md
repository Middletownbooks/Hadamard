# A fresh look at the open Hadamard and skew-Hadamard orders: what was tried, what is now closed

Date: 2026-11-09 to 2026-11-13 (second phase on 2026-11-13: §6.4 item 5, §7.2, §8a and the σ-invariant result in §7.1). Branch `claude/focused-planck-j1645f`. Answers `report/prompt_fable_fresh.md`.
Scripts: `theory/fable_*.py`. Logs and working notes: `logs/fable_fresh/` (`NOTES.md` is the session's state file).

**Tags.** [P] proved here (proof in the text); [M] exact computation (script named); [E] estimate. A cap or a timeout
is not a zero.

**Credit.** Egan, Kotsireas, London and Ó Catháin (EKLO) for the survey and Theorem 7; Alpöge, Voinov and
Reynolds-Haertle for the twelve matrices of August 2026; K. Kumar for his constructions, his bordered arrays with a character set T and
their border problem (`search/NOTES.md` §18 and §102 of his repository) and `kt4.c`; the OpenAI model for Theorem B and the compression normal
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
4. **Bordered arrays cannot reach the skew census** (§6.4) [P]. A coset-constant border (of width below √n) of a
   Goethals–Seidel core is necessarily one of Kumar's character borders (his `search/NOTES.md` §18): it is carried
   by a set T of characters of the quotient and has width 4|T|. Kumar notes that for a quotient of odd order the
   trivial character must lie in T; it follows that |T| is odd whenever the group has odd order, as a skew-type
   block requires, so the order is ≡ 0 (mod 8). The 630 open values 4v, v odd, can only come from unbordered arrays
   over a group of order v or from a different kind of array.
5. **An unused part of Kumar's character-border family** (§7). His construction allows any Galois-stable set T of
   characters of the quotient Q; the censuses (his and this repository's) take T = {1} (Theorem 5), all nontrivial
   characters (Theorem 7), all characters, or sets defined by character orders. Over Q = ℤ₂^a every subset is
   admissible, the shapes are the affine classes of subsets of F₂^a, and the border depends only on |T| integer
   quadruples [P]. One instance outside the three classical shapes was built and verified: four affinely
   independent characters of ℤ₂³ over ℤ₂³ × ℤ₃, a Hadamard matrix of the (known) order 112 [M]. For 4p such shapes
   start at ℤ₂⁴ (|T| = 7 or 11) and exist at every open prime. With a multiplier acting on S alone, the first
   instances (599, 751) die at the Q-level for every T [P]. **With a multiplier of order 3 acting on Q as well, the
   Q-level of 4·599 = 2396 is alive** (§7.1): G = ℤ₂⁴ × ℤ₃₇, |H| = 9, 78 orbits per block, ten non-classical T,
   2688 exact Q-level solutions [M]. No border was found for them (a cyclic single-character ansatz is exactly
   empty [M]; general borders are undecided), and no core was sought.
6. **The idle cores reran the one calibrated slice this repository had left unfinished** (skew 361, split torus
   over ℤ₁₉²) with fresh lists: top 2000 split values joined, no solution (§8).
7. **Reading.** Among the classical structures, everything cheap enough for 4 cores is now either searched or
   closed by a theorem of the kind above; Kumar's own ranking puts the nearest live targets at about 10³ core-hours
   (4792 = 8·599 by a Spence plug-in; skew 661, 883, 1483). The character sets of item 5 are the one place where
   unscreened cheap structures remain. §9 says what I would run next.
8. **Second phase (2026-11-13), after the first seven items were written.**
   - **Every border of a Goethals–Seidel core over an abelian group is a character border** (§6.4 item 5) [P]:
     the coset-constancy assumed in item 4 is automatic.
   - **A symmetry that every such border problem over ℤ₂^a has** (§7.2) [P]: at a real character the Goethals–Seidel
     array is a quaternion multiplication on the four blocks, so the quaternion group Q₈ acts on the other side by
     signed block permutations. Searching for borders invariant under Q₈, or under one of its cyclic subgroups of
     order 4, is an exact clique problem that takes seconds where a generic CP-SAT model took minutes or gave up.
     Verified instances: H(112), H(108), and H(116) from a five-character set that CP-SAT had left undecided [M].
   - **The 2396 slice of item 5 has no border in any of the symmetric classes tested** (§7.1, §7.2) [M, exhaustive]:
     <!-- BL8 -->none invariant under Q₈ or under any of its three cyclic subgroups of order 4; none invariant under the
     slice's own symmetry σ of order 3 (CP-SAT had been undecided); and none invariant under any other subgroup of
     order 3 of the problem's symmetry group (five classes in all). So a border, if there is one, has no symmetry
     of order 3. The Q-level solutions with a non-constant profile (representative: solution 2160) have no border
     invariant under Q₈, ⟨i⟩ or ⟨j⟩ either; ⟨k⟩ and σ are running.<!-- /BL8 --> Not excluded: a border whose only symmetries are involutions that move the
     quotient Q, or one with no symmetry at all.
   - **The tier-0 slices of the earlier screen** (§8a) [M, P]: I reran them before noticing that
     `report/h4892_answer.md` had already decided most; the reruns reproduce its counts, including the 1360 stage-1
     solutions of the 7572 slice and its emptiness. New here: the slice of 4156 = 4·1039 (undecided there) is empty
     at stage 1, so only 4748 remains undecided in that tier; and a one-line a-priori gate (a multiplier group that
     contains all of ℤ_ℓ\* on a prime factor ℓ ≥ 5 leaves no core) that kills 18 of the 81 ranked slices, 7572 among
     them, before any enumeration.
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
| 3 | Williamson pattern with T-matrices that satisfy three "half-amicability" relations (one free block, three symmetric) | exhaustive over cyclic T-matrices, t ≤ 13 | none exist [M], §3.4 |
| 4 | Theorem 7 with large k, S = ℤ₃, F₅ or F₇ (p = 4k − 1, 6k − 1, 8k − 1) | cell-energy identity, sixth-root reformulation, border analysis | closed [P, E], §6 |
| 4a | Which borders can exist at all (shape analysis of coset-constant borders) | Gram matrices of the border | Kumar's character borders only; none for skew 4v, v odd [P], §6.4 |
| 4b | Character borders with a set T of characters of ℤ₂^a outside the three classical shapes | cores by annealing, borders by CP-SAT, full matrix checked | **works**: H(112) from a non-classical T [M]; two open-order gates negative [P]; no census, §7 |
| 5 | The torus base point: G = F_p\* (the Theorem 5 group of order p − 1), blocks χ(x + c) | exact autocorrelation, energy count | no nearby solution [P, E], §5.2 |
| 6 | Norm-equation form over the fixed field of the multiplier group | reformulation only | no algorithm gained [P], §5.1 |
| 7 | Uniform cyclotomy (punctured lines of F_q²) in the fibres of Theorem 7 | energy identity | impossible for q ≥ 7 [P], §6.5 |
| 8 | Ito / dicyclic matrices (negaperiodic Golay pairs of length 2p) | sign count | negative for p > about 160 [E]; dropped |
| 9 | Legendre pairs over a group of order 2p − 1 | sign count; Kumar screened the four-block form | dropped |
| 10 | 10012 = 4·2503 over ℤ₂² × F₆₂₅ (the EKLO family at q = 25) | Kumar's notes | already examined by Kumar (F₂₅ × F₂₅ exactly empty; "2 % per run" otherwise); dropped |
| 11 | 4792 = 8·599 by Kumar's Spence plug-in ℤ₇ × ℤ₁₉ × F₉ | his calibrated cost | 1.5·10³ core-hours against 288 available; dropped |
| 12 | Skew T-structures for Paley blocks beyond t = 3 | Kumar's exhaustive t = 9, 19; count at t = 33 | about 10⁻³ expected [E]; dropped |
| 13 | Unstructured small skew values (89, 101, 107, 149) | size | 2³¹¹ sign vectors at 89 with no multiplier; dropped |
| 14 | (second phase, 2026-11-13) Borders with the block symmetry Q₈ that every Goethals–Seidel core over ℤ₂^a has ("quaternionic borders"), and the weaker cyclic ones | exact clique search, calibrated on four solved cases | a fast border finder (H(112), H(116), H(108) in under a second); for the 2396 slice: **no Q₈-invariant and no σ-invariant border** [M], §7.1, §7.2 |
| 15 | Exact stage 1 of the tier-0 slices the screen lists as "not run" (most were already decided in `report/h4892_answer.md`, which I had overlooked), then tiers 1 and 2 | the screen's own method B with planted controls | reruns agree with the earlier report; **new: 4156 empty**; an a-priori energy gate that kills 18 of the 81 tiered slices [P, M], §8a |
| 16 | Cores over non-abelian groups with class-function blocks; commutative association schemes in place of multiplier orbits | cell energy on the cosets of the commutator subgroup | Frobenius and extraspecial groups have constant cells off the centre or off G′, which the energy identity forbids [P]; no scheme on 599 − t points was found to try; dropped |
| 17 | Time–frequency (Heisenberg) symmetries x(g + l) = ±η(g)x(g) | transform computed | they are Kronecker products (order ≡ 0 mod 16) when all blocks share the symmetry, and cost a factor 2 in equations otherwise [P]; dropped |
| 18 | The Alpöge half-shift template read as a ℤ₂ split with arc supports (two aperiodic complementary quadruples of length about v/4) | count | negative beyond v ≈ 450, in line with `report/alpoge_methods_forensics.md` [E]; dropped |

Ideas 8–13 and 16–18 were dropped on the stated one-line reason and are not discussed further.

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
**none** has the three relations. t = 13: all 64 chunks, 5,710,848 T-matrices, **none** with the relations (`theory/fable_tamic13.py`,
`logs/fable_fresh/tamic13.jsonl`). Not run: non-cyclic groups; assignments outside this
pattern for t ≥ 7.

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
E Mᵀ Zᵀ = 4kK. The converse is the same computation read backwards. ∎ This is Kumar's formulation (his `search/NOTES.md` §18 and §102) with a general M in
place of his Γ ⊗ I. On W the matrix M/√b is orthogonal, and the condition says that the orthonormal bases M·E and Z of W
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

5. The hypothesis that E′ and F′ are constant on cosets is automatic [P, added 2026-11-13]. For any Hadamard matrix
   [[K, E′], [F′, C]] with C the Goethals–Seidel array of four ±1 sequences on G, F′F′ᵀ = nI − C Cᵀ = I₄ ⊗ Dev(nδ₀ − P)
   is G-invariant, so the column space of each block of F′ is the sum of the character spaces of
   T′ = {χ ∈ Ĝ : Σ_k |x̂_k(χ)|² ≠ n}. Every function in that space is constant on the cosets of S = ⟨T′⟩^⊥, so the
   columns of F′ are; the same holds for E′ through CᵀC. S may be trivial (then Q = G), and items 1–4 apply with
   that S. **So every border of a Goethals–Seidel core over an abelian group is a character border.**

These are exactly Kumar's bordered arrays with a character set T (his `search/NOTES.md` §18, which states the
construction, the orthogonal-matrix form of the border problem, and item 4). What the computation above adds is the
converse: a coset-constant border cannot have any other shape. T = Q̂ ∖ {1} with β = 4 is EKLO's Theorem 7 (so k must
be even there), T = Q̂ with β = 0 is the all-character family.

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

## 7. Character borders over ℤ₂^a with a general character set: an unused part of Kumar's family

§6.4 shows that a coset-constant border forces a core with

  Σ_k |x̂_k(χ)|² = n at every character χ of G outside T,  = 4|T| at χ ∈ T,  n = 4|G| + 4|T|,

for a set T of characters of Q = G/S. Equivalently Σ PAF(g) = −4 Σ_{ρ∈T} ρ(g) for g ≠ 0. This is Kumar's
construction (`search/NOTES.md` §18 of his repository), stated there for any Galois-stable T. The cases that have
been used are T = {1} (Theorem 5), T = Q̂ ∖ {1} (EKLO's Theorem 7), T = Q̂ (the all-character family), products of
these, and, in his censuses, T defined by character orders; on ℤ₂² he notes that {1, χ₁, χ₂} is the nontrivial set
after a twist. **Correction (2026-11-13):** his §102 ("Alternatives that fail", partial character sets) does record
one more case, t = 5 over ℤ₂³ or F₁₆, with the remark that borders exist and an estimate E ≤ 10⁻¹¹ for the target
he had in mind, and partial Teichmüller sets over ℤ₄ʳ (599 over ℤ₄² × ℤ₃₇, t = 7: no border in any symmetric class he
tested). So partial sets over ℤ₂³ are his as well. What I did not find there is a classification of the shapes or
any use of them over ℤ₂^a with a ≥ 4; the wider literature was not checked.

**Observation.** When Q = ℤ₂^a every character of Q is itself a ±1 vector, so **every subset T ⊂ Q̂ is a candidate
shape**, and:
- [P] the compressed core acts on the border space as a direct sum of 4 × 4 integer blocks: M(e_a ⊗ ρ) = Γ_ρ(e_a) ⊗ ρ
  for ρ ∈ T, with Γ_ρ the Goethals–Seidel matrix of the four Walsh values ŷ_k(ρ) of the compressions (the characters
  are real and inversion is trivial on Q), Γ_ρΓ_ρᵀ = 4|T|·I. So the existence of a border depends only on T and on
  these |T| integer quadruples, not on S or on the rest of the core. Theorem B is the case T = Q̂ ∖ {1}, Q = ℤ₂²,
  where the normal form leaves one quadruple type.
- Translating T by a character (twisting every block by it) and applying GL(a, 2) give equivalent shapes, so shapes
  are affine classes of subsets of F₂^a. A point is Theorem 5; a flat is the all-character family on a quotient; a
  flat minus a point is EKLO's family. **Every other class is unused.** There are, for instance, 2 classes of 4-subsets
  of F₂³ (one unused), 3 classes of 7-subsets of F₂⁴ (two unused) and 2 classes of 11-subsets of F₂⁴ (both unused) [M:
  `theory/fable_tclasses.py`, which also lists the open primes per residue class].

**Verified instances [M]** (cores by simulated annealing on Q × ℤ₃, borders by CP-SAT, the full matrix assembled and
H Hᵀ = nI checked; `theory/fable_border2.py`, `theory/fable_border4.py`, logs `logs/fable_fresh/border*`):

| Q | T | shape | n | result |
|---|---|---|---|---|
| ℤ₂² | {1} | point (Theorem 5) | 52 | border found, H(52) verified (control) |
| ℤ₂² | three nontrivial | flat minus point (EKLO k = 4) | 60 | border found, H(60) verified (control; several cores) |
| ℤ₂³ | three independent | flat minus point after a twist (k = 4 with even \|K\|) | 108 | border found for 2 of 2 random cores, H(108) verified |
| ℤ₂³ | {1, χ₁, χ₂, χ₄}: four affinely independent | **not classical** | 112 | single-character border found (323 s), **H(112) verified** |
| ℤ₂³ | five characters {1, χ₁, χ₂, χ₃, χ₄} | not classical (Kumar notes that borders exist for t = 5 over ℤ₂³) | 116 | CP-SAT UNKNOWN at 600 s; **border found in under 1 s by the quaternionic search of §7.2, compressed Gram identity verified** |
| ℤ₂⁴ | seven characters not in a 3-flat | not classical, the shape of 599 | 220 | <!-- T7 -->running at the time of writing<!-- /T7 --> |

The orders 112 and 116 are of course known; the point is that a border exists for a character set outside the three
classical shapes, for a core with unstructured compressions.

**Where the unused shapes sit.** For n = 4p one needs p = 2^a·s + |T|, so |T| ≡ p (mod 2^a) is determined by a. For
a ≤ 3 the shape is classical. Unused shapes start at a = 4: |T| = 7 for p ≡ 7 (mod 16) (599, 1031, 1223, 1319, 1367,
1543, 1559, 1783, 1831, 1847, 1879, …) and |T| = 11 for p ≡ 11 (mod 16) (907, 971, 1019, 1259, 1291, 1451, 1499,
1579, 1787, 1979, 2027, …); at a = 5, |T| = 15 for p ≡ 15 (mod 32) (719, 751, 1039, 1103, 1327, 1423, 1487, 1583, 1871, …),
and so on. Every open prime gets such shapes on groups ℤ₂^a × S that the censuses ran only with T = Q̂ ∖ {1}.

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

**A first table [M].** `theory/fable_tcensus.py` lists, for every open prime, a = 4, 5, 6, cyclic S and every cyclic
multiplier group on S with 24 to 130 orbits per block, whether the Parseval sum Σ y² can be written with the available
cell sums at all: 188 slices, 142 dead at this level, 46 pass (`logs/fable_fresh/tcensus.log`).
- Of the 46, eight have a multiplier group of odd order, the case with eight sign bits per equation. Five of them
  admit only a classical T (|T| = 3 at 739, 1619 and 1747; |T| = 15 over ℤ₂⁴ at 751). The other three are 599 with
  the group of order 9 and 751 over ℤ₂⁵ × ℤ₂₃ (both dead above) and 2287 over ℤ₂⁵ × ℤ₇₁, which has a single cell of
  size 69 or 71 and is dead by the same Parseval count. **So no slice with cyclic S, a multiplier group of odd
  order acting on S alone, at most 130 orbits and a non-classical T survives.**
- The other 38 have a group of even order (it contains −1 on ℤ_s when s is an odd prime, which halves the sign
  bits per equation) or an even s. A few have a real choice of cell sums, for example 599 with the group of order
  6 (112 orbits, cell sums 1, 11, 13, 23, 25), 1499 and 1579 with |T| = 11, 1847 and 1223 with |T| = 7. None was
  decided at the Q-level (`theory/fable_qlevel.py`, a CP-SAT gate for all T at once, is too weak: UNKNOWN at 300 s
  even on the 599 slice that is dead by hand).
- Not covered: non-cyclic S, multipliers that act on Q as well, non-split 2-parts (Galois rings), products of
  multiplier groups on several prime factors.

So these shapes exist, they are unexplored, and their cheapest instances at the two smallest orders tried die at the
same kind of gate as the classical shapes. A census (shapes × groups × multiplier groups with the exact Q-level
gate, then a border for the surviving Walsh quadruples) is the natural next computation; it was not done here.

### 7.1 A live Q-level at 2396 with a non-classical character set [M]

`theory/fable_tcensus2.py` repeats the gate with multiplier groups ⟨σ⟩ × U, σ an element of odd order of GL(4, 2)
acting on Q = ℤ₂⁴ (cells and characters), U a cyclic group of units of ℤ_s, T a union of σ-orbits of characters, not
classical. It is exact (block enumeration by a meet in the middle with the bound |ŷ| ≤ √(4|T|) on T, then a 2 + 2
join on energy vectors) and restart-safe per prime. Done for 599, 907, 971 and 1019 (`logs/fable_fresh/tcensus2.jsonl`,
331 records): 66 dead, 65 alive, 200 not decided (caps; not zeros).

- **The slice.** 4·599 = 2396 = 4(16·37 + 7). G = ℤ₂⁴ × ℤ₃₇, H = ⟨σ⟩ × ⟨10⟩ with σ of order 3 fixing only 0 of
  ℤ₂⁴ and 10 of order 3 in ℤ₃₇\* (|H| = 9, odd; 6 σ-orbits of Q times 13 U-orbits of ℤ₃₇ = 78 orbits per block).
  T = the union of two σ-invariant 2-flats of characters through the trivial one, e.g. {0, 1, 2, 3, 4, 8, 12} (ten
  such T, one affine class, not classical).
- **Q-level: 2688 solutions** for T = {0, 1, 2, 3, 4, 8, 12} (`theory/fable_q599.py`, every solution re-checked
  from scratch; `logs/fable_fresh/q599_solutions.json`). Cell sums are small (1, 5, 7, 11, 13, 17). At every
  character of T the four Walsh values are (±2, ±2, ±2, ±4), and in 2592 of the 2688 each block's Walsh value is
  constant on T, so the compressed core acts on the border space as Γ ⊗ I₇ with Γ a 4 × 4 matrix of entries ±2, ±4,
  ΓΓᵀ = 28I. This is the same Γ as for Kumar's spike compressions at index 8 (his d = (4, 2, 2, 2)).
- **Border: not found.** The admissible block functions (±1 functions on Q with Walsh support in T) are 30: the
  functions of the two coordinates dual to one 2-flat or to the other. Kumar's index-8 border for 3356 cannot be
  transplanted: its block rows are balanced but mostly not single characters, and its four-point Walsh supports
  meet in two points where ours meet in one. Single-character borders with the cyclic symmetry l ↦ l + 1 on the 7
  labels: **exactly none** (`theory/fable_border9.py`: 2744 orbit representatives per side, 576 compatible pairs,
  all 13,034 orthogonal E-quadruples tried). Tensor bases H₄ ⊗ (7 block functions): 0 admissible Z-rows. General
  single-character borders by CP-SAT (`theory/fable_border7.py`, one worker, run alone): **UNKNOWN at the 7200 s
  limit** (2026-11-12, 7.3 GB). General borders (`theory/fable_border5.py`): not run to a verdict (memory). **Not a zero.**
- **σ-invariant borders** (σ the order-3 element that defines the slice; the compressions are σ-invariant, so M
  commutes with σ and a border may be sought whose row sets are unions of σ-orbits; `theory/fable_border10.py`,
  CP-SAT on orbit representatives, validated on a ℤ₂² Singer control where it finds and verifies a border in 9 s):
  - [P] ⟨z, σz⟩ = Σ_c ⟨f_c, f_c∘σ⟩, and each term is 16 for the six σ-fixed block functions and 0 for the other 24.
    So rows in σ-orbits of size 3 use only non-fixed block functions, and since 28 ≡ 1 (mod 3) every σ-invariant
    border has at least one σ-fixed row.
  - [P] **No σ-invariant single-character border exists.** A single-character orbit row avoids the trivial
    character entirely, while a σ-fixed single-character row uses it in every block. Each block where both carry the
    same character adds a term ≡ 2 (mod 4) to 4K(f, z) (Γ's entries are ±2, ±4), so |K| = 1 needs exactly two such
    blocks, and there are none. (`theory/fable_border11.py` confirms: no compatible pair.)
  - [M] For general block functions this gate does not close: all 331,776 admissible orbit rows are compatible with
    some σ-fixed row (`theory/fable_border12.py`). CP-SAT with the exact domain reduction (9.5 GB on 1 worker) runs
    with 1, 4 and 7 fixed rows, 2 h each; results in `logs/fable_fresh/sigma_queue.log`.
    <!-- SIGMA -->CP-SAT (2026-11-12): NFIX = 1, 4, 10 **UNKNOWN at 7200 s** each; NFIX = 7 **INFEASIBLE** (proof, 2803 s). Fixed rows lie in the 12-dimensional σ-invariant part of the border space, so NFIX ∈ {1, 4, 7, 10}.
    **Decided on 2026-11-13 by an exact clique search** (`theory/fable_border15.py`; method of §7.2: a σ-fixed first row x₁ of E, all of Z inside the rows compatible with x₁, Z's fixed rows chosen first, then a coupled search for the orbit rows of Z that leave room for the orbit rows of E; exhaustive, no cap hit): **NFIX = 1: none** (1009 s; a second run through the coupled code path, 800 s, agrees); **NFIX = 4: none** (57 s); **NFIX = 7: none** (2 s, agreeing with the CP-SAT proof); **NFIX = 10: none** (2 s). (A posteriori remark [P]: σ permutes the rows of E, so the number of fixed rows is the trace of σ on the border space, 4 · 1 = 4; only NFIX = 4 was ever possible, and the other three searches are consistency checks.) **So no σ-invariant border exists for case 0**, and hence for none of the 32 Q-level profiles that are constant on T (they are equivalent to case 0 under signed permutations of the blocks). The control (ℤ₂², Singer σ) is found by the same program in 0 s.<!-- /SIGMA -->
- **If a border exists**, the slice is a staged slice of the repository's usual shape: G′ = ℤ₂⁴ with 6 orbits,
  staging prime 37 with classes of index 12. Its stage 2 was not costed; with 78 orbits per block (312 sign bits)
  against 36 class equations of the degree-6 fixed field it is in the "rich" regime of `report/open_screen_answer.md`
  [E, rough], so even a border would leave a large search.
- 907 and 971 (|T| = 11) also have alive Q-levels with σ of order 5 (`tcensus2.jsonl`); their U have even order.

### 7.2 A universal symmetry of the border problem over ℤ₂^a, and a fast border search [P, M; added 2026-11-13]

**Observation [P].** Over Q = ℤ₂^a every character is real and inversion is trivial, so on the character space of
ρ the Goethals–Seidel array of the compressions is the 4 × 4 matrix GS(a, b, c, d) with (a, b, c, d) = (ŷ₁(ρ), …, ŷ₄(ρ)).
That matrix is the transpose of the left multiplication by the quaternion a + bi + cj + dk on ℝ⁴ = ℍ, the four
coordinates being the four blocks. It therefore commutes with the eight right multiplications by ±1, ±i, ±j, ±k,
which are signed permutations of the blocks. **So the group Q₈, acting on the block index, is a symmetry of the
border problem for every character set T ⊂ Q̂ and every core** (checked numerically in every run below: exactly 8
signed block permutations commute with a generic M, 24 with the spike compressions (1; 3, 3, 3)). It acts freely
on rows, a row is orthogonal to its three images, and K(xP, yP) = K(x, y). Over a quotient with complex characters
(ℤ₅, ℤ₆) no signed block permutation other than ±1 commutes with a random array [M], so this is special to
elementary abelian 2-groups (and ℤ₄).

**Reformulation [P].** Write a row as a function Ξ : Q → {±1 ± i ± j ± k} (one coordinate per block) with Fourier
support in T. A border whose row sets are unions of Q₈-orbits ("quaternionic border") is a pair of families
Ξ₁, …, Ξ_t and Ξ′₁, …, Ξ′_t, orthogonal for the quaternionic Hermitian product (1/k) Σ_q Ξ̄_m(q) Ξ_m′(q), such that
every κ_mm′ = ¼ Σ_ρ ξ̄_mρ γ_ρ ξ′_m′ρ lies in {±1 ± i ± j ± k} (γ_ρ the quaternion of the Walsh values at ρ, ξ_mρ the
Fourier coefficients of Ξ_m). K is then a t × t matrix over the sixteen Hurwitz units ½(±1 ± i ± j ± k), times 2:
t quaternionic rows on each side in place of 4t real ones. The weaker ansatz of invariance under one cyclic
subgroup ⟨i⟩, ⟨j⟩ or ⟨k⟩ of order 4 gives 2t complex rows and K/(1 + i) a quaternary complex Hadamard matrix.

**Search [M].** `theory/fable_border13.py` (Q₈) and `theory/fable_border14.py` (`sub=q8|z1|z2|z3`) enumerate orbit
representatives of all admissible rows. For a first row x₁ the whole of Z must lie in the set of orbits compatible
with x₁, which is small; the program lists the cliques of t (or 2t) mutually orthogonal orbits in it, each a complete
Z, and then looks for E among the orbits compatible with all of Z. It is exhaustive when it ends below the node cap.
Every border is checked in integers (K Hadamard, the compressed Gram identity), and for a real core the full matrix
is assembled and H Hᵀ = nI checked.

| case | T | n | CP-SAT before | symmetric search |
|---|---|---|---|---|
| ℤ₂², spike compressions (1; 3, 3, 3), s = 7 | three nontrivial | 124 | found in 9 s (σ-invariant model) | Q₈: found, 0 s |
| ℤ₂³ × ℤ₃, annealed core | {1, χ₁, χ₂, χ₄} | 112 | found in 323 s | Q₈: found, 0 s; H(112) assembled and verified |
| ℤ₂³ × ℤ₃, annealed core | {χ₁, χ₂, χ₄} | 108 | found | Q₈: **none (exhaustive)**; ⟨i⟩, ⟨j⟩, ⟨k⟩: found, 0 s each; H(108) verified |
| ℤ₂³ × ℤ₃, annealed core | {1, χ₁, χ₂, χ₃, χ₄} | 116 | UNKNOWN at 600 s | Q₈: found, 0 s; H(116) assembled and verified |
<!-- QB -->
| ℤ₂⁴, the 2396 slice of §7.1, case 0 (stands for the 32 constant Walsh profiles) | two 2-flats, 7 characters | 2396 | UNKNOWN at 7200 s | Q₈: **none, exhaustive** (101,250 orbits; two runs, without and with the reduction to 1097 first rows); ⟨i⟩, ⟨j⟩, ⟨k⟩: **none, exhaustive** (202,500 orbits each; 800 to 2200 s) |
<!-- /QB -->

So the symmetric search finds in under a second the borders that the generic model needed minutes for or did not
find, and the Q₈ ansatz is not always sufficient (the 108 case).

**All symmetries of order 3 of the 2396 border problem [P for the group, M for the searches].** For a Walsh
quadruple constant on T the problem depends only on Γ = GS(4, 2, 2, 2) (every signed permutation of (4, 2, 2, 2)
gives an equivalent problem: automorphisms of ℍ, multiplications by units, transposition). Its symmetries, apart
from dualities exchanging E and Z, form the group (Q₈ ⋊ ℤ₃) × Aut_Q of order 24 · 1152: the 24 signed block
permutations commuting with Γ (Q₈ and the rotation ρ : i → j → k about the axis of the imaginary part of
4 + 2i + 2j + 2k) and the 1152 affine maps of F₂⁴ that keep T (all of them keep the table of K-values [M]). Its
order is 2¹⁰ · 27, so a border with a nontrivial symmetry (modulo the sign −1) is invariant under an element of
order 3 or under an involution modulo sign. The subgroups of order 3 fall into five classes: Q-part a 3-cycle on
one of the two flats (a) or on both (d, the σ of the slice), with block part 1 or ρ (a, d, r, ra, rd).
`theory/fable_border16.py` decides each class in the normal form, for every admissible number of fixed rows:
<!-- O3 -->
The number of fixed rows is forced: it is the trace of the element on the border space [P].

| class | block part | Q-part | fixed rows (forced) | orbits of size 3 | result |
|---|---|---|---|---|---|
| d (= σ) | 1 | 3-cycles on both flats | 4 | 8 | **none** (48 s; also none for the impossible counts 1, 7, 10; third confirmation of §7.1, in other coordinates) |
| rd | ρ | as d | 1 | 9 | **none** (61 s) |
| ra | ρ | 3-cycle on one flat | 4 | 8 | **none** (16 s) |
| r | ρ | 1 | 7 | 7 | **none** (160 s) |
| a | 1 | 3-cycle on one flat | 16 | 4 | **none** (6 s, orbit rows first: of the 290,304 sets of four mutually orthogonal orbits none has even one orbit compatible with all four) |
<!-- /O3 -->
**So a border of this shape, if one exists, has no symmetry of order 3** [M], and its stabiliser in the group
above is a 2-group that meets Q₈ in ±1 only. The involutions modulo sign are the elements i, j, k of Q₈ (one class, by ρ; the three runs in the table above
are therefore one result obtained three times) and many classes with a nontrivial Q-part, which were not run.

---

## 8. The production run: skew 361 over ℤ₁₉ × ℤ₁₉ [E; a finished run without a hit excludes nothing]

`report/skew_census_answer.md` left one calibrated slice unfinished: v = 361, H = ⟨diag(4, 7)⟩, pattern
(1; 17, 23, 25), searched to 10,000 of 33,600 ranked split values. Its lists were not kept, so the run was restarted
with new Metropolis lists (seed 2, target 3·10⁷ per block) and the same pipeline (`theory/skew_kt.py`, Kumar's
`kt4.c`); launcher `scripts/fable_run361.sh`, log `logs/fable_fresh/p361_s2.log`. Expectation: 0.5 to 7 twin
classes in the whole lists, depending on whether Kumar's measured Metropolis yield (0.07) or the Gaussian count
applies.

<!-- RUN361 -->
**Result (2026-11-12 01:07 UTC): no solution.** Seed 2: fresh Metropolis lists of 4,194,304 (block A, complete);
8,997,706, 10,196,987 and 9,818,064 (B, C, D after twin reduction); Gaussian estimate inside the lists 10^0.83 [E].
32,916 split values seen on both sides; the top 2000 by Kumar's ranking were joined by `kt4` (3.4·10¹⁰ AB pairs,
6.0·10¹¹ CD pairs, 2.6 h on 4 cores): **0 solutions**. With the previous run (seed 1, other lists, 10,000 of 33,600
split values, no solution) this is the second miss on a slice whose expectation inside the lists is 0.5 to 7 twin
classes. A miss excludes nothing [heuristic], but two independent list draws without a hit move this slice
towards the expensive end; it is no longer the natural next run on 4 cores.
<!-- /RUN361 -->

---

## 8a. The tier-0 slices again: one new decision, and an a-priori gate [M, P; added 2026-11-13]

**Correction first.** I ran the eight tier-0 slices that `report/open_screen_answer.md` lists as "not run" before
noticing that the repository's later report `report/h4892_answer.md` (2026-11-06, §4) had already decided most of
them, by two methods each, including the non-empty stage 1 of 7572 and its death at the ℤ₄₃ compression. So the
table below is a rerun, not a discovery. It reproduces that report's counts exactly (same tool, method B, `--strong`,
planted controls; `scripts/fable_tier0_queue.sh`, `logs/fable_fresh/tier0_stage1.jsonl`). **What is new is one
decision (4156) and the gate after the table.**

<!-- TIER0 -->
| 4v (v) | slice | method B here | `h4892_answer.md` §4 | status |
|---|---|---|---|---|
| 4156 (1039) | k4, ℤ₄ × ℤ₇ × ℤ₃₇, s = 4, U = −1 | **0 solutions, 11,527,476 blocks, 11 row patterns, 5110 s** (the planted control was cut by a container restart) | undecided (B timed out at 3600 s) | **new: empty at stage 1** (one method) |
| 5708 (1427) | b1, ℤ₂₃ × ℤ₂ × ℤ₃₁, s = 10, U of order 22 | 0, 94 blocks, planted 3/3 | empty, two methods | agrees |
| 5036 (1259) | b1, ℤ₁₇ × ℤ₂ × ℤ₃₇, s = 12, U of order 16 | 0, 142 blocks, planted 3/3 | empty, two methods | agrees |
| 7572 (1893) | b1, ℤ₁₁ × ℤ₂² × ℤ₄₃, s = 6, U of order 10 | 1360 solutions, 12,532 blocks, planted 3/3 | 1360 classes by two methods; empty at the ℤ₄₃ compression (energy 7564 > 7396) [P + M] | agrees; a second, a-priori proof below |
| 5692 (1423) | b1, ℤ₃² × ℤ₂ × ℤ₇₉, s = 6, U = −1 | 0, 109,832 blocks, planted 1/1 | empty, two methods | agrees |
| 9260 (2315) | b1, ℤ₁₃ × ℤ₂ × ℤ₈₉, s = 8, U of order 6 | 0, 520 blocks, planted 1/1 | empty, two methods | agrees |
| 3788 (947) | b1, ℤ₁₁ × ℤ₂ × ℤ₄₃, s = 14, U of order 10 | 0, 354 blocks, planted 1/1 | empty, two methods | agrees |
| 4748 (1187) | k4, ℤ₄ × ℤ₈ × ℤ₃₇, s = 4, U = −1 | OverflowError (block cap) at 11 GB | undecided (same cap) | still undecided |
<!-- /TIER0 -->

4156 is a k4 slice with odd |K|, so by Theorem B a core there would have completed to a matrix at once; with it,
every tier-0 slice of the screen except 4748 is now decided empty.

**An a-priori gate [P].** The earlier report kills 7572 class by class after stage 1. The same conclusion follows
before any enumeration. The multiplier group of that slice contains the whole of ℤ₁₁\*, so every block has the form
x(a, g) = X⁰(g) for a = 0 and X¹(g) for a ≠ 0, with g in G₂ = ℤ₂² × ℤ₄₃ (172 elements). At a character (α, χ) with
α ≠ 0 the transform is X̂⁰(χ) − X̂¹(χ). All 172 of these characters need Σ_k |·|² = n, so by Parseval
Σ_k Σ_g (X⁰_k(g) − X¹_k(g))² = n = 7572, while the left side is at most 4 · 172 · 4 = 2752. In general:

  *if the blocks are invariant under all of ℤ_ℓ\* on a factor ℤ_ℓ of G = ℤ_ℓ × G₂ (ℓ prime), and the border
  characters are trivial on ℤ_ℓ, then 4ℓ|G₂| ≤ n ≤ 16|G₂|, so ℓ ≤ 3.*

It is the Parseval identity for the compression to G/ℤ_ℓ, whose alphabet is {±ℓ, ±(ℓ − 2)}. Of the 81 slices the
screen placed in tiers 0–3, **18 have a full unit group on a prime axis ℓ ≥ 5 and are dead by this gate without any
enumeration**: 4 of the 10 in tier 0 (947, 1259, 1427, 1893, all four also decided by enumeration), 1 of 5 in tier 1
and 13 of 53 in tier 2 [M: key scan of `logs/open_screen/table.jsonl`]. Over the whole screen, 214 of the 1283 slices that pass its gates A and B have such an axis (b0: 2 of 48; b1: 111 of 740; k4: 101 of 495) [M: key scan of `logs/open_screen/screens_s*.jsonl`].

**Tier 1 and tier 2.** The 44 slices of tiers 1 and 2 that survive the gate are being rerun the same way
(`scripts/fable_tier12_queue.sh`, 45 minutes and 7 GB each at most; `logs/fable_fresh/tier12_stage1.jsonl`). Several
of them are in `h4892_answer.md` already (4892 with its 384 classes and no core, 3004, 5892, 7180, …).
<!-- TIER12 -->
So far: 599 (b1, ℤ₁₃ × ℤ₂, q = 23, s = 2) empty; 1487 (k4, ℤ₄ × ℤ₇ × ℤ₅₃, s = 4) empty (1,450,382 blocks; planted
control undecided); 1103 (b1, ℤ₁₉ × ℤ₂, q = 29, s = 4) empty; 751 (k4) empty and 1223 with 384 solutions, both as in
the earlier report; 1987 (k4), 1831 and 1473 capped (undecided).
<!-- /TIER12 -->

---

## 9. What I would do next

1. **With about 10³ core-hours:** Kumar's Spence plug-in for 4792 = 8·599 over ℤ₇ × ℤ₁₉ × F₉ (his `search/eightp/gs`,
   calibrated on 5912 within a factor of two), then the skew slices he ranks at 10³ (661, 883, 1483). These are the
   only live targets with a measured rate.
2. **The 2396 slice of §7.1.** Its border question is now narrow: no border exists with Q₈, with a cyclic subgroup
   of order 4 of Q₈, or with any symmetry of order 3 (§7.2). What remains is the involutions with a nontrivial Q-part (a finite list of classes; the engine of
   `theory/fable_border15.py` / `fable_border16.py` handles each in minutes to hours), and borders with no symmetry
   at all, for which I see no exact method within reach (a rough count of the same search tree without symmetry
   gives about 10² core-days [E]). I would not spend more on it: even with a border the slice has 78 orbits per
   block, and a rough application of Kumar's cost formula (work ≈ N³/(E·L), his `METHOD.md`) gives 10²² list
   operations or more, against the 10¹³ he calls reachable [E, rough; his triage code was not run on it].
3. **A census of character borders with general T over ℤ₂^a** (§7): affine classes of T in F₂^a for a = 4, 5, 6, every group
   ℤ₂^a × S (and the non-split 2-parts) of order p − |T|, every multiplier group, the exact Q-level gate, and a
   border library per Walsh-quadruple type. This is the one set of structures found here that is both unscreened
   and cheap to screen with the repository's existing gates.
4. **Borders.** §6.3 turns "does a border exist for this compression" into a concrete question about pairs of ±1
   bases. A theorem for a class of non-spike compressions at k = 6, 8 or 10 would open the mid-range indices
   (k = 10–60 with |S| = 11–59 at 599 and 719), which no census has touched because no border is known there.
5. **Formal multiplications at t ≥ 7** for the dimension-5 types outside the Williamson pattern. The count does not
   exclude a sporadic identity; t = 7 needs a three-level join (7·10⁵ block classes).
6. **Borders as a tool.** The Q₈ search of §7.2 is a border finder for any character set over ℤ₂^a. Kumar's notes
   list borders for spike compressions at index k ≥ 10 as not found (and single characters as impossible for
   k ≥ 12); a Singer cycle times Q₈ reduces that case to one quaternion-valued sequence of length k − 1 on each
   side, which is a small search I did not run.
7. **A gate to add to the screen.** Before costing a slice, test every prime axis for a full unit group (§8a), and
   more generally the Parseval value of every compression against its alphabet, not only the one along ℤ_q.
8. **Not worth repeating:** product theorems of the two kinds closed in §3 and §4; Theorem 7 with |S| ≤ 7; the
   tier-0 list of the screen (all decided empty except 4748).

## 10. Reproduce

```
python3 theory/fable_tmult_types.py 3            # all 164 types at t = 3
python3 theory/fable_tmult_h.py 5 "f s s s" "f f s =2" "f f s =0" "f f s *0"
python3 theory/fable_tmult_sat.py 5 generic 600  # CP-SAT, INFEASIBLE
python3 theory/fable_cod4.py all_sym AB_sym A_sym D_eq_Cstar generic
python3 theory/fable_octo.py                     # 2048 compositions, pattern sets E / O
python3 theory/fable_tamic.py 3 5 7 9 11
python3 theory/fable_border4.py 3 3 0,1,2,4 1 380 1   # a non-classical character set, H(112)
python3 theory/fable_border2.py 3 3 1,2,4 2 300 1       # twisted k4 with even |K|, H(108)
scripts/fable_run361.sh 4                        # needs /home/user/skew_work/bin (see logs/fable_fresh/NOTES.md)
python3 theory/fable_clique.py                   # self-test of the compiled clique routines
python3 theory/fable_border14.py sub=q8 core 3 3 0,1,2,3,4 1   # five characters: border in under a second, H(116) assembled
python3 theory/fable_border14.py sub=z2 core 3 3 1,2,4 2       # H(108): no Q8 border, a cyclic one exists
python3 theory/fable_border14.py sub=q8 q599 0                 # 2396 slice: no Q8-invariant border (about 6 minutes)
python3 theory/fable_border15.py q599 0 4                      # no sigma-invariant border (4 fixed rows; 1 minute)
python3 theory/fable_border16.py rd                            # order-3 classes in normal form: d, r, ra, rd, a
scripts/fable_tier0_queue.sh                                   # the tier-0 reruns of section 8a
```
