# Theory arm for the open orders 2092–2956

Composition closure, the Frobenius group at 669, algebraic slices, and why no ≤ 2^45 decoupling was found.

Date: 2026-10-08. Branch `claude/focused-planck-j1645f`. Code: `theory/`. Logs: `theory/logs/`.

**Tags:**
- **[M]** verified by exact computation in this session. The script is named.
- **[P]** proved here. The proof is in §6.
- **[L]** literature recollection that I could not re-check: arXiv and several other hosts are blocked by this session's network policy.
- **[E]** estimate (Gaussian model, same conventions as the repo).

---

## 0. Bottom line

1. **No Hadamard matrix of any of the eleven orders was constructed.** Nothing here closes an order.

2. **The composition closure is empty.** No multiplication, plug-in or array theorem I know produces H(4v) for any open v from
   ingredients that exist in October 2026 (§1, table 1.2). The reasons:
   - Every multiplication theorem outputs orders divisible by 8.
   - For the eight prime v, every plug-in theorem needs an object of order v itself.
   - For v = 3p (537, 573, 669), Cooper–Wallis and Baumert–Hall need Williamson-type matrices of order p, or T-matrices
     (equivalently a BH array) of order p. None is known.

3. **The 2026 ingredients do not change this [M].**
   - Aut(H(716)) = Z₂ × Z₂ and Aut(H(892)) = Z₆, computed with Traces from nauty 2.8.8; every element was checked as P·H·Qᵀ = H.
   - A matrix equivalent to a Williamson-type array must contain Q₈ in its automorphism group (Lemma A′). Neither group
     does.
   - So the matrices closed in 2026 cannot serve as the Williamson-type ingredient for 537 or 669.

4. **669 and Z₂₂₃ ⋊ Z₃: no new leverage [P].**
   - The Goethals–Seidel mechanism needs commuting blocks with transposition realised by conjugation. Over a non-abelian
     group this forces the four blocks into a commutative subalgebra of R[F].
   - In the natural maximal one, C[Z₂₂₃] + Z(C[F]), which contains every class function, a Parseval count rules out a
     solution: the degree-3 irreducibles need spectral weight 594 072 and at most 198 912 is available (Lemma D).
   - As a symmetry of an abelian core, the Z₃ action is the order-3 multiplier, which T(223, k = 3) already uses
     (≈157 bits).
   - The fixed-point-free "Frobenius" affine maps of Z₆₆₉ are row-sum infeasible, because 3 would divide every block sum
     and 9 ∤ 2676 (Lemma E2). Their block-twisted version leaves one block unconstrained, so it gives no reduction.

5. **Algebraic (high-symmetry) slices: exhaustive and empty [M].** "Slice" means the four (or two) sequences are invariant
   under a symmetry subgroup S ≤ Aut(G). Every such subgroup with ≤ 30 orbits was enumerated, using Goursat's lemma for
   the non-abelian part, and each slice was solved exhaustively. There are two exceptions: one b = 1 slice for 573 is
   undecided (join cap; Gaussian log₂E ≈ −58), and Z₂×Z₄×Z₆₇ was swept with diagonal automorphisms only.

   | family | slices | solutions |
   |---|---|---|
   | Legendre pairs at ℓ = 2v − 1 (all eleven ℓ, both abelian groups of order 1197) | 832 | 0 |
   | classical bordered GS (b = 1) over every abelian group of order v − 1 | 391 | 0 (1 slice undecided, log₂E ≈ −58) |
   | Goethals–Seidel (b = 0) over Z_v | 62 | 0 |

   - The b = 1 groups include five non-cyclic groups that the repo never costed: Z₂³×Z₆₇, Z₂²×Z₁₁×Z₁₃, Z₂²×Z₁₆₇,
     Z₂×Z₃²×Z₂₉ and Z₂×Z₃²×Z₄₁.
   - **Positive controls:** the same code finds and certifies, by assembling the full matrix and checking H·Hᵀ = N·I:
     - H(28), H(32), H(60) and H(72) via Legendre pairs at ℓ = 13, 15, 29 and 35;
     - H(36) and H(52) via b = 1 over Z₂³ and Z₂²×Z₃;
     - H(340) via b = 1 over Z₂²×Z₃×Z₇ with a symmetry group of order 6, the same shape as the 573 group;
     - H(52) and H(76) via b = 0.

6. **New or clarified lemmas [P]** (§6):
   - **A.** The Williamson array works iff three Δ-matching identities hold. This is weaker than pairwise amicability.
   - **B.** The Cooper–Wallis plug-in criterion:
     - pairwise-amicable T-matrices would accept any GS quadruple;
     - no amicable T(3) exists (all 192 T(3) checked);
     - with T(3) the requirement is a single antisymmetric circulant identity, still (p−1)/2 extra constraints.
   - **C.** Williamson-type matrices of odd order > 1 can never satisfy the T-matrix product condition, so T-type and
     W-type are disjoint.
   - **D.** The Frobenius-group impossibility in point 4.
   - **E.** The row-sum lattice lemma. It explains why LP and b = 1 tolerate high symmetry and b = 0 does not.
   - **F.** Legendre ⊕ D-optimal ⇒ GS: for prime v ≡ 3 mod 4, a circulant D-optimal pair of length v gives H(4v). Among the
     targets only 643 passes the arithmetic.
   - **G.** A support-split counting lemma that says exactly when a T(103)-style decoupling pays.
   - **H.** Order-2 block twists at primes v ≡ 3 mod 4 force the overdetermined propus shape, or an impossible
     odd-length periodic Golay pair.

7. **(C) ansatz ≤ 2^45: none found.** Lemma G together with H explains why the natural decouplings fail. The best new
   structured configurations are the non-cyclic b = 1 slices (§4.5):

   | v | group, symmetry | cost [E] | old best | status |
   |---|---|---|---|---|
   | 573 | Z₂²×Z₁₁×Z₁₃, S of order 15 | n ≈ 167 | 193 (repo) | survives every exact screen small enough to run |
   | 739 | Z₂×Z₃²×Z₄₁, \|S\| = 20 | n ≈ 187, log₂E ≈ +2.5 | ≈ 200–210 | killed by the exact Z₂×Z₃² quotient screen [M]; next survivor n ≈ 282 |
   | 537 | Z₂³×Z₆₇, Z₁₁ | n ≈ 179 | 290 (repo) | killed by an exact Z₂³ compression screen [M] |

   The 537 and 739 slices are killed by exact screens, so the only surviving improvement is 573 (167 vs 193). It is
   still far beyond the ~90–100 bits reachable with lists of 2^28.

---

## 1. Composition closure

### 1.1 Constructions that output H(4v) with v odd

An order 4v with v odd is never a product of two Hadamard orders larger than 2. So every route to it is either direct,
or a plug-in whose ingredients have orders dividing v.

| # | construction | output | hypothesis | tag |
|---|---|---|---|---|
| T1 | Paley I | H(q+1) | q ≡ 3 (mod 4) prime power; here q = 4v − 1 | [L] |
| T2 | Paley II | H(2q+2) | q ≡ 1 (mod 4) prime power; here q = 2v − 1 | [L] |
| T3 | symmetric conference matrix C of order 2v | H(4v) = [[C+I, C−I],[C−I, −C−I]] | a conference matrix of order 2v exists. Necessary: 2v − 1 is a sum of two squares. Known orders: q + 1, Mathon's q²(q+2)+1, a few sporadic | [L] |
| T4 | Williamson / Williamson-type of order v | H(4v) | families: Turyn (q+1)/2; Whiteman–Mukhopadhyay type q^r(q+1)/2; Miyamoto–Seberry–Yamada order q ≡ 1 (mod 4) | [L] |
| T5 | Miyamoto | H(4q) | q ≡ 1 (mod 4) prime power and H(q − 1) exists | [L] |
| T6 | Seberry–Whiteman | H(4(q²+q+1)) | arithmetic conditions on q | [L] |
| T7 | Goethals–Seidel | H(4v) | four ±1 sequences over an abelian group of order v with ΣPAF = 0 (SDS) | classical |
| T8 | Cooper–Wallis | H(4tw) | T-matrices of order t and Williamson(-type) of order w | [L] |
| T9 | Baumert–Hall, Welch, Ono–Sawade–Yamamoto arrays | H(4tw) | BH array OD(4t; t,t,t,t) and Williamson-type of order w | [L] |
| T10 | T-sequences ⇒ T-matrices | T(t) | sources: BS(n+1, n) gives t = 2n+1; TT(n) ⇒ BS(2n−1, n) gives t = 3n−1; Golay gives t = g+1; Yang gives t = (2n+1)·y. Known: BS(n+1,n) only for n ≲ 40, TT(n) only for even n ≤ 36 | [L] |
| T11 | Legendre pair of length ℓ (two G-developed cores) | H(2ℓ+2) | ℓ = 2v − 1. Families: ℓ a prime power, ℓ = 2^k − 1, ℓ = q(q+2) twin prime powers | [L]; builder verified [M] |
| T12 | classical bordered GS (b = 1) | H(4v) | four ±1 sequences over an abelian group of order v − 1 with ΣPAF = −4 off 0 | classical; verified for non-cyclic G [M] |
| T13 | inflated bordered GS (b ≥ 3, 2026) | H(4v) | even c ≥ 4 with c \| v + 1, and a core over Q × Z_h | repo |
| T14 | multiplication theorems: Agaian 8mn, Craigen–Seberry–Zhang 16mnpq, Seberry–Yamada "W-type u, v ⇒ W-type 2uv", Craigen signed groups, Kronecker | orders ≡ 0 (mod 8) | – | [L] |
| T15 | Seberry–Yamada M-structures | H(4mn) from an M-structure of order 4m and W-type of order n | definition not retrieved | [L] |
| T16 | cocyclic / Ito Hadamard groups Q₈ₜ | H(4t) | Q₈ₜ known to be a Hadamard group when 2t−1 or 4t−1 is a prime power, and by computer for t ≤ 46 (Schmidt) | [L] |
| T17 | Legendre ⊕ D-optimal (Lemma F, elementary) | H(4v), v prime ≡ 3 (mod 4) | a circulant D-optimal pair of length v. Needs 2v − 1 = a² + b² | [P] |

**Why the open primes are all ≡ 3 (mod 4)** [L + arithmetic]. For a prime v ≡ 1 (mod 4), Miyamoto (T5) needs only H(v − 1),
and v − 1 < 2000 is a known order. So every prime ≡ 1 (mod 4) in this range is closed, and the open primes are exactly the
primes ≡ 3 (mod 4) that no direct family covers. The open composite v = 3p have p ∈ {179, 191, 223}. All three p are
≡ 3 (mod 4), and these were the hard orders 716, 764 and 892.

### 1.2 Per-order closure table

Arithmetic is from `theory/closure_table.py` [M]. "–" means the construction does not apply on arithmetic grounds alone.

| v (4v) | T1 4v−1 | T2 2v−1 | T3 conf(2v) | T4–T6 families | T7–T10 splits t·w | T11 LP(2v−1) | T12 b=1 groups (v−1) | T13 c | T17 |
|---|---|---|---|---|---|---|---|---|---|
| 523 (2092) | 3·17·41 – | 5·11·19 – | – | – | prime | 5·11·19, no family | Z₂×Z₉×Z₂₉, **Z₂×Z₃²×Z₂₉** | 4 | – |
| 537 (2148) | 19·113 – | 29·37 – | 1073 = 7²+32² = 17²+28², no known conference matrix | – | (3,179), (179,3) | 29·37, no family | Z₈×Z₆₇, Z₂×Z₄×Z₆₇, **Z₂³×Z₆₇** | – | composite |
| 571 (2284) | 3·761 – | 7·163 – | – | – | prime | 7·163 | Z₂×Z₃×Z₅×Z₁₉ (cyclic only) | 4, 22, 26, 44, 52 | – |
| 573 (2292) | 29·79 – | 5·229 – | 1145 = 11²+32² = 19²+28², none known | – | (3,191), (191,3) | 5·229 | Z₄×Z₁₁×Z₁₃, **Z₂²×Z₁₁×Z₁₃** | 14, 82 | composite |
| 599 (2396) | 5·479 – | 3²·7·19 – | – | – | prime | 3²·7·19: Z₉×Z₇×Z₁₉ or **Z₃²×Z₇×Z₁₉** | cyclic only | many | – |
| 643 (2572) | 3·857 – | 5·257 – | 1285 = 14²+33² = 18²+31², none known | – | prime | 5·257 | cyclic only | 4, 14, 28, 46, 92 | **yes**: needs D-opt(643) |
| 647 (2588) | 13·199 – | 3·431 – | – | – | prime | 3·431 | cyclic only | many | – |
| 659 (2636) | 5·17·31 – | 3·439 – | – | – | prime | 3·439 | cyclic only | many | – |
| 669 (2676) | 5²·107 – | 7·191 – | – | – | (3,223), (223,3) | 7·191 | Z₄×Z₁₆₇, **Z₂²×Z₁₆₇** | 10, 134 | composite |
| 719 (2876) | 5³·23 – | 3·479 – | – | – | prime | 3·479 | cyclic only | many | – |
| 739 (2956) | 3·5·197 – | 7·211 – | – | – | prime | 7·211 | Z₂×Z₉×Z₄₁, **Z₂×Z₃²×Z₄₁** | 4, 10, 20, 74, 148 | – |

Notes:
- **T4–T6 [M, arithmetic].** No v equals (q+1)/2, q^r(q+1)/2, q² + q + 1 or 2q² + 2q + 1 for a prime power q. The last
  two are the D-optimal families relevant to T17.
- **T10 [M].** For every t ∈ {v, 179, 191, 223}: n = (t−1)/2 ≥ 89, which is beyond the known BS(n+1, n). t − 1 is never
  2^a·10^b·26^c. (t+1)/3 is 60, 64, 200, 216, 220 or 240, never an even n ≤ 36.
- **Bold groups** are the non-cyclic groups whose symmetry groups were costed and swept here for the first time (§4.5).
- **Hadamard difference sets (4v−1, 2v−1, v−1)** (T1 and its relatives) [M, arithmetic]:
  - Paley needs 4v − 1 prime; all eleven are composite.
  - Singer needs 4v − 1 = 2^k − 1; impossible.
  - Twin prime powers need 4v − 1 = q(q+2); no.
  - Hall's sextic needs 4v − 1 = 4x² + 27 prime; again 4v − 1 is never prime.
- **Aperiodic complementary quadruples of length v (BS(v, v))** would give a GS quadruple directly. They are known only
  for small or Golay-type lengths [L], and no v here is 2^a·10^b·26^c [M].
- **Non-abelian groups of order 2v − 1** exist for ℓ = 1045, 1197, 1317 and 1477. A two-core LP over such a group
  needs commuting cores (left times right regular) whose autocorrelations are class functions. That adds ℓ − k
  constraints per sequence (k = number of classes), so it is hopelessly overdetermined [P, same argument as §3].

#### The (t, w) splits for v = 3p, ingredient by ingredient

| route | needed ingredient | status, October 2026 |
|---|---|---|
| Cooper–Wallis or BH(12) with w = p | Williamson-type of order p | unknown [L]. The 2026 matrices H(716) and H(892) are **not** Williamson-type [M, §2]. H(764) (Djoković 2008, GS type) is not in the repo; for circulant blocks, W-type needs symmetric pairwise cross-correlations. Searching symmetric circulant Williamson directly (§3, item 6): for p = 223 with K = ⟨−1, λ⟩, n ≈ 169 bits and log₂E ≈ −34 [E], overdetermined. For p = 191 every K ∋ −1 of order 10, 38 or 190 is row-sum infeasible [M] (exact DP: \|K\| = 10 allows only \|a\| ∈ {9, 11} ≤ 2√191); K = {±1} leaves a ≈ 380-bit search. For p = 179 only K = {±1} is available (≈ 360 bits). |
| … with the weaker Lemma B condition | GS quadruple of order p with one antisymmetric identity G₁₂ + G₂₃ + G₃₁ = 0 | about 2× the GS constraint count. Overdetermined (§6, Lemma B). |
| Cooper–Wallis with t = p | T(p) or BH(4p) | T(191): dead for every odd-order multiplier group (repo [M]). T(179): no usable odd subgroup (178 = 2·89). T(223, k=3): about 157 bits [E]. No BH(4p) is known without T-structure [L]. |
| M-structure × W(3) | M-structure of order 4p | [L]. Unresolved here, because the definition could not be retrieved. If every GS-type H(4p) were an M-structure, Djoković's H(764) would have closed 573 in 2008, so presumably they are not. A plug-in of the form Σ_k T_k ⊗ A_k with W-type A_k provably needs (T_k) to be an OD(4m; m,m,m,m) (§6). |
| Ito / cocyclic Q₈·₆₆₉ | – | no theorem [L] |

**The inference "v open, therefore ingredient Y unknown", re-checked:**
- **179, 223:** the 2026 matrices exist, but they are the wrong kind of object (GS / bordered, not W-type or BH).
- **103:** T(103) divides no open v.
- **H(2060), H(2884):** 515·w and 721·w exceed 750 for every w > 1.
- So no 2026 ingredient feeds any composition theorem for the eleven v.

---

## 2. Re-checked ingredients: automorphism groups of the 2026 matrices [M]

**Method** (`theory/hadaut_traces.c`, `run_aut.py`, `check_autgens.py`):
- The graph has vertices r±, c± and an edge r^ε – c^δ iff ε·δ·H_rc = +1, with rows and columns in different colour classes.
- Its automorphism group is the group of signed-permutation pairs (P, Q) with P·H·Qᵀ = H.
- Traces from nauty 2.8.8 was built from source. Plain nauty did not finish within 12 minutes on H(716).

| matrix | source | \|Aut\| | elements | time |
|---|---|---|---|---|
| H(716) | `data/hadamard_716.csv.gz` (b = 1 core, Z₁₇₈) | 4 | τ = (−I, −I); one involution fixing 4 rows up to sign (the order-2 translation); their product | 402 s |
| H(892) | `data/hadamard_892.csv.gz` (b = 0 core, Z₂₂₃, d = 3) | 6 | τ; the order-3 multiplier, fixing the 4 index-0 rows; their products | 1203 s |

Every element was re-verified as P·H·Qᵀ = H (`check_autgens.py`, generators in `theory/logs/aut*_traces_gens.txt`).

**Consequence.** A Williamson-type matrix (A, B, C, D pairwise amicable) gives an H(4n) that commutes with
I ⊗ ρ(i) and I ⊗ ρ(j), right multiplication by the quaternion units. So its automorphism group contains a Q₈ whose central
involution is τ (Lemma A′). Neither group above has order divisible by 8. **Neither H(716) nor H(892) is equivalent to a
Williamson-type matrix.**

Note that a Goethals–Seidel matrix does not even carry the group translations as automorphisms. Only elements g with
2g = 0 survive, because the R-twist in columns 2–4 turns g into −g. That is why |Aut| is so small.

---

## 3. 669 = 3·223 and the Frobenius group F = Z₂₂₃ ⋊ Z₃

1. **What a GS-type construction over F would need.** The Goethals–Seidel orthogonality uses two facts:
   - the four blocks commute;
   - there is an R with R·Xᵀ·R⁻¹ = X on the blocks.

   If R·L(x)ᵀ·R⁻¹ = L(x) on a set closed under products, that set is commutative, because transposition reverses
   products and conjugation does not. So the four F-developed blocks must lie in a commutative ∗-subalgebra of R[F].
   - With R the inversion permutation, this forces the blocks to be central (class functions).
   - The largest "obvious" commutative subalgebra is A = C[Z₂₂₃] + Z(C[F]), of dimension 225 and containing every class
     function. It contains only 2^225 ±1 functions, against 2^669 for Z₆₆₉.

2. **Lemma D [P].** No four ±1 functions in A satisfy the GS condition. On each degree-3 irreducible ρ_β the outer
   class sums act as 0, so ρ_β(x) = diag(f̂(β), f̂(λβ), f̂(λ²β)), where f = x restricted to Z₂₂₃. The GS condition needs
   Σ_k Σ_{β≠0} |f̂_k(β)|² = 74·3·4v = 594 072. Parseval on Z₂₂₃ caps the left side at 4(p² − 1) = 198 912.

   The general (noncommutative) Parseval budget is exactly consistent: 4v² = 12v + 74·36v. So the obstruction comes
   precisely from the commuting requirement.

3. **The Z₃ action as a symmetry of an abelian core.**
   - The multiplier λ (λ³ = 1 mod 223) is the symmetry already used by T(223, k = 3), at about 157 bits.
   - The fixed-point-free affine maps (x, y) ↦ (x+1, λ^j y) of Z₃ × Z₂₂₃ have every orbit of size divisible by 3. That
     makes every block sum divisible by 3, and 9 ∤ 4v = 2676, so they are infeasible (Lemma E2). This holds for 537 and
     573 as well, since 9 ∤ 12p.
   - With the block 3-cycle twist B → C → D, row sums become feasible: (33, 23, 23, 23) or (51, 5, 5, 5), from
     3a′² + a_B² = 892. But B is then constrained only by the identity: one free block of 669 bits and about 111
     constraint classes, so n ≈ 550 bits. No reduction.

4. **Bordered and cocyclic options.**
   - b = 1 needs a core of order 668. The only odd symmetries are Z₃ (on Z₂²) and Z₈₃. Every slice has r ≤ 12 or
     r ≥ 168.
     - All r ≤ 30 slices are exhaustively empty [M].
     - The best mid-size slice has log₂E = −14 [E].
   - b ≥ 3 needs c ∈ {10, 134}; the repo's cheapest is about 590 bits.
   - The Ito/cocyclic route has no theorem at t = 669 [L].

5. **Symmetric Williamson quadruples of order 223** (amicable, so they would give 669 via Cooper–Wallis with T(3)).
   - With K = ⟨−1, λ⟩ (order 6), the row sums must be ≡ ±5 (mod 12). The only profile is (17, 17, 17, 5).
   - Gaussian cost (`theory/cost_b0.py`): n ≈ 169 bits, log₂E ≈ −34 for ordered tuples. This is overdetermined; no
     solutions are expected.
   - Calibration: the same code gives n ≈ 152, log₂E ≈ +130 for the non-symmetric d = 3 GS profile (23, 11, 11, 11) of
     the published 223 core. The repo's figure is 156.6.

6. **Verdict.** Every known GS or bordered mechanism needs commuting blocks, so it cannot use F, and F supplies no
   symmetry beyond the order-3 multiplier. The best 669 route remains T(223, k = 3) with Cooper–Wallis(W(3)), at about
   157 bits [E]. The prompt's alternative, H(892) as a Williamson-type ingredient, is closed by §2.

---

## 4. Algebraic families and the structure of the search space

### 4.1 Legendre pairs at ℓ = 2v − 1

- **No known family applies** [M, arithmetic]: no ℓ is a prime power, 2^k − 1, or q(q+2).
- **Overdetermination** [E]. A pair has 2ℓ bits and (ℓ−1)/2 constraints, each costing about ½·log₂(4πℓ) − 2 ≈ 4.8 bits
  at ℓ ≈ 1100. Then log₂E ≈ 2ℓ − 2.4ℓ < 0, and the crossover is at ℓ ≈ 330.
  - A K-invariant slice divides both terms by |K|, so it stays overdetermined. LPs at ℓ = 1045–1477 can only come from
    algebra.
  - This is consistent with the computational frontier: open LP lengths start at 115.
- **Exhaustive slices [M]** (`sweep_b1.py … lp`):
  - Every subgroup K ≤ Aut(Z_ℓ), and for ℓ = 1197 every subgroup of GL(2,3) × (Z/7)* × (Z/19)* acting on Z₃²×Z₇×Z₁₉,
    with ≤ 30 orbits was solved exactly: PSD filter, then an exact hash join.
  - **0 Legendre pairs in all 832 slices** (table in §7).
  - Positive controls: ℓ = 13, 29 (Legendre), 15 = 3·5 and 35 = 5·7 (twin primes) give certified H(28), H(60), H(32)
    and H(72).
- **Quadratic generalised cyclotomy for ℓ = pq** [P, sketch in §6]. With A, B built from χ_p·χ_q on units and constants
  on the axes, matching the axis lags forces |q − p| ≤ 2. This recovers exactly the twin-prime construction
  (q − p = 2). Our ℓ have prime gaps of at least 8, so none qualifies.

### 4.2 Conference matrices of order 2v

- Only 537, 573 and 643 pass the two-squares condition.
- A group-developed symmetric core would be a Paley-type PDS in an abelian group of non-prime-power order (1073, 1145,
  1285). As I recall, Ma's results confine Paley-type PDS in abelian groups to prime-power orders or orders n⁴ and 9n⁴,
  and these three are neither [L].
- Two-circulant cores of order v − 1 (even) are overdetermined like LPs [E].
- No family applies.

### 4.3 Legendre ⊕ D-optimal (Lemma F) and 643

- For prime v ≡ 3 (mod 4), the Legendre sequence L has PAF ≡ −1 off zero, whatever the value at 0.
- So GS(L, L, C, D) is Hadamard as soon as (C, D) is a circulant D-optimal pair: PAF_C + PAF_D = 2 off 0, which forces
  s_C² + s_D² = 4v − 2.
- **The only target with 2v − 1 a sum of two squares is v = 643.**
- **Status of D-optimal(643):**
  - It is overdetermined like an LP [E].
  - Neither known family (v = q² + q + 1 or v = 2q² + 2q + 1) contains 643 [M, arithmetic].
  - Symmetric versions are row-sum infeasible for |K| ≥ 107 (Lemma E). For |K| = 3 the pair has 428 free bits against
    roughly 640 constraint bits [E].
- This is a clean reformulation, not a route.

### 4.4 Row sums decide which architectures tolerate symmetry (Lemma E) [P]

| architecture | row-sum condition | consequence |
|---|---|---|
| b = 0 | Σ_k a_k² = 4v, so typical \|a_k\| ≈ √v | Orbits of size D outside the fixed points make a_k ≡ (fixed part) mod 2D (for an even orbit count). Once 2D exceeds about 2√v − f (f = number of fixed points), the a_k collapse to ±1 or ±3 and Σa² ≪ 4v. |
| b = 1 | (±2, 0, 0, 0), or all ±1 | no such barrier |
| LP | ±1 each | no such barrier |

So b = 0 high-symmetry slices die on row sums at d ≳ 25–40 for v ≈ 520–740. That is the pattern in the repo's `nbits` ∞
entries and in the 669 analysis. b = 1 and LP are the natural homes for algebraic, high-symmetry constructions, which is
why both were swept exhaustively here.

### 4.5 b = 1 over the non-cyclic groups of order v − 1

**Exhaustive slices with ≤ 30 orbits [M]** (`sweep_b1.py`):

| group | target | slices | solutions |
|---|---|---|---|
| Z₂³×Z₆₇ | 537 | 80 | 0 |
| Z₂²×Z₁₁×Z₁₃ | 573 | 48 | 0 |
| Z₂²×Z₁₆₇ | 669 | 7 | 0 |
| Z₂×Z₃²×Z₂₉ | 523 | 55 | 0 |
| Z₂×Z₃²×Z₄₁ | 739 | 65 | 0 |
| Z₂×Z₄×Z₆₇ (diagonal automorphisms only) | 537 | 7 | 0 |

The cyclic groups were also swept (§7). One cyclic slice (Z₈×Z₆₇, S of order 22, r = 29) exceeded the join cap. It is
**decided infeasible** by the exact Z₈-compression screen [M] (`screen_undecided_537.py`), in agreement with the repo's
earlier verdict.

**Mid-size slices: Gaussian cost** [E] (`cost_sweep.py`, same model as `b1_multiplier_enum.py`), then the exact quotient
screen [M] (`compress_screen.py`):

| v | G | \|S\| | r | L | n [E] | log₂E [E] | exact screens | repo best for this v |
|---|---|---|---|---|---|---|---|---|
| 537 | Z₂³×Z₆₇ | 11 | 56 | 31 | 179 | +24.9 | **killed** by the Z₂³ compression (also by hand: ΣY² = 1880 over 32 entries forces ≤ 4 large entries, so some lag cannot reach −268) | 290 (b = 1, d = 6) |
| 537 | Z₂³×Z₆₇ | 3 | 180 | 91 | 440 | +259 | survives Z₂³ | – |
| 573 | Z₂²×Z₁₁×Z₁₃ | 15 | 54 | 27 | 167 | +20.8 | **survives** the quotients Q = Z₂², Z₁₁ and Z₂²×Z₁₁; Q = Z₁₃, Z₂²×Z₁₃ and Z₁₁×Z₁₃ exceed the enumeration cap (undecided) | 193 (b = 1, d = 15) |
| 523 | Z₂×Z₃²×Z₂₉ | 14 | 54 | 29 | 191 | +3.5 | survives Q = Z₂×Z₂₉ and Q = Z₂₉; the others exceed the cap | 142 (b = 0, d = 9) |
| 739 | Z₂×Z₃²×Z₄₁ | 20 | 54 | 29 | 187 | +2.5 | **killed** by the Z₂×Z₃² quotient (H = Z₄₁) [M]; next survivor n ≈ 282 | ≈ 200–210 (b = 0, d = 9) |

**Screen validation [M]** (`screen_control.py`): 34 screens on slices with known solutions, 0 false negatives.

**Reading.** The non-cyclic groups improve the best known configuration for 573 only: 167 against 193 bits, and the slice
survives every exact screen small enough to run. Their best slices are killed by exact screens for 537 and 739, and
523's best slice (191) is worse than its b = 0 route. Nothing comes near 100 bits.

---

## 5. Is there a structured ansatz at ≤ 2^45? What fails and why

1. **The cost model.** With memory capped at lists of 2^28, the bucketed 4-list k-tree costs about 2^(n − 54)
   operations: one bucket handles about 84 bits, and each extra bit doubles the work. That matches the repo (2^35 per
   hit at the flat-restricted n′ ≈ 91). So 2^45 means n ≈ 99 effective bits. Every configuration above is at 140 bits or
   more; only the repo's flat restriction (523) and the twisted 571 d = 15 (94.6 bits, but lopsided) are in range.

2. **When does a T(103)-style split pay? (Lemma G.)** Suppose a problem has B free bits and n constraint bits
   (log₂E = B − n), and fixing F bits (a "support") leaves independent lists. Then each support choice carries
   2^(B − F − n) solutions in expectation, so the split pays only if F ≤ B − n.

   | problem | B | n | F | B − F − n | outcome |
   |---|---|---|---|---|---|
   | T(103), k = 3 | 105 (3 bits per orbit) | ≈ 63 | 35 | ≈ +7 | works, as observed |
   | b = 1 over Z₂×Z_h, P/R split (P = (X₀+X₁)/2 and R = (X₀−X₁)/2 are ternary with complementary supports on Z_h; given the supports, the problem splits exactly into two independent 4-list problems, ΣPAF_P = −2 off 0 and ΣPAF_R = v·δ) | – | – | 4h | ≈ −0.3v to −0.45v [E] | fails; symmetry divides every term by d but raises n per class, so it stays negative |
   | b = 0, symmetric/antisymmetric split (PAF_B = PAF_{sym B} + PAF_{antisym B} exactly, giving 8 lists), at 523/9 | ≈ 236 (59 orbits per block) | 142 | 116 | ≈ −22 [E] | fails; splitting only two blocks (F = 58) keeps E′ ≈ 2^32, but the pieces (2^15) are far too small for an 8-list k-tree at 142 bits |

3. **Order-2 twists are what would decouple, and the primes forbid them (Lemma H) [P].**
   - With a block swap (B ↔ C under a multiplier u), (B ± C)/2 have disjoint supports, which gives exact 4-list
     independence.
   - But for prime v ≡ 3 (mod 4), every even-order multiplier contains −1. That forces A and D to be symmetric and
     C = B̄ (propus shape, about 2v bits against about 2.45v constraint bits: overdetermined [E]). A double swap forces
     PAF_A + PAF_B = 0, a periodic Golay pair of odd length, which is impossible.
   - All eight open primes are ≡ 3 (mod 4). At 537, 573 and 669 (≡ 1 mod 4), the multiplier u ≡ −1 mod 3,
     u ≡ w mod p avoids −1. That gives about 7/12 of the bits for 1/2 of the constraints, but the row-sum lemma kills
     every w of large order.
   - 3-cycle twists (571, d = 15) do not decouple: the Z₃ components E₀ and E₁ are coupled pointwise.

4. **Z₃-compression for v = 3p.** The GS condition over Z₃ × Z_p splits exactly into:
   - a compression layer Y_k ∈ {±1, ±3} with ΣPSD = 12p;
   - an Eisenstein layer Z_k = −Y_k·ω^{j_k} on the support where |Y_k| = 1, with total support Σ|S_k| = 3p (Parseval) and
     ΣPSD = 3p.

   The compression layer is underdetermined (about 7.2p bits against 2.9p). The Eisenstein layer is overdetermined
   (about 4.75p bits against 5.4p). So a random compression lifts with probability about 2^(−0.65p) [E]. Fixing an
   algebraic Eisenstein layer instead (for example sextic characters) leaves only p free bits for about 2.9p constraint
   bits. Neither order of stages gives a cheap search.

5. **Conclusion for (C).** Among GS b = 0, b = 1 (cyclic and non-cyclic), LP, T-type, twisted and compression-lift
   ansätze, none reaches ≤ 2^45 at any of the eleven v. The repo's two near-misses are the only configurations within
   about 10 bits:
   - 523 flat-restricted b = 0 d = 9 (n′ ≈ 91), already running;
   - 571 twisted d = 15 (n ≈ 95, but lopsided).

   A genuine ≤ 2^45 route would need either new symmetry with feasible row sums at d ≈ 15–25, which does not exist for
   these v (§4.4), or an algebraic family.

---

## 6. Lemmas and proofs

**Lemma A (Williamson array).** Let W = [[A,B,C,D],[−B,A,−D,C],[−C,D,A,−B],[−D,−C,B,A]] with real n×n blocks, and write
Δ_XY = XYᵀ − YXᵀ.
- The diagonal blocks of WWᵀ all equal ΣXXᵀ.
- The first block row has off-diagonal blocks Δ_BA + Δ_DC, Δ_CA + Δ_BD and Δ_DA + Δ_CB. The other rows are signed
  permutations of these.

So W is Hadamard iff ΣXXᵀ = 4nI and Δ_AB = Δ_DC, Δ_AC = Δ_BD, Δ_AD = Δ_CB. Pairwise amicability is sufficient but
not necessary.

*Proof.* Write W = A⊗1 + B⊗λ(i) + C⊗λ(j) + D⊗λ(k) with λ the left regular representation of the quaternions. Then use
λ(q)ᵀ = λ(q̄) and the multiplication table. Checked numerically on random blocks (`theory/lemma_checks.py`). ∎

**Lemma A′.** If W is as in Lemma A, then W commutes with I⊗ρ(i) and I⊗ρ(j), the right multiplications. These are
anticommuting signed permutations with square −I. Hence Aut(W) contains a Q₈ whose centre is (−I, −I), and so does the
automorphism group of every matrix equivalent to W. Only this direction is used above. Conversely (sketch), a free Q₈ of
signed permutations with centre (−I, −I) splits the coordinates into orbits {±e, ±Je, ±Ke, ±JKe}. In that basis H takes
the form of Lemma A. ∎

**Lemma B (Cooper–Wallis plug-in).** Let T₁..T₄ be circulant T-matrices of order t and A..D commuting matrices with
ΣXXᵀ = 4wI. Put X_r = Σ_k T_k ⊗ N_rk, with N the Williamson pattern. Then

  Σ_r X_r X_rᵀ = 4tw·I + Σ_{k≠l} T_k T_lᵀ ⊗ G_kl,  where G_kl = Σ_r N_rk N_rlᵀ is antisymmetric.

Consequences:
- If the T_k are pairwise amicable, the cross terms cancel for **any** commuting A..D (GS quadruples included).
- **There is no amicable T(3)**: all 192 T-matrices of order 3 are signed shifts, and none is amicable [M].
- With T(3) = (I, S, S², 0), the cross term is (S⁻¹ − S) ⊗ (G₁₂ + G₂₃ + G₃₁). So one extra antisymmetric circulant
  identity in A..D replaces amicability. That is (w−1)/2 extra constraints on top of GS: about 4w free bits against
  roughly 5w constraint bits, so overdetermined [E].

Verified numerically (`lemma_checks.py`). ∎

**Also.** Any plug-in N = Σ_k P_k ⊗ A_k with A_k Williamson-type and ±1 entries forces the P_k to have disjoint supports.
NNᵀ = 4mn·I then requires P_k·P_kᵀ = mI and P_k·P_lᵀ + P_l·P_kᵀ = 0, so (P_k) is an OD(4m; m,m,m,m). In this form the
ingredient for v = 3p is a BH array of order 4p or 12, nothing weaker. ∎

**Lemma C (product condition vs amicability).** Let Q be a W-type quaternion matrix of odd order n > 1 (entries
a+bi+cj+dk with a, b, c, d = ±1). Then not all entries can have the same coordinate product abcd.

*Proof.* q/2 lies in the binary tetrahedral group, and its two half-integer Q₈-cosets are exactly abcd = ±1. If all
entries are in one coset, each term q_ij·q̄_kj is 4 times an element of {±1, ±i, ±j, ±k}. A sum of n such units vanishes
only if each axis count is even, which forces n even. ∎

Corollaries:
- T-type (the product condition) and W-type matrices are disjoint classes at odd order.
- The coset-preserving Kronecker multiplication of W-type matrices is vacuous at odd order.

**Lemma D** (Frobenius group): §3, item 2.

**Lemma E (row sums).**
- (E1) If an S-invariant block has fixed set F and all other orbits have sizes divisible by D, then
  a ≡ Σ_{x∈F} X(x) (mod D), and mod 2D when the number of such orbits is even. The exact feasibility test is the
  signed-orbit-size DP (`cost.py::rowsum_dp`).
- (E2) If every orbit has size divisible by 3, then 3 | a_k for all k. So 9 | Σa_k², which is impossible at b = 0 when
  9 ∤ 4v. This applies to v = 3p, p ≠ 3, i.e. 537, 573 and 669. ∎

**Lemma F (Legendre ⊕ D-optimal).**
- For prime v ≡ 3 (mod 4), PAF_L(s) = −1 + L(0)(χ(s) + χ(−s)) = −1, because χ(−1) = −1.
- So GS(L, L, C, D) is Hadamard iff PAF_C + PAF_D = 2 off 0. That is a circulant D-optimal pair of length v, which needs
  2v − 1 = a² + b².
- More generally, any LP(v) together with a D-optimal pair of length v gives H(4v). ∎

**Lemma G (support splits).** §5, item 2. The count 2^(B − F − n) per fixed support follows from the same Gaussian model
as n.

**Lemma H (order-2 twists at primes ≡ 3 mod 4).** If v is prime ≡ 3 (mod 4) and u has even order 2e, then u^e = −1, since
the 2-part of (Z/v)* is {±1}.
- Suppose u maps A ↦ A, D ↦ D and B ↔ C. Then the odd power u^e = −1 gives A∘(−1) = A, D∘(−1) = D and C = B∘(−1).
- Suppose instead u maps A ↔ D and B ↔ C. Then D = Ā and C = B̄, so ΣPAF = 2(PAF_A + PAF_B). That would make (A, B) a
  periodic Golay pair of odd length, impossible since PAF ≡ v (mod 4). ∎

**Quadratic generalised cyclotomy for ℓ = pq (sketch).** Take A = f + g, B = f′ + g′, with f = χ_pχ_q on units and g
constant (μ on the p-axis, ν on the q-axis, c at 0).
- If χ_p(−1)χ_q(−1) = −1, the cross-correlation of f with an axis-supported g vanishes. Otherwise it is
  f(s)(1 + f(−1))(c − μ − ν), and the same equations below follow once that term is matched.
- The axis-lag equations then become c_Aμ_A + c_Bμ_B = p − q and c_Aν_A + c_Bν_B = q − p, with all terms ±1. This
  forces |q − p| ≤ 2. With q − p = 2 the unit-lag equation μ_Aν_A + μ_Bν_B = −2 is consistent (μ = −c, ν = c), which
  is the twin-prime construction. Our ℓ have gaps of 8 or more.
- The full index-4 slice (9 orbits) is covered exhaustively by the LP sweep [M].

---

## 7. Computations: scripts, controls, results

All under `theory/`. Every "0 solutions" is an exhaustive result for the stated slice family.

| script | what it does | control |
|---|---|---|
| `abgroup.py`, `subgroups.py` | abelian groups G, Aut(G) = GL-parts × units, all subgroups up to GL-conjugacy (Goursat) | subgroup counts match hand counts: Aut(Z₅×Z₇) has 16; S₃×Z₂ up to conjugacy has 10; Aut(Z₈×Z₆₇) has 64 |
| `slicesolve.py` | exhaustive solver for S-invariant LP / b=1 / b=0 slices: numba Gray code, PSD and row-sum filters, exact integer PAF, exact hash join, memory caps that report "undecided" | `test_controls.py`: certified H(28), H(60), H(32), H(72), H(36), H(52), H(52), H(76); `control_340.py`: certified H(340) from b = 1 over Z₂²×Z₃×Z₇ with \|S\| = 6 |
| `certify.py` | builds the full matrix (two-circulant LP, bordered GS for any abelian G, GS) and checks H·Hᵀ = N·I | same |
| `sweep_b1.py`, `run_lp_all.sh`, `run_b0b1_cyclic.sh` | all slices with r ≤ 30 per group | – |
| `cost.py`, `cost_sweep.py` | Gaussian n and log₂E (exact uniform row-sum sampling) | matches the repo's model conventions |
| `compress_screen.py`, `screen_top.py` | exact quotient screens (budgeted enumeration, exact join) | `screen_control.py`: 34 screens on solvable slices, 0 false negatives |
| `hadaut_traces.c`, `run_aut.py`, `check_autgens.py` | Aut(H) via Traces | every element verified |
| `closure_table.py` | all arithmetic in §1.2 | – |
| `lemma_checks.py` | Lemmas A, B, C numerically and exhaustively | – |

**Building the automorphism tool.** nauty is not installed system-wide. Its source ships inside the PyPI sdist:

```
pip download pynauty==2.8.8.1 --no-deps --no-binary :all:
tar xzf pynauty-2.8.8.1.tar.gz
cd pynauty-2.8.8.1/src/nauty2_8_8 && ./configure && make nauty.a
gcc -O3 -o theory/hadaut_traces theory/hadaut_traces.c -I<that dir> <that dir>/nauty.a
python3 theory/run_aut.py 716 out.txt
python3 theory/check_autgens.py 716 out.txt
```

**Exhaustive sweep results** (from `summarize_sweeps.py`; r ≤ 30 orbits per sequence):

| kind | group G | slices (r<=30) | max r | blocks surviving PSD+row-sum (total) | solutions | undecided (join cap) |
|---|---|---|---|---|---|---|
| LP (ℓ = 2v−1) | Z5xZ11xZ19 | 38 | 30 | 504 | 0 | 0 |
| LP (ℓ = 2v−1) | Z29xZ37 | 33 | 30 | 942 | 0 | 0 |
| LP (ℓ = 2v−1) | Z7xZ163 | 40 | 30 | 16060 | 0 | 0 |
| LP (ℓ = 2v−1) | Z5xZ229 | 26 | 26 | 6992 | 0 | 0 |
| LP (ℓ = 2v−1) | Z3xZ3xZ7xZ19 | 474 | 30 | 14048 | 0 | 0 |
| LP (ℓ = 2v−1) | Z9xZ7xZ19 | 116 | 30 | 108 | 0 | 0 |
| LP (ℓ = 2v−1) | Z5xZ257 | 24 | 27 | 11020 | 0 | 0 |
| LP (ℓ = 2v−1) | Z3xZ431 | 9 | 22 | 2012 | 0 | 0 |
| LP (ℓ = 2v−1) | Z3xZ439 | 10 | 21 | 2860 | 0 | 0 |
| LP (ℓ = 2v−1) | Z7xZ191 | 15 | 24 | 1888 | 0 | 0 |
| LP (ℓ = 2v−1) | Z3xZ479 | 5 | 9 | 20 | 0 | 0 |
| LP (ℓ = 2v−1) | Z7xZ211 | 42 | 30 | 12202 | 0 | 0 |
| b = 0 (G = Z_v) | Z523 | 6 | 19 | 44376 | 0 | 0 |
| b = 1 (\|G\| = v−1) | Z2xZ9xZ29 | 10 | 30 | 125064 | 0 | 0 |
| b = 0 (G = Z_v) | Z3xZ179 | 5 | 9 | 40 | 0 | 0 |
| b = 1 (\|G\| = v−1) | Z8xZ67 | 27 | 30 | 898016 | 0 | 1 |
| b = 0 (G = Z_v) | Z571 | 8 | 20 | 63221 | 0 | 0 |
| b = 1 (\|G\| = v−1) | Z2xZ3xZ5xZ19 | 13 | 30 | 83380 | 0 | 0 |
| b = 0 (G = Z_v) | Z3xZ191 | 9 | 22 | 229088 | 0 | 0 |
| b = 1 (\|G\| = v−1) | Z4xZ11xZ13 | 23 | 30 | 31152 | 0 | 0 |
| b = 0 (G = Z_v) | Z599 | 4 | 24 | 639407 | 0 | 0 |
| b = 1 (\|G\| = v−1) | Z2xZ13xZ23 | 12 | 30 | 116300 | 0 | 1 |
| b = 0 (G = Z_v) | Z643 | 4 | 7 | 22 | 0 | 0 |
| b = 1 (\|G\| = v−1) | Z2xZ3xZ107 | 5 | 18 | 1560 | 0 | 0 |
| b = 0 (G = Z_v) | Z647 | 4 | 20 | 93254 | 0 | 0 |
| b = 1 (\|G\| = v−1) | Z2xZ17xZ19 | 13 | 30 | 170882 | 0 | 1 |
| b = 0 (G = Z_v) | Z659 | 4 | 15 | 3266 | 0 | 0 |
| b = 1 (\|G\| = v−1) | Z2xZ7xZ47 | 9 | 28 | 7644 | 0 | 0 |
| b = 0 (G = Z_v) | Z3xZ223 | 10 | 21 | 32468 | 0 | 0 |
| b = 1 (\|G\| = v−1) | Z4xZ167 | 5 | 12 | 132 | 0 | 0 |
| b = 0 (G = Z_v) | Z719 | 2 | 3 | 2 | 0 | 0 |
| b = 1 (\|G\| = v−1) | Z2xZ359 | 2 | 6 | 8 | 0 | 0 |
| b = 0 (G = Z_v) | Z739 | 6 | 19 | 43350 | 0 | 0 |
| b = 1 (\|G\| = v−1) | Z2xZ9xZ41 | 10 | 30 | 119626 | 0 | 0 |
| b = 1 (\|G\| = v−1) | Z2xZ2xZ2xZ67 | 80 | 30 | 292852 | 0 | 0 |
| b = 1 (\|G\| = v−1) | Z2xZ2xZ167 | 7 | 12 | 132 | 0 | 0 |
| b = 1 (\|G\| = v−1) | Z2xZ3xZ3xZ29 | 55 | 30 | 752816 | 0 | 2 |
| b = 1 (\|G\| = v−1) | Z4xZ2xZ67 | 7 | 30 | 16544 | 0 | 0 |
| b = 1 (\|G\| = v−1) | Z2xZ2xZ11xZ13 | 48 | 30 | 322271 | 0 | 1 |
| b = 1 (\|G\| = v−1) | Z2xZ3xZ3xZ41 | 65 | 30 | 577240 | 0 | 0 |

totals per kind: {'lp': 832, 'b0': 62, 'b1': 391}

Notes on the table:
- "blocks surviving" counts S-invariant sequences that pass the row-sum and PSD filters, before the exact join.
- Six slices exceeded the join cap ("undecided"). Five were then **decided** by exact quotient screens
  (`decide_undecided.py`, `screen_undecided_537.py`, logs in `theory/logs/`):
  - Z₈×Z₆₇, |S| = 22;
  - Z₂×Z₁₃×Z₂₃, |S| = 33;
  - Z₂×Z₁₇×Z₁₉, |S| = 36;
  - both Z₂×Z₃²×Z₂₉ slices with |S| = 28.
- **One slice stays undecided:** Z₂²×Z₁₁×Z₁₃ (573), |S| = 45, r = 30. It survives the Z₂², Z₁₁ and Z₂²×Z₁₁ quotients,
  and its Gaussian log₂E ≈ −58 [E].
- Z₄×Z₂×Z₆₇ was swept with diagonal (unit) automorphisms only; its full Aut(Z₂×Z₄) part was not enumerated.
- Totals: 832 LP slices, 62 b = 0 slices and 391 b = 1 slices, with no solutions in any of them.

---

## 8. Recommendations

1. **Keep 523 (b = 0, d = 9, flat k-tree) as the primary search.** Nothing here beats it. It is the only configuration
   with n′ < 100 and four independent lists.

2. **573: if search time becomes available, the best new target is b = 1 over Z₂²×Z₁₁×Z₁₃ with |S| = 15 (r = 54).**
   - Cost: n ≈ 167 [E], against 193 for the repo's cyclic version.
   - First computation: decide the three undecided quotient screens (Q = Z₁₃, Z₂²×Z₁₃, Z₁₁×Z₁₃). That needs a DP over
     S-orbits keyed on the compressed PAF, instead of the current enumeration with its 1.2M-vector cap. Then, only if the
     slice survives, a flat-pool measurement of n′ as in the repo's S1 step (`s1_nprime.py`), generalised to this G.
   - Positive control [M] (`theory/control_340.py`): b = 1 over Z₂²×Z₃×Z₇ (v = 85), the same Z₂²×Z_p×Z_q shape, with
     symmetry groups of order 6 (r = 20 and 21).
     - The solver finds solutions in three such slices, assembles H(340) and certifies H·Hᵀ = 340·I.
     - Every decidable quotient screen passes on them.
     - Saved in `theory/logs/controls/`.
     - The smaller Z₂²×Z₃×Z₅ (v = 61) has no solvable slice with |S| ≥ 2 and r ≤ 26, so it is not a usable control.
   - Expectation: at n ≈ 167 the bucketed k-tree costs about 2^113. The slice is worth anything only if a flat
     restriction brings n′ below about 100, a cut of about 40%; the repo's flat restriction cut 36% at 523. So treat this
     as a cheap screening question (two computations), not as a search to launch.

3. **739: no improvement.** The best non-cyclic b = 1 slice (|S| = 20, n ≈ 187) is killed by an exact quotient screen
   [M]. Keep the repo's b = 0 d = 9 configuration (about 210 bits) as the reference.

4. **Algebraic lottery, cheap and parallel.**
   - Extend the exhaustive LP and b = 1 sweeps to r ≤ 34. Each slice is about 2^33 Gray-code steps, a few minutes in
     numba, and the families have tens of slices per group.
   - Add block-twisted b = 1 slices (signed block permutations combined with Aut(G)). These are not covered yet.
   - The prior is low, but any hit closes an order with a certificate.

5. **Do not spend compute on:**
   - Williamson-type ingredients from the 2026 matrices (§2);
   - Frobenius / cocyclic structure at 669 (§3);
   - T(3)-based relaxations (Lemma B);
   - Z₃-compression lifts (§5, item 4);
   - support or symmetric splits of b = 0 or b = 1 (Lemma G);
   - LP or D-optimal pairs at random (overdetermined, §4).

6. **What would change the picture:**
   - an algebraic family for 4-SDS over Z_p with p ≡ 3 (mod 4) beyond Paley-type defect −4 (the b = 1 case);
   - an M-structure theorem whose hypothesis the existing GS-type H(4p) satisfy (check Seberry–Yamada 1990 directly);
   - a Williamson-type matrix of order 179, 191 or 223 found by other means: each gives an open v = 3p immediately via
     Cooper–Wallis with T(3).
