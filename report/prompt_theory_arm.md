# Prompt: theory arm for the open Hadamard orders 2000–3000 (self-contained)

*For a maximum-reasoning session, with or without tools. If you have Claude Code, the repository is
`middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`, and the cited files exist there. Tags:
[M] verified in exact arithmetic, [E] estimate, [L] literature recollection to be checked.*

## 0. Goal

Be the first to construct a Hadamard matrix of order 4v for at least one of

    v ∈ {523, 537, 571, 573, 599, 643, 647, 659, 669, 719, 739}     (orders 2092 … 2956).

Constructing a single matrix is enough, and a matrix certifies itself (H·Hᵀ = 4v·I in integers).

**What I want from you is not more search.** It is one of:
- **(A)** an algebraic family or explicit construction covering at least one of these v;
- **(B)** a composition (plug-in, multiplication or array theorem) that reaches one of them from ingredients that
  exist now;
- **(C)** a structured ansatz whose effective search cost is at most 2^45 elementary operations, with a positive
  control.

Competitors are active. In 2026 one group closed every open order below 2000, and an independent construction closed
2060 (§2). A theorem covering a whole class of v is the kind of result a large-lab sweep would produce, and it would
be the best outcome here too.

## 1. Arithmetic of the targets

Every open v is squarefree.

| v | factorization | v mod 8 | 2v − 1 | v + 1 |
|---|---|---|---|---|
| 523 | prime | 3 | 5·11·19 | 2²·131 |
| 537 | 3·179 (179 ≡ 2 mod 3) | 1 | 29·37 | 2·269 |
| 571 | prime | 3 | 7·163 | 2²·11·13 |
| 573 | 3·191 (191 ≡ 2 mod 3) | 5 | 5·229 | 2·7·41 |
| 599 | prime | 7 | 3²·7·19 | 2³·3·5² |
| 643 | prime | 3 | 5·257 | 2²·7·23 |
| 647 | prime | 7 | 3·431 | 2³·3⁴ |
| 659 | prime | 3 | 3·439 | 2²·3·5·11 |
| 669 | 3·223 (223 ≡ 1 mod 3) | 5 | 7·191 | 2·5·67 |
| 719 | prime | 7 | 3·479 | 2⁴·3²·5 |
| 739 | prime | 3 | 7·211 | 2²·5·37 |

- The only group of order v is cyclic, **except 669**, which also admits the non-abelian Frobenius group
  Z₂₂₃ ⋊ Z₃.
- None of the 2v − 1 is a prime power. None is of the form 2q + 1 with q a prime power, or p(p + 2) with p, p + 2 twin
  primes. None of the v is of the form (q + 1)/2 with q a prime power.

## 2. Ingredient library as of October 2026

**Hadamard matrices.** Every order 4n < 2000 is now known. The twelve formerly open orders were closed by Alpöge,
Voinov and Reynolds-Haertle (August 2026): 668, 716, 892, 1132, 1244, 1388, 1436, 1676, 1772, 1916, 1948 and 1964.
- Their cores (in `data/core_v*.txt`) are Goethals–Seidel (GS) quadruples of three kinds:
  - **b = 0:** four ±1 circulants with ΣPAF = 0, at v = 223, 283, 311, 487, 491;
  - **b = 1:** one border row, at v = 167, 179, 419, 443;
  - **b ≥ 3:** an "inflated" bordered architecture, at v = 347, 359, 479. See §4 and `report/bge3_construction.md`.

**T-matrices of order 103.** These are K-invariant for K = {1, 46, 56} ≤ (Z/103)*. They come from L. Schneider's
H(2060) (schneiderlo/hadamard-2060), built by Cooper–Wallis with T(103) and Williamson matrices W(5).
- We found a second T(103) blind (`cw/T103_blind.json`).
- We used Cooper–Wallis with W(7) to build **H(2884)**, v = 721 [M], `cw/RESULT_2884.md`.

**Other ingredients:**
- The classical libraries [L]: Williamson and Williamson-type orders, T-sequences and base sequences BS(n+1, n) known
  for small n, Golay and Yang numbers, orthogonal designs, and Baumert–Hall, Welch and Ono–Sawade–Yamamoto arrays.
- **Circulant Hadamard conjecture:** reportedly proved (OpenAI, Sept 2026). This is irrelevant here, since nothing
  needs a circulant Hadamard matrix.

## 3. Closed or priced-out routes (do not re-derive; challenge them if you think one is wrong)

| route | status | reason |
|---|---|---|
| T-matrices + Cooper–Wallis for 573 = 3·191, 537 = 3·179 | dead | T-form row sums: three ≡ 0, one ≡ ±1 mod k, Σs² = t, impossible for t = 191 with every odd k; (Z/179)* has no usable odd subgroup [M] |
| T(223) for 669 | 2^157 density | about 94 bits beyond the solved T(103); twisted symmetries (multiplier plus signed row permutation) give no gain [M] |
| H(892) as a Williamson-type ingredient for 669 | no | no pairwise cross-correlation of the 223 core is symmetric under any shift or reversal [M] |
| H(716) as a Williamson-type ingredient for 537 | **untested** | would need an automorphism of order 179 (nauty unavailable) |
| T-matrices for prime v | no gain | equivalent to a b = 0 GS core with the product condition a·b·c·d = +1 |
| b = 0 GS core at 523, d = 9 | **live, being searched** (Kaggle k-tree, flat pools) | n ≈ 142 bits (flat-restricted ≈ 91) |
| twisted b = 0 at 571, d = 15 (block 3-cycle) | feasible at 94.6 bits but impractical | reduces to a lopsided 2-list problem, 2^35 × 2^110 |
| original b ≥ 3 architecture, all open v | ≥ 2^80 or dead | exact screens kill every option at or below 300 bits for 599, 647, 739; best survivor 571/c = 22 needs ≳ 2^80 and a border matched to the compression (none known for c ≥ 10) [M] |
| generalized-Ψ bordered architecture (Ψ ⊂ Q̂, \|Ψ\| = b odd) | best 523/c = 12 at about 2^62, then killed by the exact screen [M]; others ≥ 2^80 | `report/bge3_answer_review.md` |
| algebraic slices with large multiplier symmetry (182 exhaustive slices) | empty | `report/bge3_answer.md` §3.3 |

**Our practical compute limit.**
- Lists of about 2^28 per slot, i.e. 2^30–2^40 operations per hit.
- The k-tree 4-list join has solved instances at about 63 bits (T(103)) and about 83 bits flat-restricted (v = 487).
- Local search scales as about 2^(0.55·n).

## 4. The architectures, briefly (exact statements)

- **GS b = 0:** four ±1 functions on a group G of order v with Σ_k PAF_k(g) = 0 for g ≠ 0. They plug into the
  Goethals–Seidel array, giving H(4v).
- **Bordered b ≥ 3:** v = c(h + 1) − 1 and G = Q × Z_h with |Q| = c = b + 1. The core satisfies
  Σ_k PAF_k = −4b on H ∖ {0} and +4 off H. The border is an H(4c) with the four rows H₄ ⊗ 1_c removed, inflated by h.
  The corner is forced: C = −K·Cbᵀ·Kcᵀ/(4c) must be ±1.
- **Generalized-Ψ border:** the target is Σ_k PAF_k(q, z) = 4v·δ − 4·Σ_{ψ∈Ψ} ψ(q).

## 5. What to think about (ranked by my guess at upside; reorder freely)

1. **Composition closure, done rigorously.** Enumerate every plug-in or multiplication theorem you know that outputs
   H(4v) with v odd, from inputs that exist as of October 2026. For each open v say which inputs it would need, and
   whether they exist.
   - In particular: T-matrices, T-sequences or orthogonal designs OD(4t; t, t, t, t) of order t | v; Williamson-type
     matrices of order w | v; Baumert–Hall or Welch arrays.
   - Also the newly available ingredients: T(103), H(2060), H(2884), and the 2026 cores as possible Williamson-type,
     T-type or SDS-type objects after equivalence.
   - The inference "v is open, therefore ingredient Y is unknown" expired for every order closed in 2026. Re-check it
     for each such inference.
2. **669 = 3·223 and the non-abelian group Z₂₂₃ ⋊ Z₃.**
   - Are there GS-type or Ito/cocyclic constructions over non-abelian groups, for example a group-developed core or a
     difference family in a Frobenius group, that give H(4·669)?
   - Is the Z₃ action usable as an extra symmetry that cuts the search bits for a structured core?
   - A cyclotomic construction in GF(223) with the 3-fold structure would be ideal.
3. **Algebraic families matching an open v.** These include:
   - relative or divisible difference families;
   - cyclotomic constructions over GF(q) × Z_m;
   - Legendre pairs or two-circulant cores of length ℓ = 2v − 1. These are composite for every open v, so they admit
     compression and subgroup structure.

   For each, state the exact construction, which v it hits, and the theorem or computation it relies on. Flag anything
   that is conjecture.
4. **A structured ansatz at ≤ 2^45**, with an argument for why its constraint count is low. For example, a symmetry
   group larger than multipliers, or exact factorization of the problem into independent lists. The T(103) success
   came from a support split that turned coupled sequences into four independent lists. A similar "decoupling" for
   any architecture at an open v would be valuable.
   - Include a positive control: a smaller instance of the same ansatz that must have solutions.

## 6. Rules for the answer

- **Certificates.** A construction counts only if it can be verified by exact computation (H·Hᵀ = 4v·I), or follows
  from a cited theorem whose hypotheses you check explicitly. "Expected to exist by a heuristic" is not a construction.
- **Separate** what is a theorem (with citation), what is new and proved by you (give the proof), and what is a guess.
- **For every proposed route, give:**
  - the exact objects to construct;
  - their size or bit cost;
  - the smallest checkable instance;
  - the first computation to run.
- **A useful first milestone:** a complete, rigorous composition-closure table for the eleven open v, listing the
  ingredients each route needs and their existence status as of October 2026. Even if no route closes, this table is
  valuable.
