# Alpöge-style builders (b = 0 and b ≥ 3) with the staged method

Date: 2026-10-14. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_alpoge_design.md`.
Scripts: `theory/ab_*.py`. Logs: `logs/ab/`. Bundle: `pub/H3732/`.

**Tags.** [M] exact computation (script named); [P] proved here (proof in the text); [E] Gaussian estimate (estimator
and calibration named); [L] literature recollection, unchecked.

---

## 0. Bottom line

1. **(A) H(3732) is built: v = 933, b = 0 (Goethals–Seidel, array `g`) [M].**
   - Four ±1 circulants of length 933 with Σ_k PAF_k = 3732·δ. Over Z₉₃₃ = Z₃ × Z₃₁₁, each sequence is constant on the
     sets {g} × C_j, where the C_j = 17^j·M are the ten cosets of the index-10 subgroup M = ⟨260⟩ (order 31) of
     F₃₁₁^*. That is 33 sign orbits per block.
   - Found by the staged search `theory/ab_b0z3.py`:
     - stage 1 on Z₃ (36 exact solutions);
     - a dual level on Z₃₁₁ (the compression along Z₃);
     - an exact cross join over the fibres of each dual solution.
   - The first core came after 173 dual solutions: about 2^25.7 hashed join entries and about 3 minutes on one
     core, after 20 s of setup.
   - H Hᵀ = 3732·I was checked three ways:
     - by the search (float64 BLAS);
     - by the independent `scripts/verify_h3732.py` (integer arithmetic, no shared code);
     - by the standalone `pub/H3732/verify.py`, which also checks the published normalized matrix.
   - sha256 of the packed rows: `8fea9ea2e581cdaa76dbac0af7779d9bdefac8edb2eda932aa29eae5da1a437e`.
   - 30 cores in all (10 orbits under a Z₃ shift) came from two stage-1 solutions in 25 minutes (§4).
2. **The method is validated end to end at the same shape [M], and its Gaussian is calibrated against exact counts.**
   - The same code builds H(372) (v = 93 = 3·31) and H(852) (v = 213 = 3·71), both verified by the independent
     verifier.
   - The cross-level Gaussian (lattice-exact, `theory/ab_b0z3_calib.py`) is within about 1 bit of the exact core
     counts, in both directions, at q = 71 and at q = 311 (§4.5).
3. **b ≥ 3: the staged pipeline reproduces the published 347 slice end to end (positive control) [M].**
   - In 347's own slice (b = 5, c = 6, G = Z₂ × Z₃ × Z₃ × Z₁₉, index 6 on Z₁₉) the hierarchical stage-2 solver
     (S/D split, Wagner middle level, compiled D-check) found a **new** core:
     - it differs from the published core in 540 of 1368 entries;
     - it is completed with the transferred border and forced corner to a new H(1388), verified exactly
       (`scripts/verify_packed_h.py`).
   - Cost: 11 middle values, 2^14.2 S-solutions, 28 minutes on one core, well under the 2^40 asked for. The Gaussian
     D-deficit (16.2 bits [E]) predicted about 2^16 S-solutions per core.
4. **Borders [P, M].**
   - c = 4, 6, 8: each published frame borders every open order with the same c, by the transfer theorem (§3).
   - **Every even c with a known H(2c)**: a frame exists, inside H₂ ⊗ H(2c) (§3.1) [P].
   - What remains for c ≥ 10 is not the frame but the forced corner. The corner is ±1 exactly when the compressed
     core, restricted to the block-balanced space, has a flat matrix in the two frames (§3.1).
5. **(C) Exact obstructions [M]** (§7). Among them:
   - b = 0:
     - these GP slices are empty: 537, 669 and 789 at every s ≤ 6; 933 at s = 2, 5; 917 (q = 131) at s = 2, 5, and at
       s = 10 with any multiplier on Z₇; 955 (q = 191) at s = 2, 5, 10;
     - 955 at s = 19 with −1 on Z₅, and 917 at s = 13 with −1 on Z₇, have no stage-1 solution.
   - b ≥ 3:
     - every GP slice of 571 at c = 4 (q = 71, s = 2, 5, 7, 10, 14) dies at the Q-level;
     - 863, c = 4, s = 6, U = −1: exact stage 1 is empty (31,188 α-quadruples);
     - 863, c = 6, in 347's border class: empty;
     - 599, c = 6, q = 11, s = 2, H' = Z₉, with five multiplier groups (§7): empty.
6. **(D)** The ranked table is §8.
   - The other composite orders reach b = 0 GP slices whose stage 2 is Gaussian-positive, but whose fibres are too
     large for the 933 pipeline (917 at s = 26 with an order-3 multiplier: +52 bits per y over 38,380 y's; 955 at
     s = 38: +113 bits per y). They need one more level (§6).
   - The best b ≥ 3 survivor is 659 at c = 4 (q = 41, s = 8), where stage 2 has C₂ ≈ 285 bits (§7).

---

## 1. The staged systems [P], checked on the published cores [M]

### 1.1 Setting

- b ≥ 3: c = b + 1 even, c | v + 1, h = (v + 1)/c − 1, G abelian of order m = ch, H ≤ G of index c, Q = G/H.
  Core: four ±1 functions x_k on G with Σ_k PAF_k = T, T(0) = 4m, T = −4b on H \ 0, T = +4 off H.
  Character form: Σ_k |X̂_k(χ)|² = 4v, except 4b when χ ≠ 1 is trivial on H.
- GP-type slice: G = G' × Z_q with q a prime, Z_q ≤ H (automatic when q ∤ c), so H = H' × Z_q, H' ≤ G' of index c.
  The slice group is S = U × M: M ≤ F_q^* of index s (order f = (q − 1)/s) acts on Z_q, U is a group of
  per-axis multipliers on G' preserving H'. On the fibre over b ∈ G': x(b, 0) = α(b), x(b, z) = c_j(b) on the class
  C_j = γ^j M.
- b = 0 at composite v = F·q: the same with G' = Z_F and T = 4v·δ.

### 1.2 Stage 1 [P]

Compressing along Z_q: y_k(b) = Σ_z x_k(b, z) = α_k(b) + f·σ_k(b), σ = Σ_j c_j, and

    Σ_k PAF_{G'}(y_k)(b) = T1(b) := Σ_z T(b, z) = 4m − 4b(q − 1) (b = 0), −4bq (b ∈ H' \ 0), +4q (b ∉ H').

Proof: (b, z) ∈ H ⟺ b ∈ H' because Z_q ≤ H; sum T over the fibre. In characters: Σ_k |ŷ_k(φ)|² = 4v, except 4b
when φ ≠ 1 is trivial on H'. So stage 1 is itself a "b ≥ 3 core" on the smaller group G' with integer letters.

Consequences used by the screens:
- **Lag 0 [P].** Σ_{k,b}(y² − 1) = 4fs(|G'| − b), i.e. with y = α + 2ft (s even): Σ (αt + ft²) = s(|G'| − b)
  (b = 0: s·F; b = 1: s(F − 1)).
- **Q-level [P].** The Q-compression z_k(q) = Σ_{b ∈ q+H'} y_k(b) satisfies Σ_k PAF_Q(z_k) = T_Q, i.e.
  Σ_k |ẑ_k(ψ)|² = 4b (ψ ≠ 1), Σ_k a_k² = 4v; hence every block has |ẑ_k(ψ)|² ≤ 4b.
- **α-congruence [P].** Every letter is ≡ ±1 (mod 2f₀) for the largest such f₀; with y = α + 2f₀t,
  PAF(y) = PAF(α) + 2f₀·(CC(α, t) + CC(t, α)) + 4f₀²·PAF(t) and the middle term is 2f₀ times an even number, so
  Σ_k PAF(α_k) ≡ T1 (mod 4f₀) at every lag (the b ≥ 3 targets enter only through T1), rows of α ≡ rows of y (mod 2f₀).
- **S/D [P].** On G' = Z₂ × G'': σ = (y(0, ·) + y(1, ·))/2 and δ = (y(0, ·) − y(1, ·))/2 satisfy
  Σ_k PAF(σ_k) = (T1(0, ·) + T1(1, ·))/4 and Σ_k PAF(δ_k) = (T1(0, ·) − T1(1, ·))/4.

### 1.3 Stage 2 [P]

For ψ = ψ_a ≠ 1 on Z_q, X̂_k(φ, ψ) = Σ_b φ(b) Z_{k,a}(b), Z_{k,a}(b) = α_k(b) + Σ_j c_{k,j}(b) η_j(ψ_a) (Gauss
periods). (φ, ψ) is nontrivial on Z_q ≤ H, hence on H, so Σ_k |X̂_k(φ, ψ)|² = 4v for every φ, i.e.

    Σ_k Σ_b Z_{k,a}(b) · conj(Z_{k,a}(b + b')) = 4v·δ_{b',0}    (in O_K, K = Q(ζ_q)^M),

exactly as for b = 0 and b = 1. Stage 1 + stage 2 ⟺ the core equations (Fourier inversion on G). **So the b = 1
stage-2 machinery (S/D split on a Z₂ factor of G', middle levels, D-check) applies unchanged; only the group G' and
its size differ.** In PAF form (used by the solver): with x_S = x(0, ·) + x(1, ·) and x_D = x(0, ·) − x(1, ·) on
G'' × Z_q, Σ_k PAF(x_{S,k}) = T(0, ·) + T(1, ·) and Σ_k PAF(x_{D,k}) = T(0, ·) − T(1, ·) [P].

### 1.4 Checks on the published cores [M] (`theory/ab_verify_published.py`, `logs/ab/verify_published.jsonl`)

| core | staging prime | slice found in the core | stage 1 exact | letters | α-congruence (f₀) | Q-compression | stage 2 (max error) |
|---|---|---|---|---|---|---|---|
| 479 (b = 3) | 7 | S = ⟨×13⟩: −1 on Z₇, order 4 on Z₁₇ | yes | ±1, ±3, …, ±7 | yes (2) | yes, ΣΣz² = 488 | 9·10⁻¹³ |
| 479 | 17 | same | yes | α + 4σ / α + 2·(8 signs) | yes (2) | yes | 9·10⁻¹³ |
| 347 (b = 5) | 19 | ⟨−1 on both Z₃⟩ × ⟨×7⟩ (index 6 on Z₁₉) | yes | α + 3σ | yes (3) | yes, 248 | 7·10⁻¹³ |
| 347 | 3 | same | yes | ±1, ±3 | yes (2) | yes | 7·10⁻¹³ |
| 359 (b = 7) | 11 | ⟨×9⟩: squares on Z₁₁ (Paley, s = 2) | yes | ±1, ±9, ±11 | yes (5) | yes, 204 | 5·10⁻¹³ |
| 892 (b = 0) | 223 | cube roots of unity, G' trivial | rows² = 892 | α + 3σ | – | – | 9·10⁻¹³ |

All three published b ≥ 3 cores are multiplier slices that act on a prime factor of h by a cyclotomic subgroup:
479 is doubly staged with a diagonal multiplier (like 647 ×183), 347 is a GP slice with U = −1 on the Z₃'s, and 359
is a pure Paley(11) GP slice.

---

## 2. Exact screens and their validation [M]

(`theory/ab_sweep.py`, `theory/ab_exact_small.py`, `theory/ab_s1sd.py`)

- **Screen A (Q-level, exact).** Per block, the Q-compression z and the lag-0 energy E over all letter assignments
  (DP per U-orbit of Q, energies as bitsets), the per-character bound |ẑ(ψ)|² ≤ 4b, and a hashed 4-block join of
  (PAF_Q(z), E) against (T_Q, T1(0)). For b = 0 it is lag 0 + rows.
- **Screen B (α-congruence, exact),** with row residues from the row quadruples admitted by A.
- **Exact stage 1** (`theory/ab_s1sd.py`): α-quadruples → S-join over u = t(0, ·) + t(1, ·) → D-join over
  w = t(0, ·) − t(1, ·) → verification; energy-split hashed 4-list joins.
- **Validation.**
  - Brute force on 100 small slices (79 decided, 22 with solutions, up to 8·10¹⁰ quadruples): **0 false negatives**
    for A and B (`logs/ab/exact_small.jsonl`).
  - The published 347 slice passes A and B; its exact stage 1 contains 347's own y (the y's α-quadruple solved alone
    returns 8 stage-1 solutions including it) [M].

---

## 3. Borders for c = 4, 6, 8 at every open order [P, M]

(`theory/ab_assemble.py`, `theory/ab_template.py`, `theory/ab_border_transfer2.py`, `logs/ab/border_transfer.json`)

- **Templates.** The three published matrices were decomposed into core functions F_k, layout, border K, Kc and
  corner C (G = Q × H split for 347 and 479; 359 is non-split with Q = G/(2Z₈ × Z₁₁) = Z₂³ and a real-character row
  twist). Each is rebuilt exactly by the assembler, and the forced-corner formula C = −K·Cb(z)ᵀ·Kcᵀ/(4c) reproduces
  each published corner [M].
- **Theorem (border transfer) [P].** Cb(z) is linear in z, and every row of K and Kc is balanced inside each block, so
  the block means do not enter the corner. For each template the kernel of z ↦ K·Cb(z)ᵀ·Kcᵀ is, block by block,
  one-dimensional (w = 1 for 347; a real character of Q for 479 and 359, whose templates are character twists of the
  natural cores) [M]. Moving block k's compression by s_k·w keeps the corner and moves the natural row sum a_k by
  c·s_k. So the template's border fits order v' (same c) whenever Σ_k (a_k + c s_k)² = 4v' has a solution with
  s_k ≡ h' − h (mod 2).
- **Result [M].** For every open order with c = 4 (18 orders), c = 6 (599, 659, 719, 839, 863, 917, 947, 971) and
  c = 8 (599, 719, 751, 823, 839, 863, 919, 991) the row equation has 18–106 solutions, and every one checked
  (up to 50 per order) gives a ±1 corner with C Cᵀ = 4b·I.
- **What this means for a search.** Restrict stage 0 (the Q-compression) to a template's class. The direct test is
  cheap (`theory/ab_borderfit.py`): a found y is bordered by a template frame if the forced corner is ±1 for some block
  order and signs; it accepts more than the shift class (e.g. 347's block B translated inside Z₃ also fits).
- **c ≥ 10:** no published template. Frames exist for every even c (§3.1); the open part is the flatness of the
  forced corner.

### 3.1 Frames for every even c, and what the corner really requires [P]

*Frame* = a 4b × 4c ±1 matrix K with K Kᵀ = 4c·I whose rows sum to 0 on each of the four blocks of c columns.
Equivalently, an H(4c) containing H₄ ⊗ 1_cᵀ as four of its rows: the other rows are orthogonal to those four, and H₄
is invertible.

**Lemma 1 (frames exist).** Let c be even with an H(2c), normalized so that its first row r is all ones, and let s be
any other row. In H₂ ⊗ H(2c), the four rows (r, r), (s, s), (r, −r), (s, −s) take two kinds of values:
- in column j of the first half: (1, s_j, 1, s_j);
- in column j of the second half: (1, s_j, −1, −s_j).

s is balanced, so each of the four patterns occurs exactly c times. Grouped by pattern, the four rows are
H₄ ⊗ 1_cᵀ, and the remaining 4c − 4 rows form a frame. H(2c) is known for every multiple of 4 below 668 [L] (668 is the smallest open order), so frames exist for
every even c ≤ 332. Checked numerically for c = 6 and c = 10 (Paley I H(12) and H(20)) [M]. ∎

**Lemma 2 (the corner is always orthogonal; only its flatness is in question).**
- Let V ⊂ R^{4c} be the block-balanced subspace (dimension 4b), Π_V the projection onto it, and K, Kc frames, so
  that Kᵀ K = Kcᵀ Kc = 4c·Π_V.
- The compressed GS array Cb of the Q-compressions z_k satisfies Cb Cbᵀ = I₄ ⊗ M(Σ_k PAF_Q(z_k)), M(·) the group
  matrix (the GS identity over Q). Its eigenvalues are Σ_k |ẑ_k(ψ)|² = 4b on V (ψ ≠ 1, the Q-level equation of §1.2)
  and 4v on V^⊥.
- Cb and Cbᵀ map block-constant vectors to block-constant vectors, so they preserve V.
- The same holds for Cbᵀ Cb (the transposed array is again of GS type).
- Hence the forced corner C = −K Cbᵀ Kcᵀ/(4c) (the convention of §3 and `theory/ab_assemble.py`) satisfies
  C Cᵀ = K Cbᵀ (4c·Π_V) Cb Kᵀ/(16c²) = 4b·K Kᵀ/(4c) = 4b·I
  for every Q-level solution z and every pair of frames.

So C is always a real orthogonal matrix times √(4b). It is a ±1 matrix exactly when all its entries have absolute value
1, i.e. when the matrix of Cb|_V in the two frame bases is flat. Equivalently, for a chosen Kc and a chosen corner
C ∈ H(4b), the other frame is forced to be K = −C·Kc·Cb/(4b), and the condition is that this K be ±1. ∎

Scope: the untwisted layout (blocks ±group-developed or ±back-developed on Q), as in 347 and 479. Every identity
above is checked on the 347 template [M]:
- K Kᵀ = Kc Kcᵀ = 24·I, with block-balanced rows;
- Cbᵀ Cb = Cb Cbᵀ = 20 on V;
- C Cᵀ = 20·I;
- K = −C·Kc·Cb/20 exactly.

359's character twist multiplies Cb by diagonal ±1 matrices and changes V accordingly.

**Consequences.**
- c ≥ 10 is not blocked by frames. The open question is the flatness of the corner for some Q-level compression z
  and some pair of frames.
- The published templates are flat solutions for c = 4, 6, 8, and §3's transfer theorem carries them to every open
  order with the same c.
- For c = 10 (the only b ≥ 3 option of 669 and 789) the first computation is a search for a flat pair: the Q-level
  solutions z on Z₁₀ or Z₂ × Z₅, against frames from H₂ ⊗ H(20) and their signed column permutations within blocks
  (§9).

---

## 4. H(3732): b = 0 at v = 933 = 3·311 [P, M]

(`theory/ab_b0z3.py`, `theory/ab_b0_stage1.py`, `theory/ab_b0_dual.py`, `theory/ab_b0z3_calib.py`; logs
`logs/ab/b0z3/`; bundle `pub/H3732/`; verifiers `scripts/verify_h3732.py`, `pub/H3732/verify.py`)

### 4.1 Slice

- G = Z₃ × Z₃₁₁ ≅ Z₉₃₃ by the CRT, n ↦ (n mod 3, n mod 311).
- M = ⟨260⟩ = ⟨17¹⁰⟩ is the index-10 subgroup of F₃₁₁^* (order f = 31, odd, so −1 ∉ M); C_j = 17^j·M, j = 0..9.
- Per block k: x_k(g, 0) = α_k(g) and x_k(g, z) = c_{k,j}(g) for z ∈ C_j. That is 3 + 30 = 33 signs per block, 132 in
  all. This is the GP slice "U = 1, s = 10" of the b = 0 sweep; it passes screens A and B [M].

### 4.2 Levels (exact identities) [P]

- **Stage 1** (compression along Z₃₁₁).
  - y_k(g) = α_k(g) + 31·σ_k(g), with σ_k(g) = Σ_j c_{k,j}(g) even, |σ| ≤ 10.
  - Σ_k PAF_{Z₃}(y_k) = 3732·δ.
  - Exact enumeration gives **36 solutions**, up to block order and sign [M] (`ab_b0_stage1.py`).
- **Dual level** (compression along Z₃).
  - x̄_k = A_k at 0 and T_{k,j} on C_j, with A_k = Σ_g α_k(g) ∈ {±1, ±3} and T_{k,j} = Σ_g c_{k,j}(g) ∈ {±1, ±3}.
  - Σ_k PAF_{Z₃₁₁}(x̄_k) = 3732·δ: six integer equations, at lag 0 and at the five lag classes of ⟨M, −1⟩.
  - Spectral bound per block: |A + Σ_j T_j η_j(ψ)|² ≤ 3732 for the ten Gauss-period vectors.
- **Energy lemma.** Σ_k A_k² + 31·Σ_{k,j} T_{k,j}² = 3732, with A² ∈ {1, 9} and Σ T² = 40 + 8·n₃ (n₃ = number of
  columns with |T| = 3, i.e. constant columns). Since 3732 ≡ 12 (mod 31) and 2492 = 10·248 + 12, this forces
  Σ_k A_k² = 12 and **exactly ten constant columns in all four blocks together**.
- **Fibres.** Given y and T_k, block k's signs form a 3 × 10 sign matrix with row sums σ_k(g) and column sums
  T_{k,j}. Their number m(T) ≤ 3^{10 − n₃} (at most 4,653 here).
- **Cross level.** G has 17 lag classes under ⟨M, −1⟩:
  - lag 0;
  - (±1, 0);
  - (0, C), five classes;
  - (1, C_j), ten classes, since −1 ∉ M.
  
  Given stage 1 and the dual level, 10 independent equations remain. Stage 1 + dual + cross ⟺ Σ_k PAF_G(x_k) =
  3732·δ (all 17 classes are checked).

### 4.3 Search [M]

1. **Dual lists.** All column-sum vectors T with multiplicities m(T), by a DP over columns, under the spectral bound:
   about 75,000 per block.
2. **Selection.** Per block and energy class, keep the K = 200 vectors of largest m(T): 1,802–1,870 per block.
3. **Dual join.** All pairs (0, 1) against all pairs (2, 3), with a 64-bit linear hash of the 6-vector of dual PAF
   values and an exact re-check of every match. It gives **64,894 dual solutions** in 2.3 s; their fibre products go
   up to 2^36.4, median 2^34.6.
4. **Order.** Dual solutions are taken in decreasing order of log₂(m₀m₁m₂m₃) − log₂(join cost + overhead). The yield
   is proportional to the fibre product and the cost to the best 2 + 2 split. So balanced fibres (36, 36, 4242, 4653)
   come before the larger but lopsided (1, 4242, 4653, 4653): their cost is 3.2·10⁵ instead of 2·10⁷ join entries.
5. **Cross join.**
   - The four fibres are listed once and cached, with the full PAF over G at the 17 lag classes, hashed.
   - 2 + 2 join, with the smaller side sorted and the larger streamed.
   - Every match is re-checked exactly (Σ_k PAF_G = 3732·δ).
6. **Assembly.** Each core is mapped to Z₉₃₃ by the CRT, re-checked cyclically, written as ± sequences, and the GS
   array is checked (H Hᵀ = 3732·I).

### 4.4 Same-shape controls [M]

The same code at q = 31 and q = 71 has the same combinatorics: 3 × 10 sign matrices, five ⟨M, −1⟩ classes, f odd.

| v | q | f | stage-1 solutions | dual solutions (K) | result |
|---|---|---|---|---|---|
| 93 | 31 | 3 | 110,266 | 14,178 (50) | cores at once; H(372) verified |
| 213 | 71 | 7 | 1,008 | 678,523 (200) | 570 cores from 300 dual solutions (41 s); H(852) verified |
| 933 | 311 | 31 | 36 | 64,894 (200) | 12 cores from 3,520 dual solutions (y₀), 18 from 2,598 (y₁); **H(3732)** |

`scripts/verify_h3732.py` checks H(372) and H(852) independently (integer H Hᵀ).

### 4.5 Calibration of the cross level [E vs M]

Estimator: the lattice-exact Gaussian (`theory/q523_gauss.gauss_bits`), fitted per dual solution to 20,000 random
4-tuples of fibre elements (`theory/ab_b0z3_calib.py`). The lattice has rank 9 at every dual solution tested.

| q | dual solutions | predicted cores (Gaussian) | exact cores | ratio |
|---|---|---|---|---|
| 71 | top 6 by fibre product (2^36.0) | 6 × 2^1.65 = 18.8 | 6 | 0.32 (−1.6 bits) |
| 71 | 10 by fibre product, first run | 31 | 18 | 0.58 (−0.8 bits) |
| 311 | 6,118 (y₀ and y₁, efficiency order) | 14.9 (43.2 bits per tuple) | 30 | 2.0 (+1.0 bit) |

- The cores come in orbits of 3. The block with y = (1, 1, 1) can be shifted along Z₃: its fibre is closed under the
  shift and its PAF is unchanged.
- So the event count is a third of the core count. The Gaussian is right to within about ±1.5 bits in both
  directions, consistent with the 1–2 bits seen for cyclic stage 2 at 523 and 573.
- The stage-2 total per y (H₂ = 95.5, C₂ = 67.5, so +28 bits) splits into C_dual = 25.6 and C_cross ≈ 42. The cross
  bits measured per dual solution (42.6–43.5) agree.

### 4.6 Cost [M]

- Setup (stage 1, dual lists): 20 s. Dual join: 2.3 s.
- First core: after 173 dual solutions, 5.5·10⁷ ≈ 2^25.7 join entries, about 3 minutes on one core.
- Steady state: about 0.25 s per dual solution once the fibre cache is warm, and about one core orbit per 600 dual
  solutions: 2–3 minutes per orbit.
- **≈ 2^26 operations per solution, far inside the 2^45 budget (B).**

### 4.7 Verification and bundle [M]

| check | tool | result |
|---|---|---|
| Σ_k PAF_G = 3732·δ over Z₃ × Z₃₁₁ | search, exact integer | yes |
| cyclic PAF on Z₉₃₃ | search (FFT, rounded) | yes |
| H Hᵀ = 3732·I | search (float64 BLAS) | yes |
| PAF by integer dot products; GS array; H Hᵀ = 3732·I in int64 blocks | `scripts/verify_h3732.py` (no shared code) | yes; sha256 of packed rows `8fea9ea2…a437e` |
| same, plus the published normalized matrix equals the normalization of the built H | `pub/H3732/verify.py` (numpy only) | yes |
| the sequences really are constant on the stated slice orbits | inline check in the commit session | yes |

Bundle: `pub/H3732/sequences_Z933.txt` (header gives the slice and the generators), `pub/H3732/H3732_normalized.txt.gz`,
`pub/H3732/verify.py`. b = 0 has no border, so there is no `border.json`. Packed matrix: `data/H3732_packed.npz`.
Row sums (A, B, C, D) = (−61, −1, −1, 3), and 61² + 1 + 1 + 9 = 3732.

---

## 5. Positive control for b ≥ 3: a new H(1388) in 347's slice [M]

(`theory/ab_prod347.py`, `theory/ab_stage2.py`, `theory/ab_assemble.py`; `logs/ab/prod347/`)

- **Slice.** G = Z₂ × Z₃ × Z₃ × Z₁₉; S = ⟨−1 on both Z₃⟩ × ⟨×7 on Z₁₉⟩ (index 6, f = 3). Stage 1 is fixed to the
  published y, so the Q-compression, and hence the template border, is fixed.
- **Stage 2.** The hierarchical solver of `theory/ab_stage2.py`:
  - S/D split on the Z₂ axis;
  - Wagner middle level along the H-part Z₃, with energy-pinned row values chosen by weighted sampling;
  - S-join, then the vectorized exact D-check.
- **Run.**
  - setup 256 s;
  - 11 middle values; 2^14.2 S-solutions in total (431–5,686 per middle value);
  - 1,420 s of joins and D-checks;
  - one D-solution, verified (Σ PAF = T exactly).
  
  It was assembled with the template border and the recomputed forced corner into a 1388 × 1388 matrix, which is
  Hadamard (search, and independently `scripts/verify_packed_h.py`; sha256 `27dd9ae6…0e07`). The core differs from
  the published one in 540 of 1,368 entries, with the same compression.
- **Calibration.** The level Gaussians on 347's own y gave D_def = 16.2 bits, i.e. about 2^16 S-solutions per core
  [E]. One core arrived at 2^14.2. That is one event, consistent to within about 2 bits.
- **Cost.** About 28 minutes on one core, roughly 2^35 elementary operations. This meets `report/bge3_answer.md` §7's
  request to reproduce 347's slice in ≲ 2^40.

---

## 6. b = 0 at the other composite orders

(`logs/ab/sweep_b0.jsonl` (screens), `logs/ab/b0_rank2.log` (stage 1 + stage-2 Gaussians, `theory/ab_b0_stage1.py`
with the DP fibre sampler `ab_gauss.stage2_sampled`, checked against the enumerating estimator on 933: 67.42 vs
67.47 bits))

| v | slice (staging prime q, index s, f, multiplier U on Z_F) | screens | stage 1 [M] | stage 2 per y [E] | status |
|---|---|---|---|---|---|
| 933 | q = 311, s = 10, f = 31 | pass | 36 | +28.0 (C₂ 67.5) | **built** (§4) |
| 933 | q = 311, s = 31, f = 10 | pass | 14,184 | +9.8 (C₂ 327.7) | not needed |
| 669 | q = 223, s = 37, f = 6 | pass | 131,006 | +30.8 (H₂ 407, C₂ 376.5) | fibres 3 × 37; the period field has prime degree 37, so there is no intermediate level |
| 917 | q = 131, s = 26, f = 5, U = ⟨2⟩ (order 3) on Z₇ | pass | 38,380 | +52 (H₂ 278, C₂ 225) | fibres 3 weighted rows × 26; dual vectors in {±1, ±5, ±7}^26 are too many to list |
| 955 | q = 191, s = 38, f = 5, U = −1 on Z₅ | pass | 274 | +113 (H₂ 419, C₂ 305) | same obstacle (3 rows × 38) |
| 917 | q = 131, s = 13, f = 10, U = ⟨2⟩ | pass | 96 | −27.8 (mahal 4.2) | dead [E] |
| 917 | q = 131, s = 26, U = Z₇^* (order 6) | pass | 5 | −411 | dead [E] |
| 917 | q = 131, s = 13, U = −1 | pass | **0** | – | **empty [M]** |
| 955 | q = 191, s = 19, U = −1 | pass | **0** | – | **empty [M]** |
| 917 | q = 131, s = 10 or 13, U = 1 | pass | not enumerated (7 orbits × 22 letters) | – | fibres 7 × 10 or 7 × 13 (≈ 2^50): the cross join needs another level |
| 955 | q = 191, s = 19, U = 1 | pass | running | – | fibres 5 × 19 |
| 917, 955 | q = 7 or 5 (index 2 or 3) with an order-13 or order-19 multiplier on Z₁₃₁ or Z₁₉₁ | pass | not run | – | stage 1 is itself a b = 0 problem on Z₁₃₁ or Z₁₉₁ |
| 537, 669, 789 | every s ≤ 6 | **A fails** | – | – | **empty [M]** |
| 917, 955 | s = 2, 5; 955 also s = 10 (all U); 917 s = 10 with U ≠ 1 | **A fails** | – | – | **empty [M]** |
| 537, 789 | s = 89, 131 (f = 2) | pass | – | – | 270- and 396-orbit slices, out of reach |

**Why 933 and nothing else.**
- The 933 pipeline needs three things:
  - a small G' (3 rows);
  - an index s ≈ 10, so that the dual vectors (4^s) can be listed and ranked;
  - fibres (≤ 3^s) that can be enumerated.
- Among the open composite v, only 933 has a prime factor q with a divisor s ≈ 10 of q − 1 that survives the
  screens. 669, 537 and 789 have q − 1 = 2·3·37, 2·89 and 2·131: the mid-size indices are missing or dead.
- 917 and 955 have F = 7 and 5. With U = 1 their fibres are too large; with the U that make 3 rows, either stage 1
  is empty (U = −1) or stage 2 is negative (index 13), or s is 26 or 38 and the dual level has too many vectors.
- The next level these need is a split of the dual level. For 917 at s = 26 the period field of index 26 has a
  degree-13 subfield (26 = 2·13), which gives a natural two-step dual. This is the first computation proposed in §9.

---

## 7. b ≥ 3 at the open orders

(`theory/ab_sweep.py` screens, `logs/ab/sweep_bge3.jsonl` (384 configurations); exact stage 1 `theory/ab_s1sd.py`;
bordered batches `theory/ab_cbatch.py` → `logs/ab/c4_batch.jsonl`, `logs/ab/c6_bordered.jsonl`; ranking
`theory/ab_rank.py` → `logs/ab/rank_out.jsonl`)

**Screens [M].** Screens A and B have 0 false negatives on 100 brute-forced small slices.

| v, c | configurations | A fails (empty) | pass A and B | other |
|---|---|---|---|---|
| 659, 4 | 54 | 27 | 27 | – |
| 659, 6 | 20 | 14 | 6 | – |
| 719, 4 | 6 | 3 | 3 | – |
| 719, 6 | 31 | 5 | 3 | 23 (B undecided or memory) |
| 751, 4 | 79 | 16 | 18 | 45 |
| 823, 4 | 57 | 30 | 13 | 14 |
| 839, 4 | 68 | 16 | 16 | 36 |
| 839, 6 | 12 | 8 | 4 | – |
| 863, 4 | 57 | 25 | 16 | 16 |

**Exact stage 1 [M].**
- 863, c = 4, G = Z₂² × Z₅ × Z₄₃, q = 43, s = 6, U = −1 on Z₅: all 31,188 α-quadruples, **0 solutions** (6,941 s).
- 863, c = 6, the 347-shaped configuration restricted to 347's border class: 272 admissible compressions, 2,072
  α-quadruples, 224 S-solutions, **0 stage-1 solutions**.
- 599, c = 6, q = 11, s = 2, H' = Z₉:
  - U = ⟨2⟩ (order 6) on Z₉: 80,910 α-quadruples, **0**;
  - U = −1 on Q's Z₃ (with or without −1 on H'), combined with ⟨2⟩ on Z₉: 29 α-quadruples each, **0**;
  - U = −1 on Q's Z₃ and on H', combined with ⟨4⟩: 1,446 α-quadruples, **0**.
- 571, c = 4, q = 71, H' = Z₂, s = 2, 5, 7, 10, 14: screen A fails (no admissible row quadruple), **empty**. The
  ranking's positive stage-1 Gaussians for s = 10 and 14 (+2.3 and +18.0) were wrong: one more instance of
  "whole-slice estimates decide nothing".

**Survivors and stage-2 size [E]** (`logs/ab/rank_out.jsonl`).

| v, c | slice | stage 1 [E] | stage 2 per y [E] | compare |
|---|---|---|---|---|
| 659, 4 | q = 41, s = 8, f = 5, H' = Z₄ | +54.8 | H₂ 385, C₂ 285, +100.7 | 347 control: C₂ 118; 523: C₂ 146 |
| 659, 4 | q = 41, s = 8, H' = Z₂² | +27.3 | H₂ 385, C₂ 274, +111 | |
| 659, 4 | q = 41, s = 10, f = 4, H' = Z₄ | passes A | (estimator ran out of memory) | |
| 719 / 863 / 917, 6 | best ones | – | whole-slice −13 to −270, except 863 s = 4 ⟨3⟩ (+37) and 917 s = 6 −1 (+133) | not yet exact |

- The best b ≥ 3 survivor, 659 at c = 4, is a stage 2 roughly 140 bits larger than anything solved (C₂ 285 against
  118 and 146). The 347 and 523 cost model scales with the D-deficit. It is not computed here, but at this size it is
  certainly above 2^45.
- So no b ≥ 3 configuration of an open order meets (B) yet. The exact c = 4 and c = 6 batches (all configurations
  with ≤ 26 G'-orbits) were still running when this was written; their records are in the logs named above.

---

## 8. Ranked table (D)

Status: **built** / exact obstruction **[M]** / estimate [E] / not reached. Higher rows are better.

| rank | v | builder | configuration | status | next computation |
|---|---|---|---|---|---|
| 1 | 933 | g (b = 0) | Z₃ × Z₃₁₁, index 10 | **built: H(3732)** | – |
| 2 | 917 | g | Z₇ × Z₁₃₁, s = 26, U = ⟨2⟩ | 38,380 y's, +52 per y [E] | a two-step dual through the degree-13 subfield (§9.2) |
| 3 | 955 | g | Z₅ × Z₁₉₁, s = 38, U = −1 | 274 y's, +113 per y [E] | as above, with the degree-19 subfield |
| 4 | 669 | g | Z₃ × Z₂₂₃, s = 37 | 131,006 y's, +30.8 per y [E] | no subfield; a norm-equation approach, or c = 10 b ≥ 3 |
| 5 | 659 | u/v/w, c = 4 | q = 41, s = 8 | passes A, B; stage 1 +55 [E]; C₂ ≈ 285 | exact stage 1 (c = 4 batch), then the level decomposition |
| 6 | 669, 789 | c = 10 | – | frames exist (Lemma 1) | corner-flatness search (§9.3) |
| 7 | 917 | g | s = 10 or 13, U = 1 | passes A, B | a smarter stage 1 (22^7 letters); fibres need a split |
| 8 | 863 | c = 6 | q = 13, s = 4, U = ⟨3⟩ on Z₁₁ | whole +37 [E] | exact bordered stage 1 |
| 9 | 917 | c = 6 | q = 19, s = 6, −1 | whole +133 [E] | exact bordered stage 1 (L₁ = 26) |
| – | 933 | g | s = 31 | 14,184 y's, +9.8 [E] | not needed |
| – | 917 | g | s = 13, U = ⟨2⟩; s = 26, U = Z₇^* | −27.8, −411 [E] | dead |
| – | 917, 955 | g | s = 13, U = −1; s = 19, U = −1 | **stage 1 empty [M]** | – |
| – | 537, 669, 789, 917, 933, 955 | g | the A-failing indices listed in §6 | **empty [M]** | – |
| – | 571 | c = 4 | every GP slice (q = 71) | **empty [M]** (screen A) | 571's other c (22, 26) |
| – | 863 | c = 4 | q = 43, s = 6, U = −1 on Z₅ | **stage 1 empty [M]** | – |
| – | 863 | c = 6 | 347-shape, 347 border class | **stage 1 empty [M]** | – |
| – | 599 | c = 6 | q = 11, s = 2, the five U's above | **stage 1 empty [M]** | the remaining 599 c = 6 configurations |
| – | 659, 719, 751, 823, 839, 863 | c = 4, 6 | 144 configurations failing screen A (§7) | **empty [M]** | – |

---

## 9. Next computations

1. **More H(3732)s and a census.** `theory/ab_b0z3.py --q 311 --y I` for the other 34 stage-1 solutions; each
   gives cores in minutes. It is not needed for (A).
2. **917, s = 26, U = ⟨2⟩: a two-step dual.**
   - Fix the class sums over the 13 cosets of the order-10 subgroup ⟨M, −1⟩ ⊃ M. These are pair sums
     T_j + T_{j'} over C ∪ −C.
   - Those sums satisfy the equations of the characters that are trivial on ⟨M, −1⟩ (the real Gauss periods): an
     exact sub-level with 13 columns.
   - Then the split T_j − T_{j'} meets the imaginary part.
   - First computation: the exact list of 13-column sub-level vectors per block under the spectral bound, and the
     Gaussian of the two dual sub-levels, calibrated on the same shape at v = 7·q for small q ≡ 1 (mod 26).
3. **c = 10 corner flatness (669, 789).**
   - Enumerate the Q-level solutions z on Z₁₀ and Z₂ × Z₅ for h = 66 and 78 (screen A's DP lists them).
   - For each z, search frames K, Kc from H₂ ⊗ H(20) (Lemma 1), with block-internal signed column permutations, for a
     flat corner. Lemma 2 reduces this to ±1-ness of K = −C·Kc·Cbᵀ/(4b).
   - A hit gives a c = 10 template, and §3's transfer then gives both orders.
4. **659, c = 4.** Finish exact stage 1 (`logs/ab/c4_batch.jsonl`). Then compute the level decomposition (H_S, C_S,
   h_D, C(D|S)) with `theory/ab_levels.py` on its y's, and cost the hierarchical solver.

---

## 10. Reproduce

```
python3 theory/ab_verify_published.py                 # section 1.4
python3 theory/ab_exact_small.py --sweep 131 --rmax 22 --out FILE      # screen validation
python3 theory/ab_sweep.py --v 659 --c 4,6 --out FILE  # screens A, B per configuration
python3 theory/ab_sweep.py --b0 --v 917 --out FILE
python3 theory/ab_s1sd.py @logs/ab/cfg_347_s6.json --max-sol 300 --out FILE   # exact stage 1
python3 theory/ab_levels.py @logs/ab/cfg_347_s6.json core:data/core_v347_b5.txt --amap '[[0,2],[1,3],[2,3],[2,19]]'
python3 theory/ab_assemble.py --check-template        # H(1388) rebuilt from the 347 template
python3 theory/ab_template.py                         # templates of 347, 479, 359 (rebuilt, corner checked)
python3 theory/ab_border_transfer2.py                 # border transfer to every open order with c = 4, 6, 8
python3 theory/ab_prod347.py --tag a1 --R 800 800 --K 4 128 --stop-on-core   # positive control (347 slice): new H(1388)
python3 scripts/verify_packed_h.py logs/ab/prod347/H1388_a1_438c58a201c4084b_1.npz 1388
python3 theory/ab_b0z3.py --q 71 --y 0 --K 200 --max-dual 300 --min-score 28 --tag ctl71       # same-shape control
python3 theory/ab_b0z3.py --q 311 --y 0 --K 200 --min-score 30 --tag p311e_y0 --stop-on-core   # H(3732): minutes
python3 scripts/verify_h3732.py logs/ab/b0z3/seq_p311e_y0_1_Z933.txt data/H3732_packed.npz      # independent check
cd pub/H3732 && python3 verify.py                     # standalone bundle check
python3 theory/ab_b0z3_calib.py --q 71 --y 0 --K 50 --n 6 --exact                              # calibration
python3 theory/ab_b0_stage1.py '<cfg json>' --max 3 --N 8000                                   # b = 0 stage 1 + stage 2 [E]
```
