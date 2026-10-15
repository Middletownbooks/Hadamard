# Alpöge-style builders (b = 0 and b ≥ 3) with the staged method

Date: 2026-10-14. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_alpoge_design.md`.
Scripts: `theory/ab_*.py`. Logs: `logs/ab/`. Bundle: `pub/H3732/`.

**Tags.** [M] exact computation (script named); [P] proved here (proof in the text); [E] Gaussian estimate (estimator
and calibration named); [L] literature recollection, unchecked.

---

## 0. Bottom line

1. **(A) H(3732) is built: v = 933, open on the prompt's list; b = 0 (Goethals–Seidel, array `g`) [M].**
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
   - The cross-level Gaussian (lattice-exact, `theory/ab_b0z3_calib.py`) is within about 1.5 bits of the exact core
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
   - **c = 10 (669, 789) [P, M] (§3.2).**
     - 359's c = 8 border is "Singer-type": monomial compressions and an order-7 automorphism.
     - With one constant compression block the c = 10 core is fixed by all of Sym(Z₁₀), so the same ansatz exists with a
       group of order 9.
     - An exact search (Z₉ and Z₃ × Z₃, all 5,198 first-blocks, 811,680 frame × Q corner joins) finds **no c = 10
       template**.
     - The same search finds c = 4, 6, 8 templates in seconds. That gives new border classes for c = 4 and c = 6, fitting
       every applicable open order.
     - c = 10 borders with smaller or no symmetry remain open.
5. **(C) Exact obstructions [M]** (§7). Among them:
   - b = 0:
     - these GP slices are empty: 537, 669 and 789 at every s ≤ 6; 933 at s = 2, 5; 917 (q = 131) at s = 2, 5, and at
       s = 10 with any multiplier on Z₇; 955 (q = 191) at s = 2, 5, 10;
     - 955 at s = 19 with −1 on Z₅, and 917 at s = 13 with −1 on Z₇, have no stage-1 solution.
   - b ≥ 3:
     - every GP slice of 571 at c = 4 (q = 71, s = 2, 5, 7, 10, 14) dies at the Q-level;
     - 863, c = 4, s = 6, U = −1: exact stage 1 is empty (31,188 α-quadruples);
     - 863, c = 6, q = 13, s = 4, H' = Z₁₁, U = ⟨−1 on Z₃ × Z₁₁, ×3 on Z₁₁⟩, in 347's border class: empty;
     - 599, c = 6, q = 11, s = 2, H' = Z₉, with five multiplier groups (§7): empty;
     - **the small Paley-type family** (index 2, 3 or 4 on a prime q ∥ h plus one multiplier group on the
       cofactor, ≤ 16 G'-orbits), decided exactly (§7):
       - at c = 4, 53 of the 61 configurations are decided empty:
         - 48 at the Q-level;
         - 4 with empty stage 1: 751 q = 11, 839 q = 19, 955 q = 17, 991 q = 19;
         - 839 at q = 11 with ⟨4⟩, empty after exact stage 2 on all 40,960 y's.

         Two more are dead in stage 2 [E] (751 q = 17, 991 q = 13, about −72 and −78 bits per y), and 6 are
         undecided (over the caps);
       - at c = 6 all 19 are empty (347's border class).
6. **(D)** The ranked table is §8.
   - The other composite orders reach b = 0 GP slices whose stage 2 is Gaussian-positive, but whose fibres are too
     large for the 933 pipeline. They need one more level (§6):
     - 917 at s = 26 with an order-3 multiplier: +52 bits per y, over 38,380 y's;
     - 955 at s = 19: +35.9 bits per y, over 3.27·10⁶ y's;
     - 955 at s = 38 with −1: +113 bits per y;
     - 669 at s = 37: +30.8 bits per y, over 131,006 y's.
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
H₄ ⊗ 1_cᵀ, and the remaining 4c − 4 rows form a frame. H(2c) is known for every multiple of 4 below 668 [L]; 668
is the smallest open order. So frames exist for every even c ≤ 332. Checked numerically for c = 6 and c = 10 (Paley I H(12) and H(20)) [M]. ∎

**Lemma 2 (the corner is always orthogonal; only its flatness is in question).** This restates
`report/bge3_answer.md` §4.1 ("C is automatically √(4b) times an orthogonal matrix") with a proof, and adds the
explicit inverse K = −C·Kc·Cb/(4b).
- Let V ⊂ R^{4c} be the block-balanced subspace (dimension 4b), Π_V the projection onto it, and K, Kc frames, so
  that Kᵀ K = Kcᵀ Kc = 4c·Π_V.
- The compressed GS array Cb of the Q-compressions z_k satisfies Cb Cbᵀ = I₄ ⊗ M(Σ_k PAF_Q(z_k)), M(·) the group
  matrix (the GS identity over Q). Its eigenvalues are Σ_k |ẑ_k(ψ)|² = 4b on V (ψ ≠ 1, the Q-level equation of §1.2)
  and 4v on V^⊥.
- Cb and Cbᵀ map block-constant vectors (V^⊥) to block-constant vectors. Each is the adjoint of the other, so each
  also preserves V = (V^⊥)^⊥.
- Likewise Cbᵀ Cb = I₄ ⊗ M(Σ_k PAF_Q(z_k)), since the transposed array is again of GS type, so Cbᵀ Cb = 4b on V.
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

**Validation of the search formulation [M]** (`theory/ab_border10.py --c 6 --published`). 347's compression,
shifted block by block to ±1 values, is an h = 1 Q-level solution (Σ_k PAF = 24, 4, 4, 4, 4, 4). With the published
Kc, the admissible set {k block-balanced ±1 : −k·Cbᵀ·Kcᵀ/24 ∈ {±1}^20} has exactly 40 elements: ±(the 20 rows of
the published K). So given (z, Kc), the other frame is rigid. A template is a pair of frames that are mutually
unbiased after the rotation Cb|_V/√(4b). Consistent with `report/bge3_answer.md` §4.1, a random frame has no
partner, and borders must be built together with the compression.

**Consequences.**
- c ≥ 10 is not blocked by frames. The open question is the flatness of the corner for some Q-level compression z
  and some pair of frames.
- The published templates are flat solutions for c = 4, 6, 8, and §3's transfer theorem carries them to every open
  order with the same c.
- For c = 10 (the only b ≥ 3 option of 669 and 789) the first computation is a search for a flat pair: the Q-level
  solutions z on Z₁₀ or Z₂ × Z₅, against frames from H₂ ⊗ H(20) and their signed column permutations within blocks
  (§9). §3.2 does this for the compressions that carry the most symmetry.

### 3.2 c = 10: the Singer-type ansatz of 359's border [P, M]

(`theory/ab_singer.py`, `theory/ab_singer_layout.py`, `logs/ab/singer/`; the automorphism computations use
`theory/ab_hadaut.py`'s graph with weighted edges.)

**What the c = 8 border of 359 is [M].**
- Untwisted (K's columns multiplied by the character ψ₆ of Q = Z₂³, Kc's by ψ₁, in the stored coordinates) its frames
  are block-balanced and its compressions are *monomial*: z_k = s_k + a_k·δ_{q*} with a = (2, −2, 2, 4), Σ a_k² = 28 = 4b.
- The template [[C, K], [Kcᵀ, Cb]] has an automorphism of order 7: on every core block it fixes one point and
  cycles the other seven (a Singer cycle of GF(8), linear on Z₂³), and the 28 border rows and columns form four
  regular orbits. Up to ±1 the template's automorphism group has order 14.
- In these "Singer coordinates" the corner is a 4 × 4 array of 7 × 7 circulants, and K, Kc are [e ⊗ 1₇ | a 4 × 4 array of
  circulants] with e, e′ ∈ H(4). Only three sequence classes occur: Legendre(7) `+++-+--`, the interval `++++---`, and
  1 − 2δ.
- The c = 4 and c = 6 templates (479, 347) have automorphism groups of order 2 (mod ±1); 347's is an H(44) that is
  not Paley-equivalent.

**Lemma 3 (monomial cores are Sym(Q)-symmetric) [P].** Let z_k = a_k·δ₀ (constants do not matter on balanced frames)
and suppose one a_k = 0, or Q is an elementary abelian 2-group. Then for every permutation σ of Q the core array Cb (347
layout) is fixed by σ acting on two blocks and ρσρ on the other two (ρ = negation on Q), on rows and columns alike.
- *Proof.* A diagonal block a₀·δ₀(h − g) is a₀·I, fixed by (π, π). An off-diagonal block a_k·δ₀(±(g + h)) is a_k·ρ,
  fixed by (π, π′) exactly when π′ = ρπρ.
- With the zero block in position D, the nonzero off-diagonal blocks join {0, 3} to {1, 2}. So π = σ on blocks 0, 3 and
  ρσρ on 1, 2 satisfies every block.
- With the zero in A, rows take σ and columns ρσρ.
- With Q = Z₂^m, ρ = 1.
- Checked: |Aut| of the c = 10 core is 2·10!·8 (igraph), and 20 random σ fix it [M]. ∎

So any group Γ of order b = c − 1 acting regularly on Q∖{0} fixes Cb. V_Q is then the regular representation of Γ,
frames can be unions of four regular orbits, and in Singer coordinates (block J = {∞} ∪ Γ) Cb = Q_S ⊗ I_c. Here Q_S
is a 4 × 4 integer matrix with Q_S Q_Sᵀ = 4b·I. For c = 10 the monomial profile must be a = (4, 4, 2, 0) up to order
and signs: the only other profile, (6, 0, 0, 0), makes Cb|_V six times a signed permutation, and the corner entries
±40/6.

**Lemma 4 (Γ-invariant templates) [P].** With Γ acting regularly on four border orbits, write
K = [e ⊗ 1 | circ(f)], Kc = [e′ ⊗ 1 | circ(g)]: f, g are 4 × 4 arrays of ±1 functions on Γ, and e_IJ = −Σf_IJ.
- K is a frame iff every row-block (f_I0, …, f_I3) has Σ_J PAF(f_IJ) = −4 off 0, and every two row-blocks have
  vanishing cross-correlation Σ_J [e_IJ e_I′J + PCF(f_IJ, f_I′J)(d)] for all d. The second condition is
  Fourier-orthogonality at every character, so it is invariant under shifting either row-block.
- The corner block C_{II′} = −(1/4c)·K_I·(Q_Sᵀ ⊗ I)·Kc_{I′}ᵀ depends only on row-block I of K and row-block I′ of Kc.
- Hence a Γ-invariant template is a K_{4,4} biclique in the bipartite "flat corner" relation between row-blocks,
  whose two sides are 4-cliques of the orthogonality relation. ∎

**Search [M]** (`theory/ab_singer.py`):
- (A) The first Kc row-block, up to affine maps of Γ and signed column permutations: 4,518 classes for Z₉, 680 for
  Z₃ × Z₃.
- (B) Its orthogonal partners by a hashed join on (PAF sums, cross-correlation), then 3-cliques give all Kc frames
  through it.
- (C) For each frame and each of the 24 right classes of Q_S (orthogonal, rows signed permutations of (4, 4, 2, 0); all
  one two-sided class, which the layout realises), every K row-block with a ±1 corner against the frame, found by a
  hashed join on the corner entries mod 8c and the PAF sums.
- Then 4-cliques of mutually orthogonal row-blocks give K. Each frame is processed once, from its minimal row-block.
- A hit is assembled and checked exactly. `theory/ab_singer_layout.py` maps it to the 347 layout and re-checks it with
  `ab_assemble.forced_corner`.

**Control [M]: c = 8, Γ = Z₇, a = (2, 2, 2, 4).**
- The search finds templates from the first first-block on (2 in the first 40 of 185 first-blocks, 77 s; 8 in about
  120 with the first version of the code).
- Each maps to the standard 347 layout over Z₂³ with z_k = a_k·δ₀ and passes `forced_corner` exactly
  (`logs/ab/singer/c8_control/*_layout.npz`).
- So the method rediscovers c = 8 borders from nothing, in the plain layout rather than 359's twisted one.

**c = 10 [M]: no Singer-type border.**
Profile a = (4, 4, 2, 0); Γ = Z₉ and Γ = Z₃ × Z₃, the two groups of order 9. Both acting regularly on Q∖{0} means
fixing one point and transitive on the other nine (`logs/ab/singer/c10_*.jsonl`).

| Γ | first-blocks | Kc frames processed | (frame, Q) corner joins | templates |
|---|---|---|---|---|
| Z₉ | 4,518 (all) | 28,942 | 694,608 | **0** |
| Z₃ × Z₃ | 680 (all) | 4,878 | 117,072 | **0** |

- Frames exist in abundance, and what fails is always the corner:
  - Z₉: 78,084 frame instances through the first-blocks, with 1,134 to 4,626 orthogonal partners per first-block
    (249 first-blocks have no frame);
  - Z₃ × Z₃: 14,915 frame instances, with 1,944 to 3,690 partners per first-block.
- Spike positions do not matter. With z_k = s_k + a_k·δ_{q_k}, every cycle of the nonzero-block graph composes to the
  identity: e.g. P_{q₀} R_{q₁} P_{q₀} = R_{q₁} and R_{−q₁} R_{−q₂} R_{q₁} = R_{q₂}, with P translations and R reflections
  of Q. So Cb is again Q_S ⊗ I in Singer coordinates. And Q_S ranges over the whole (single) two-sided class, so every
  placement of the zero block is covered.
- Order 9 forces monomial compressions. Let a block-preserving automorphism group of order 9 act on each core block
  as a 9-cycle or as regular Z₃ × Z₃, fixing one point p. Putting g = p in z(π′h − πg) = z(h − g), or in its g + h
  analogue, shows that every z_k is constant on the other nine points of Q, so z_k is monomial.
- **So no c = 10 border (any compressions) has a block-preserving automorphism group of order 9 acting regularly on
  Q∖{0} in every core block [M].** This includes every template invariant under a larger group containing such a
  subgroup, e.g. AGL(1, 9), PSL(2, 9) ≅ A₆ and PΓL(2, 9) on the ten points of PG(1, 9).
- Trend: Singer borders are plentiful at c = 4 and 6, rare at c = 8 (two in the first 40 of 185 first-blocks), and
  absent at c = 10.
- c = 12 probe (b = 11 prime, profile (6, 2, 2, 0), Γ = Z₁₁): no template in the first 250 of 75,573 first-blocks
  (517 frames × 24 Q). c = 4, 6, 8 all produced templates within their first few first-blocks. The probe stopped
  when the container was reclaimed and was not resumed (about 30 core-hours for the rest; restart-safe:
  `logs/ab/singer/c12_cyc_a6220_p*of4.jsonl`).
- The cyclic orders 3, 5, 7 are primes, and 7 = 2³ − 1 carries a Singer difference set. Order 9 is neither prime nor of the
  form 2^m − 1, and the c = 8 templates are built from Legendre(7) and an interval (above).

What is left open for c = 10:
- smaller symmetry groups, e.g. Z₃ acting as three 3-cycles: 12 orbits of size 3, so frames are 12-cliques and the
  method above does not scale;
- non-monomial compressions (h = 1: the SDS profiles (5,4,2,2) and (4,4,4,1), i.e. H(76) completions; even h: integer
  u with Σ_k |û_k(ψ)|² = 9), which have no symmetry beyond the multipliers and translation by 5 of Z₁₀.

**c = 4 and c = 6 [M]: Singer borders are plentiful.**
- c = 6, Γ = Z₅, a = (4, 2, 0, 0): 68 templates from 10 of the 15 first-blocks (30 s).
- c = 4, Γ = Z₃, a = (2, 2, 2, 0): 40 templates from 3 of the 4 first-blocks.
- Samples map to the 347 layout over Z₆ (a = (0, 0, 2, 4)) and Z₄ (a = (0, 2, 2, 2)) and pass `forced_corner`
  (`logs/ab/singer/*_layout.npz`).
- The c = 6 monomial profile is bge3's compression type (0, 0, 5, 20), for which `report/bge3_answer.md` §4.1 found no
  border among sampled frames. It has borders.
- Both profiles fit the row equation Σ_k (a_k + c·s_k)² = 4v (s_k ≡ h mod 2) at every applicable open order:
  - c = 4: 572–1,152 solutions per order;
  - c = 6: 528–1,056 solutions at 599, 659, 719, 839, 863, 917, 947, 971.
- So these are new border classes for c = 4, 6, beyond the published ones. A slice whose stage-1 Q-compressions are
  monomial can now be bordered.
- In Singer coordinates the c = 4 and c = 6 templates use a single "interval" class (`++-`, `+++--`) for every
  sequence of K and Kc. The corner uses that class and one near-constant class.

**Methods that failed the control [M].**
- *CP-SAT* (`theory/ab_border_sat.py`: entries of C, K, Kc as Booleans, products as XNORs, every frame, linking and
  corner equation linear). It is correct: with the published solution's entries fixed it reports feasible in 0.5 s.
  But it cannot rediscover the c = 8 border even with that border's full automorphism group imposed (184 free
  variables, 300 s, UNKNOWN). It needs about half the variables given: with 5% of the entries fixed it solves in
  0.6 s, with 2% it times out at 120 s.
- *Douglas–Rachford / RRR* projections between "±1 with the core fixed" and "√n·orthogonal"
  (`theory/ab_border_dr.py`) do not complete even the c = 6 (H(44)) case.

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
| the sequences really are constant on the stated slice orbits | `scripts/verify_h3732.py --slice 311 10` | yes |

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
| 955 | q = 191, s = 19, U = 1 | pass | 3,265,000 | +35.9 (H₂ 329, C₂ 293) | fibres 5 × 19 (≈ 2^53 per block): blocked at the cross level |
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
- **Cost rule of the dual/cross pipeline [E, checked on 933].**
  - A dual solution with fibre product 2^F yields about 2^{F − C_cross} cores. Its cross join costs at least
    2^{F/2} when the fibres are balanced.
  - So the cost per core is about 2^{C_cross − F/2} on the best dual solutions, and it needs F ≲ C_cross plus about
    15 for the dual join to supply enough dual solutions.
  - 933: C_cross ≈ 42 and F ≈ 34.6 give 2^{24.7} (measured 2^{25.7}).
  - 917 at s = 10, U = 1: F ≈ 4·45 and C_cross ≈ 30·4.2 = 126. The yield is ample, but each fibre is ≈ 2^45, so the
    2 + 2 cross join is ≈ 2^90 per dual solution.
  - 955 at s = 19: the same, with fibres ≈ 2^53.
  - The twisted multiplier ⟨(−1, m₀)⟩ with m₀ of order 2·19 on Z₁₉₁ (955), or 2·13 on Z₁₃₁ (917), gives a 933-like
    shape: rows {0}, {±1}, {±2} (and {±3}), 10 columns, and only 6⁵ or 8⁵ dual vectors per block. But −1 already lies
    in the slice group, so the cross equations are not halved: 20 cross equations (≈ 85 bits) against fibres of
    about 2^12 per block. The yield is about 2^{−40} per dual solution, so it is not viable [E, counting].
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

**Exact batches [M]** (`theory/ab_cbatch.py`: screen A, then exact stage 1, restricted for c = 6 and c = 8 to the
border class of the 347 and 359 templates; configurations with ≤ 26 G'-orbits; 1,200 s and a 4–4.5 GB cap each).
Final state: the c = 8 batch ran to completion; the c = 4 and c = 6 batches were stopped at 19:00 UTC, after about 3 hours:

| batch | configurations | empty (A fails) | stage 1 empty in the border class | survivors | memory or time cap | > 26 orbits |
|---|---|---|---|---|---|---|
| c = 4 (571, 599, 643, 659, 719, 739) | 64 | 21 | – | 2, both f = 2 slices of about 300–360 orbits: 599 (q = 149, s = 74: 1,500 y's) and 719 (q = 179, s = 89: 1,800 y's) | 10 | 31 |
| c = 6 (599, 659, 719) | 125 | 29 | 10 | 5 (659, see below) | 37 | 44 |
| c = 8 (599, 719, 751, 823, 839, 863, 919, 991) | 81 | 10 | – | 0 | 29 | 42 |

**659, c = 6 (h = 109 prime, G' = Q = Z₂ × Z₃).** Exact stage 1 in 347's border class
(`logs/ab/small/s1_659c6_*.json`):

| index s on Z₁₀₉ (f) | U | stage-1 y's | in the border class | stage 2 per y [E] (`ab_gauss.stage2_sampled`) |
|---|---|---|---|---|
| 27 (4) | 1 | 216 | 16 | H₂ 582, C₂ 522, **+60** (mahal 0.1) |
| 27 (4) | −1 on Z₃ | 8 | 8 | H₂ 388, C₂ 566, −178 |
| 36 (3) | 1 | 108 | 8 | H₂ 793, C₂ 523, **+270** (mahal 0.1) |
| 36 (3) | −1 on Z₃ | 4 | 4 | H₂ 529, C₂ 398, +130 (mahal 8.4: far tail, unreliable) |
| 54 (2) | −1 on Z₃ | 58 | 58 | not estimated (a 220-orbit slice) |

- These are the first b ≥ 3 configurations of an open order with exactly verified, bordered stage-1 solutions and a
  positive stage 2.
- Their stage 2, though, is 3–4 times the 347 control's (C₂ 118). Compressing along Q gives a b = 0-type dual level
  on Z₁₀₉ (target 4v·δ, since Σ_q T(q, t) = −4b + 4(c − 1) = 0 off 0). With 27–36 columns, that dual level is
  as unlistable as 669's at s = 37. **Not tractable with the present pipeline.**

**Small Paley-type c = 4 slices: an exactly decidable family [M]** (`theory/ab_paley_c4.py`, `logs/ab/paley_c4.jsonl`;
the H(2292) and 359 pattern). Q = Z₂², GP index s = 2 (the Paley classes; also s = 3, 4) on a prime q ∥ h, and one
cyclic multiplier group on a cyclic factor of H' = Z_{h/q}, keeping ≤ 16 G'-orbits.
- Then every block fibre of stage 2 has at most 2^16 elements, and the configuration is decided completely by:
  - screen A;
  - exact stage 1, either by direct enumeration under the energy bound with a symmetric energy-split join
    (`theory/ab_s1_direct.py`) or by the S/D solver;
  - the one-shot exact stage-2 join (`theory/ab_small_exact.py`).
- Any core is then bordered by `theory/ab_border_c4.py`. Its self-test rebuilds H(1916) from 479's core, untwisted.
- Validation: on the v = 35, c = 6 control, direct stage 1 plus exact stage 2 give 31,104 canonical cores, the
  brute-force count (11,943,936 = 384 · 31,104).
- Typical-y stage-2 Gaussians [E]: 751 (q = 11, ⟨2⟩ of order 8 on Z₁₇) about −5 bits per y; 839 (q = 11, ⟨4⟩ of
  order 9 on Z₁₉) about −11; 839 (q = 19, ⟨3⟩ of order 5 on Z₁₁) about 0.
- Results:
  - **c = 4 sweep over all open orders** (`logs/ab/paley_c4.jsonl`): 61 configurations with ≤ 16 G'-orbits.
    - 48 fail screen A: **empty [M]**. That covers every configuration of 571, 599, 659, 719, 823, 859, 863, 907,
      919 and 947.
    - 751, q = 11, s = 2, ⟨2⟩ of order 8 on Z₁₇: **stage 1 empty [M]**, decided two independent ways:
      - the direct method with Q-level pruning: 3,590 block vectors, 0.4 s;
      - the S/D solver: 12,992 α-quadruples, 0 S-solutions, 2,146 s.
    - 839, q = 19, s = 2, ⟨3⟩ of order 5 on Z₁₁: **stage 1 empty [M]** (2.6M block vectors, complete join).
    - 751, q = 17, s = 4, ⟨3⟩ on Z₁₁: 311,296 stage-1 solutions [M] in 102 PAF-signature classes. Stage 2 is
      about −72 bits per y [E] (mahal 5.8), so about 2^−54 cores are expected in all: **dead**.
    - 839, q = 11, s = 2, ⟨4⟩ of order 9 on Z₁₉: 40,960 stage-1 solutions [M] (12 PAF-signature classes). Exact
      stage 2 on every one: **0 cores, the slice is empty [M]** (2,340 s). The Gaussian had predicted −21 bits per y,
      about 0.02 cores in all.
    - 839, q = 11, s = 2, ⟨8⟩ of order 6 on Z₁₉ (16 orbits): undecided. The block lists exceed 6·10⁷ vectors even
      with Q-level pruning.
    - 955, q = 17, s = 2, ⟨3⟩ of order 6 on Z₇ (H' = Z₂ × Z₇): **stage 1 empty [M]** (230,422 block vectors).
    - 991, q = 19, s = 2, ⟨5⟩ of order 4 on Z₁₃: **stage 1 empty [M]** (72,552 block vectors).
    - 991, q = 13, s = 3, ⟨4⟩ of order 9 on Z₁₉: 98,304 stage-1 solutions [M] in 18 signature classes. Stage 2 is
      about −78 bits per y [E] (mahal 5.1), about 2^−61 cores in all: **dead**.
    - Undecided (over the 6·10⁷ block-vector cap, or out of memory at 7 GB): 839 q = 11 with ⟨8⟩; 955 q = 7,
      s = 2 and 3; 955 q = 17, s = 4; 991 q = 13, s = 4 with ⟨8⟩ and with ⟨4⟩.
    - **Totals for the 61 configurations: 53 empty [M] (48 at screen A, 4 with empty stage 1, 1 with exact stage 2
      over all its y's), 2 dead [E], 6 undecided.**
  - **c = 4 with 17–24 G'-orbits** (`logs/ab/paley_c4_o24.jsonl`): 20 configurations (739, 751, 823, 839, 863, 955).
    10 fail screen A (**empty [M]**); 10 exceed the direct stage-1 cap (6·10⁷ block vectors) and are undecided.
  - **c = 6 sweep (347's border class)** (`logs/ab/paley_c6.jsonl`): every configuration with ≤ 16 orbits (659,
    719, 839, 863, 947, 971: 19 in all) **fails screen A, empty [M]**. 599 and 917 have none that small.
  - Validation of the c = 6 path on the controls v = 35 and v = 47: stage 1, then the 347 border fit (8 bordered
    y's each), then exact stage 2. No bordered y of the controls carries a core, so the c = 6 assembly step itself is
    exercised only on the 347 positive control (§5).
- **Speed-up that made this family cheap [P].** For b ≥ 3, every block satisfies |ẑ_k(ψ)|² ≤ 4b for ψ ≠ 1, so its
  coset sums satisfy |z(q) − z(q′)| ≤ (2/c)(c − 1)√(4b) (≈ 5.2 for c = 4). Enumerating the letters coset by coset
  with this prefix bound cuts the block lists by factors of 10³–10⁴. The final exact bound |ẑ(ψ)|² ≤ 4b is applied
  at the end.

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
| 5 | 659 | u/v/w, c = 4 | q = 41, s = 8 | passes A, B; stage 1 +55 [E]; C₂ ≈ 285 | exact stage 1 ran out of memory (4 GB); then the level decomposition |
| 5′ | 659 | u/v/w, c = 6 | h = 109, s = 27 or 36 | **bordered stage-1 y's exist [M]**; stage 2 +60 to +270 [E], C₂ 400–520 | a two-level dual on Z₁₀₉ (27 = 3³, 36 = 2²·3², so the period field has subfields) |
| 6 | 669, 789 | c = 10 | – | frames exist (Lemma 1); **no border with an order-9 symmetry [M]** (§3.2) | borders with smaller symmetry or non-monomial compressions (§9.3) |
| 7 | 917 | g | s = 10 or 13, U = 1 | passes A, B | a smarter stage 1 (22^7 letters); fibres need a split |
| 8 | 863 | c = 6 | q = 13, s = 4, U = ⟨3⟩ on Z₁₁ | whole +37 [E] | exact bordered stage 1 |
| 9 | 917 | c = 6 | q = 19, s = 6, −1 on Z₃ and Z₈ (347's shape with H' = Z₈) | screen A passes (29 row quadruples); whole +133 [E] | exact bordered stage 1 or screen B: 2^26 α-vectors need a streaming residue join (this run ran out of memory) |
| – | 933 | g | s = 31 | 14,184 y's, +9.8 [E] | not needed |
| – | 917 | g | s = 13, U = ⟨2⟩; s = 26, U = Z₇^* | −27.8, −411 [E] | dead |
| – | 917, 955 | g | s = 13, U = −1; s = 19, U = −1 | **stage 1 empty [M]** | – |
| – | 537, 669, 789, 917, 933, 955 | g | the A-failing indices listed in §6 | **empty [M]** | – |
| – | 571 | c = 4 | every GP slice (q = 71) | **empty [M]** (screen A) | 571's other c (22, 26) |
| – | 863 | c = 4 | q = 43, s = 6, U = −1 on Z₅ | **stage 1 empty [M]** | – |
| – | 863 | c = 6 | 347-shape, 347 border class | **stage 1 empty [M]** | – |
| – | 599 | c = 6 | q = 11, s = 2, the five U's above | **stage 1 empty [M]** | the remaining 599 c = 6 configurations |
| – | 659, 719, 751, 823, 839, 863 | c = 4, 6 | 144 configurations failing screen A (§7) | **empty [M]** | – |
| – | 751 | c = 4 | Paley q = 11, ⟨2⟩ (order 8) on Z₁₇ | **stage 1 empty [M]** (two independent exact methods) | – |
| – | 839 | c = 4 | Paley q = 19, ⟨3⟩ (order 5) on Z₁₁ | **stage 1 empty [M]** | – |
| – | 839 | c = 4 | Paley q = 11, ⟨4⟩ (order 9) on Z₁₉ | 40,960 y's, **stage 2 empty on all [M]** | – |
| – | 955, 991 | c = 4 | q = 17 with ⟨3⟩ on Z₇; Paley q = 19 with ⟨5⟩ on Z₁₃ | **stage 1 empty [M]** | – |
| – | 991 | c = 4 | q = 13, s = 3, ⟨4⟩ on Z₁₉ | 98,304 y's [M]; stage 2 −78 bits per y [E] | dead |
| – | 839, 955, 991 | c = 4 | 6 Paley-type configurations over the stage-1 caps | undecided | a coset-wise meet in the middle for stage 1 |
| – | 751 | c = 4 | q = 17, s = 4, ⟨3⟩ on Z₁₁ | 311,296 y's [M]; stage 2 −72 bits per y [E] | dead |
| – | all open v | c = 4, 6 | Paley-type, ≤ 16 orbits: 48 + 19 configurations | **empty [M]** (screen A) | – |

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
3. **c = 10 borders (669, 789).** Done (§3.2): no border with a block-preserving order-9 symmetry, whatever the
   compressions. What is left:
   - [P] Structural constraint: √(4b) = 6 is an integer, so Cb|_V = 6·O with O a rational orthogonal map. A template
     is a pair of ±1 frames that are *mutually unbiased* after O: |cos| = 1/6 between all 36 × 36 pairs (§3.2's MUB
     form, here with one basis rotated by O).
   - Smaller symmetry. Z₃ acting as three 3-cycles plus a fixed point is compatible with frames (trace 0): 12 regular
     orbits of size 3. Its frames are 12-cliques, so §3.2's row-block method needs a different way to enumerate them
     (e.g. per-character Fourier blocks: V = V₁ ⊕ V_ω ⊕ V_ω̄, each of dimension 12).
   - Non-monomial compressions. At h = 1 a template is an H(76) containing the 40 × 76 partial matrix [Kcᵀ | GS(z)]
     (z an SDS of profile (5,4,2,2) or (4,4,4,1)); at even h, u is integral with Σ_k |û_k(ψ)|² = 9. Douglas–Rachford
     and CP-SAT both failed the c = 6 and c = 8 controls (§3.2), so a new method is needed.
   - If a border is found, its row equations are solvable at both orders: for the monomial profile, 192 solutions at
     669 and 120 at 789 (e.g. row sums (0, 16, 22, 44) and (0, 2, 4, 56)). The core search is then a c = 10 slice;
     669 with p = 11 has n′ = 60 (`report/unattempted_2000_3000.md`), far beyond the annealing wall.
   - c = 12 (an option at 599, 659, 719, 839, 863, 947, 971; not 669, 789): the Singer probe found nothing in the
     first 250 of 75,573 first-blocks (§3.2). Resuming it (four parts, about 30 core-hours) needs a session that
     stays active, or a machine of its own.
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
python3 scripts/verify_h3732.py pub/H3732/sequences_Z933.txt data/H3732_packed.npz --slice 311 10   # independent check
cd pub/H3732 && python3 verify.py                     # standalone bundle check
python3 theory/ab_b0z3_calib.py --q 71 --y 0 --K 50 --n 6 --exact                              # calibration
python3 theory/ab_b0_stage1.py '<cfg json>' --max 3 --N 8000                                   # b = 0 stage 1 + stage 2 [E]
python3 theory/ab_s1_direct.py '<cfg json>' --out FILE        # exact stage 1, direct (energy + Q-level pruning)
python3 theory/ab_small_exact.py '<cfg json>' FILE --tag T    # exact stage 2 for small fibres (one 2+2 join per y)
python3 theory/ab_paley_c4.py [--c 6] --max-orbits 16          # the Paley-type family, decided exactly, all open v
python3 theory/ab_border_c4.py                                # c = 4 border builder self-test (H(1916) from 479's core)
python3 scripts/make_bundle_bge3.py CORE.npy MODS NQ LAYOUT.json K.npy KC.npy C.npy OUTDIR   # b >= 3 bundle writer
python3 theory/ab_border10.py --c 6 --published               # corner-flatness search formulation, validated on 347
python3 theory/ab_singer.py --c 8 --group cyc --a 2,2,2,4 --max-first 40   # Singer borders, c = 8 control (seconds)
python3 theory/ab_singer.py --c 6 --group cyc --a 4,2,0,0                  # c = 6 (68 templates), likewise --c 4 --a 2,2,2,0
python3 theory/ab_singer.py --c 10 --group cyc --a 4,4,2,0 --part i/3      # c = 10, Z9: all 4,518 first-blocks, no template
python3 theory/ab_singer.py --c 10 --group z3z3 --a 4,4,2,0                # c = 10, Z3 x Z3: all 680 first-blocks, no template
python3 theory/ab_singer_layout.py logs/ab/singer/template_c6_cyc_a4200_f0_1.npz --Q 6   # map to the 347 layout, check
```
