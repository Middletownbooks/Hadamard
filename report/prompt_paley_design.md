# Prompt: algebraic design inside the Paley-coset ansatz (self-contained, for a fresh Claude Code session)

*For a maximum-reasoning Claude Code session with the repository `middletownbooks/hadamard`, branch
`claude/focused-planck-j1645f`. Tags: [M] verified in exact arithmetic in the repo, [E] Gaussian estimate,
[L] literature recollection (check it; web or paper search may be blocked, so do not depend on it).*

## 0. Goal and what I want back

Construct a Hadamard matrix of order 4v for one of the open v ∈ {523, 537, 571, 573, 599, 643, 647, 659, 669, 719, 739}.
The route here is the **Paley-coset ansatz** (§2), which the published v = 359 core satisfies exactly. Its search
version is about 2^300 for every open order, which is out of our reach (§4). I want **structure that removes most of
that search**: an algebraic or recursive design for some or all of the components (w, α, u) defined in §2.

Deliverables, in order of value:
- **(A)** an explicit construction that produces a core at an open v, verified by the repo's exact checks;
- **(B)** a design that reduces the remaining search for some open-order option to at most about 2^45 operations,
  with a positive control: the same design must reproduce a core at a smaller order, ideally a new 359 core;
- **(C)** a proof that a natural design family cannot work (an exact obstruction), so we stop pursuing it.

Write your results to `report/paley_design_answer.md` and your scripts to `theory/paley_*.py`. Commit and push.
- Every claim must be either [M] (with the script that checks it), a cited theorem with its hypotheses checked, or
  explicitly marked as a conjecture or estimate.
- Do not kill or modify running background jobs: `scripts/run_paley_bge3_screen.sh` writes
  `logs/paley_bge3_screen.jsonl`.

## 1. Architectures (exact)

**Goethals–Seidel (GS) core.** Four ±1 functions x_k on an abelian group G of order m. Write PAF_k(g) = Σ_y x_k(y) x_k(y+g)
and T(g) for the target of Σ_k PAF_k(g):
- **b = 0:** m = v, and T(g) = 0 for g ≠ 0.
- **b = 1 (bordered):** m = v − 1, T(g) = −4 for g ≠ 0, and the row sums are (±2, 0, 0, 0).
- **b ≥ 3 (inflated border, the 347/359/479 architecture):**
  - c = b + 1 is even, v + 1 = c(h + 1), m = v − b = c·h, and H ≤ G has index c;
  - T = −4b on H ∖ {0} and +4 off H; row sums satisfy Σ r_k² = 4v.
  - The border is an H(4c) with the four rows H₄ ⊗ 1_c removed, inflated by h. The 4b × 4b corner is forced by the
    border and by where the core's compression onto G/H sits, and must be ±1.
  - Borders always exist for c = 4. For c = 6 and 8 the 347 and 359 borders exist but must match the compression.
  - Details are in `report/bge3_construction.md` and `scripts/bge3_border_check.py`.

## 2. The Paley-coset ansatz and its exact reduction [M]

**The ansatz.**
- Take G = G′ × Z_p, with p prime, p ∥ m, and p ∤ c, so that for b ≥ 3, H = H′ × Z_p with [G′ : H′] = c.
- Let χ be the Legendre symbol mod p, δ the indicator of 0 in Z_p, and ε = 1 − δ.
- Every block restricted to every coset g′ + Z_p has the form

      x_k(g′, ·) = α·δ + u·χ + w·ε,   with α = ±1 and exactly one of (u, w) nonzero and ±1,

  i.e. "twisted" (u = ±1, w = 0) or "constant" (u = 0, w = ±1).
- That gives 3 bits per coset. The repo's code calls these labels (α, type, sign); see `scripts/paley_ansatz.py`.

**Reduced equations.** PAF at lag (h′, s) depends on s only through the class s = 0, s a square (QR) or s a non-square
(NQR). So the core condition is exactly 3 equations per h′ ∈ G′. For b ≥ 3:
- T(h′, ·) = −4b if h′ ∈ H′ (except T(0, 0) = 4m), and +4 otherwise;
- for b ≤ 1, T = −4b everywhere off (0, 0).

**Decomposition for p ≡ 3 mod 4** [M, derived from the correlation table and checked on 359 and on the found
solutions]. Write CC(f, g)(h) = Σ_x f(x) g(x + h), summed over the four blocks, and R_QR = Y + X, R_NQR = Y − X. Then

    R0 = PAF_α + (p−1)(PAF_u + PAF_w)
    Y  = −PAF_u + (p−2)·PAF_w + CC(α,w) + CC(w,α)
    X  = CC(α−w, u) − CC(u, α−w)

and R0 + (p−1)·Y = PAF(y), where y = α + (p−1)·w is the vector of Z_p-coset sums (y ∈ {±1, ±(p−2), ±p}). The core
condition is therefore equivalent to:
- **Stage 1 (compressed core on G′):**
  Σ PAF_y = T₀ + (p−1)·T_s, i.e. 4m − 4b(p−1) at 0, −4bp on H′ ∖ 0, and 4p off H′.
  For b = 1 it is −4p everywhere off 0.
- **Stage 2 (twisted signs u, given y):** u = ±1 on the twisted cosets, with
  Σ PAF_u = (p−2)·PAF_w + CC(α,w) + CC(w,α) − T_s, and the linear condition X = 0.

**Exact invariants.**
- h′ = 0: (p−1)·K + 2·Σ_const α·w = 4n′ − 4b, where K = #constant cosets over all four blocks and n′ = |G′|.
- Summing over h′ (compression onto Z_p), per block with A = Σα, U = Σu, W = Σw:
  Σ_k (−U_k² + (p−2)W_k² + 2A_kW_k) = Σ_{h′} T_s(h′).

For p ≡ 1 mod 4 the cross terms change sign behaviour (χ is even). Derive the analogue if you use such p.

## 3. Anatomy of the published 359 core [M]

v = 359, b = 7, c = 8, h = 44. G′ = Z₂² × Z₈ (coordinates (a, b, e)), H′ = ⟨(0, 0, 2)⟩ of order 4, p = 11, K = 12.
The core is in `data/core_v359_b7.txt` (blocks A–D on 2×2×8×11, last index fastest).

**Every component is separately two-valued off H′:**

| quantity | value off H′ | values on H′ = (0, 2, 4, 6) |
|---|---|---|
| Σ PAF_u | −4 | [116, 0, −4, 0] |
| Σ PAF_α | 44 = 4p | [128, 12, 52, 12] |
| Σ PAF_w | 0 | [12, −4, −4, −4] |
| Σ CC(α,w) + CC(w,α) | 0 | [−20, 8, 4, 8] |

X_αu = −X_uw ≠ 0 (they cancel lag by lag).

**w is algebraic:**
- Block A: w(a, b, 0) = s(a, b) and w(a, b, 6) = −s(a, b), with s = (+, −, −, −) bent on Z₂².
- Block B: w = (−, +) at (1, 0, 1) and (1, 0, 5), which differ by (0, 0, 4) ∈ H′.
- Block D: w = (−, +) at (0, 0, 0) and (0, 0, 4).
- Block C: no constant cosets.
- 11 of the 12 constant cosets have w = −α.

**α:** the −1 sets have sizes 9, 7, 1 and 14. Block C's α has a single −1, so its PAF is 28 everywhere. Block D's α has
PAF ≈ 0 off H′ and ±16 on H′.

**u ≈ a b = 1 core over G′ relative to H′:** Σ PAF_u = −4 off H′.

## 4. What is already established (do not redo)

**Cost of 359 in the ansatz** [E] (`scripts/paley_bge3.py estimate`):
- 384 free bits, 303.6 constraint bits, log₂E = +67.7.
- Given the published y, stage 2 alone has 116 free bits and 119.6 constraint bits (log₂E = −3.6).
- So the authors solved an approximately 300-bit system, either by design or with large compute.

**Our solvers** (`scripts/paley_ls.py` tabu and `--sa 1` simulated annealing, `scripts/paley_staged.py`):
- They solve bordered instances with n′ ≤ 10 in seconds: v = 15, 27, 31, 39, 47, 59 (labels in `logs/paley_ls_v*.json`).
- They fail at n′ = 12–28 (v = 51, 63, 65, 87), even though log₂E ≈ +60 there.
- Staged SA on 359 stalls in stage 1.
- Brute search is not the route.

**Exact screens** (`scripts/paley_bge3.py`, `scripts/paley_ansatz.py`):
- h0 = within-coset equations plus row sums.
- qnorm = compression onto G′/H′.
- For v = 63/65/87, the Z_p and Q-compressions both pass.

**Open-order options (p ≡ 3 mod 4 unless noted):**
- b = 1: these pass the exact screen.
  - G′ × Z_p: 599 Z₂₆×Z₂₃ (264 bits, log₂E +28; forced K = 5, all with w = −α); 571 Z₃₀×Z₁₉ (298, +41; K ∈ {6, 7});
    647 Z₃₄×Z₁₉; 573 Z₅₂×Z₁₁; 659 Z₉₄×Z₇; 643 Z₂₁₄×Z₃.
  - Also 647 Z₃₈×Z₁₇ and 573 Z₄₄×Z₁₃, with p ≡ 1 mod 4.
  - Anti-symmetric multiplier slices are impossible for b = 1 [M]. The best symmetric multiplier slices have
    log₂E = −3 to −10.
- b ≥ 3: dead are every option with p ≥ 17 and the c = 4 options for 647 and 739 (p = 23).
  - Undecided (screen running): 599 c = 6 p = 11 (n′ = 54, K ∈ 17–24); 669 c = 10 p = 11 (n′ = 60, K ∈ 17–25);
    719 c = 16 p = 11 (n′ = 64); 643 c = 28 p = 11 (n′ = 56); 647 c = 54 p = 11; and options with p ∈ {3, 7} or with
    p = 5, 13 (≡ 1 mod 4: 523 c = 4, 573 c = 14).
  - Full list: `logs/paley_bge3_jobs.txt`.

**The rest of the landscape:**
- `report/unattempted_2000_3000.md` covers the T-matrix routes, generalized borders and the 523 k-tree.
- `report/theory_arm_2000_3000.md` is a previous theory session; its composition-closure results stand.

## 5. Questions (reorder freely)

1. **Iterated Paley.** Several b = 1 options have G′ = Z₂ × Z_q or Z₂·q × …, with q prime: Z₂₆ = Z₂×Z₁₃,
   Z₃₄ = Z₂×Z₁₇, Z₃₈ = Z₂×Z₁₉, Z₉₄ = Z₂×Z₄₇, Z₂₁₄ = Z₂×Z₁₀₇.
   - Apply the ansatz again on the second prime, with x a combination of δ, χ, ε on Z_q tensor (δ, χ, ε) on Z_p:
     9-term coset types, a few bits per Z₂×(class) cell.
   - Derive the reduced equations exactly and count free versus constraint bits.
   - Does a doubly-twisted 599 = 1 + 2·13·23 (or 571 = 1 + 2·3·5·19, or 647 = 1 + 2·17·19) slice fall below about
     2^45, with exact screens passing?
   - Give positive controls at small orders.
2. **Designing w and α.** 359's w is a bent-function / H′-coset gadget with Σ PAF_w = 0 off H′.
   - Characterise the w-gadgets for general (G′, H′, K, p), including the forced K and Σαw.
   - Then: given w, is the α-stage (a ±1 quadruple with Σ PAF_α = 4p off H′, plus linear conditions from
     CC(α, w)) an instance of a known object: a relative or divisible difference family, a supplementary difference
     set relative to H′, or a cyclotomic family?
   - Construct it algebraically where you can.
3. **Designing u.** Stage 2 is Σ PAF_u = −4 off H′ (in 359), plus the linear condition X = 0.
   - This is a relative b = 1 core over G′ with a forced zero pattern.
   - Are there algebraic families (Paley or Szekeres on a prime-power factor of G′, product of characters, Golay or
     Turyn-type) that satisfy the quadratic part, leaving only X = 0 (linear in u) and a few H′ equations to fix by
     small search?
   - What freedom (automorphisms of G′, translations, signs) can be used to satisfy X = 0?
4. **Explain 359.** Find a design route that reproduces a 359 core (a new one, or the published one up to
   equivalence) at a total cost of at most 2^45. This is the positive control for any route you propose for an open
   order.
5. **b ≥ 3 open options with p = 11** (599 c = 6, 669 c = 10, 719 c = 16, 643 c = 28) mirror 359 most closely.
   - If the screens pass, which design transfers?
   - Note the border must match the compression for c ≥ 6. Address this too: which compressions does a c = 6 or
     c = 10 border admit?

## 6. Tools and conventions in the repo

- `scripts/paley_ansatz.py`:
  - `LABELS` / `LC` (label → (α, u, w)), `corr_table(p)`, `reduced_paf`, `exact_h0_screen_fast`;
  - symmetric-slice estimators;
  - `search` (exhaustive, for small n′).
- `scripts/paley_bge3.py`:
  - general abelian G′ and H′ (`AG`, `index_c_subgroups`);
  - `control` (the 359 check), `options v…`, `screen1 v c p`, `estimate v c p mods…`.
- `scripts/paley_ls.py`:
  - `Prob(v, b, p, mods, Hgen)` builds the reduced target;
  - `verify(P, lab)` expands to the full G′ × Z_p and checks Σ PAF by FFT, which is the acceptance test;
  - `search` (tabu) and `anneal` (SA).
- `scripts/paley_staged.py`:
  - stage 1 and stage 2 annealers;
  - `stage2_setup(P, y)` returns α, w, β = α − w and the u-target.
- A core is only a step. A Hadamard matrix needs the border as well; use `scripts/bge3_border_check.py` and the 347/359
  templates. For b = 1 the border is the standard bordered GS construction.
- Python with numpy, scipy and numba is available. There is no SageMath or GAP. nauty/Traces can be built from the
  pynauty sdist (see `theory/` from the previous session).
