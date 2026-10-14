# Alpöge-style builders (b = 0 and b ≥ 3) with the staged method

Date: 2026-10-14. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_alpoge_design.md`.
Scripts: `theory/ab_*.py`. Logs: `logs/ab/`.

**Tags.** [M] exact computation (script named); [P] proved here (proof in the text); [E] Gaussian estimate (estimator
and calibration named); [L] literature recollection, unchecked.

*This report is written while runs are in progress; §0 is the current state.*

---

## 0. Bottom line (current state)

1. **(A) No new open-order matrix yet.** (Status of the runs in §6–§7.)
2. **The staged method transfers to b ≥ 3, and the published 347 core lives in a staged slice smaller than 523's [M, E]**
   (§1, §6).
   - Stage 2 of every b ≥ 3 generalized-Paley (GP) slice is *identical in form* to b = 0 and b = 1: for every
     nontrivial character ψ of the staging prime Z_q, Σ_k |X̂_k(φ, ψ)|² = 4v for all φ on G' [P]. Only stage 1
     changes: its target is two-valued, 4m − 4b(q − 1) at 0, −4bq on H' \ 0, +4q off H' [P].
   - Verified exactly on the three published b ≥ 3 cores (479 along Z₇ and Z₁₇, 347 along Z₁₉ and Z₃, 359 along Z₁₁)
     and on the b = 0 core of H(892) [M] (`theory/ab_verify_published.py`).
   - 347's own slice (Z₂ × Z₃ × Z₃ × Z₁₉, index 6 on Z₁₉, −1 on the two Z₃'s; 10 stage-1 orbits) has stage 1 at
     log₂E ≈ +27.7 and, on 347's own y, stage 2 with H₂ = 164.4, C₂ = 117.8, C_S = 69.1, C(D|S) = 48.7, h_D = 32.5,
     **D_def = 16.2 bits** [E]. The 523 slice that built H(2092) had C₂ = 146.2 and D_def = 26.1.
3. **Borders: settled for c = 4, 6 and 8 at every open order, by transfer from the published matrices [P, M]** (§3).
   The forced corner depends only on the deviations of the Q-compression from its block means; each published frame
   (479: c = 4, 347: c = 6, 359: c = 8) borders every compression in a one-parameter-per-block class, and that class
   meets the row equation of every open order with the same c (all corners checked ±1).
4. **(C) Exact obstructions** (§4, §5): see the ranked table (§8).
5. **(D)** §8.

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
- **c ≥ 10:** no published template. Open (§8).

---

## 4. b ≥ 3 GP-type slices: sweep and staged estimates

*(to be completed: sweep table, ranking, exact stage-1 results, 863 c = 4 exhaustive stage 1)*

## 5. b = 0 at composite v

*(to be completed)*

## 6. Positive control: 347's slice end to end

*(running)*

## 7. Open-order candidates and costs

*(to be completed)*

## 8. Ranked table (D)

*(to be completed)*

## 9. Reproduce

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
python3 theory/ab_prod347.py --tag a1 --R 800 800 --K 4 128 --stop-on-core   # positive control (347 slice)
```
