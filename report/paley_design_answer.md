# Algebraic design inside the Paley-coset ansatz: exact reformulation, new obstructions, and the irreducible core

Date: 2026-10-09. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_paley_design.md`.
Scripts: `theory/paley_*.py`. Logs: `theory/logs/paley_*`, `theory/logs/stage1_*`.

**Tags.**
- **[M]** verified by exact computation here; the script is named.
- **[P]** proved here; the proof is in the text.
- **[E]** Gaussian estimate (same conventions as `scripts/paley_bge3.py`, recalibrated on 359 below).
- **[L]** literature recollection, not re-checked.

---

## 0. Bottom line

1. **(A) No Hadamard matrix of an open order was constructed.**

2. **(B) No design brings an open-order option to ≤ 2^45 operations.**
   - There is an exact reason why the natural designs fail (§2, §6).
   - After every exact reduction found here, the problem is a "stage 2": a critically constrained complex-complementary
     ±1 problem on G′ with about 4n′ free bits and about n′ equations.
   - Every decomposition I could make still costs about 2^(2n′) for it: 2^55–2^58 per candidate at 359, more at the open
     orders.

3. **(C) New exact obstructions.**
   - **599 at b = 1 (G′ = Z₂₆, p = 23), the prompt's top b = 1 option, contains no Paley-coset core at all** (§4, [P]).
   - The same "w-magnitude" screen also kills the b = 1 Paley options of 523 (p = 29), 537, 643, 659 (p = 47), 669, 719 and
     739 (p = 41), and many b ≥ 3 options [M].
   - The doubly-twisted ("iterated Paley") b = 1 slices of 599, 647, 659 and 643 are exactly empty, killed by their compression
     onto G″ [M].
   - The doubly-twisted slices of 571, 573 and 599 (c = 6) are predicted empty (log₂E = −11 to −46) [E].
   - Every "X-killing" design is predicted empty or arithmetically impossible: even slices (log₂E = −22 at 359, −50 to −59
     at 571/647/573), swap pairs, convolution designs and odd/odd slices (§6).

4. **New structure [P, M].**
   - The Paley-coset ansatz splits exactly into two independent quadratic systems on G′:
     - stage 1, a rational GS-type system for the Z_p-coset sums y;
     - stage 2, a GS system over the ring Z[(1+√−p)/2] for Z = β + √−p·u.

     They are coupled only through the ±1 entry constraint.
   - Consequences:
     - the norm identity (forced K);
     - the w-magnitude screen;
     - an exact description of 359's "separately two-valued" components and of its w-gadget;
     - the fact that the border-forced compression pins α to about 50 bits at 359.

5. **Positive controls [M].**
   - Exact stage 1 at 359, given the published compression and gadget, enumerates all 4096 stage-1 solutions in 25 min;
     the published α is among them.
   - The published 359 border, inflated, serves the whole c = 8, p = 11 family. For v = 183 exactly 2 single-deviation
     compressions are admissible.
   - The exact two-stage pipeline at v = 183: see §7 for the outcome.

---

## 1. Notation

- G = G′ × Z_p with p prime, p ≡ 3 (mod 4), n′ = |G′|, m = n′p, v = m + b.
- Each block, on each Z_p-coset, is x_k(g′, ·) = α·δ + u·χ + w·ε, with exactly one of u, w nonzero.
- Twisted cosets (w = 0) have entries ±1 in α and u. Constant cosets (u = 0) come in two types:
  - **J′**: w = −α;
  - **J**: w = +α.

  K = N₋ + N₊ is the number of constant cosets; N₋ counts the J′ cosets and N₊ the J cosets.
- Derived sequences on G′:
  - β = α − w;
  - y = α + (p−1)w, the Z_p-coset sums;
  - Z = β + √−p·u;
  - Z′ = Z/2 = π·u + e, with π = (1+√−p)/2, N(π) = (p+1)/4, and e = (α − u − w)/2 ∈ {0, ±1}.
- For b ≥ 3: H′ ≤ G′ has index c = b + 1, and Q = G′/H′. The target is T = 4m at 0, −4b on H∖0 and +4 off H, with
  H = H′ × Z_p.
- T₀ and T_s are the targets at lags (h′, 0) and (h′, s ≠ 0). T_J = T₀ + (p−1)T_s.

## 2. The exact reformulation [P, M]

**Theorem 1 (field split).**
- Let A, U, W ∈ Z[G′], and let X = A⊗δ + U⊗χ + W⊗ε ∈ Z[G′ × Z_p]. Write * for inversion.
- For p ≡ 3 (mod 4): χ* = −χ, χ² = Z_p − p, χ·Z_p = 0 and ε = Z_p − 1. Hence

      X X* = (ββ* + p UU*) ⊗ 1  +  (Uβ* − βU*) ⊗ χ  +  (AW* + WA* − UU* + (p−2)WW*) ⊗ Z_p,   β = A − W.

- Since 1, χ, Z_p are linearly independent, Σ_k X_k X_k* = T is equivalent to three equations in Z[G′]:
  - **(S2 real)** Σ_k (β_k β_k* + p U_k U_k*) = 4v. For b ≥ 3: T₀ − T_s = 4v at 0 and 0 elsewhere, with the same right
    side for b = 0 and b = 1.
  - **(S2 imaginary, "X = 0")** Σ_k (U_k β_k* − β_k U_k*) = 0.
  - **(S1)** Σ_k y_k y_k* = T_J. This is p × (the Z_p-coefficient) plus (S2 real), using (p−1)² = p² − 2p + 1.
- Equivalently, in characters ψ of G′:
  - S2: Σ_k |β̂_k(ψ) + i√p·û_k(ψ)|² = 4v for every ψ.
  - S1: Σ_k |ŷ_k(ψ)|² = 4v off H′^⊥, 4b on H′^⊥∖1 (b ≥ 3), and Σ r_k² = 4v at ψ = 1.
- In the ring O = Z[π]: Σ_k Z′_k Z′_k^† = v·δ, where † is inversion composed with complex conjugation.
  - Entries are ±π, ±π̄ (twisted, norm (p+1)/4), ±1 (J′) or 0 (J).
- This is the algebra decomposition span{I, Q_p, J} ⊗ Q ≅ Q × Q(√−p) of the Paley tournament scheme, applied blockwise.
- **Verified [M] on the published 359 core** (`theory/paley_algebra.py check359`):
  - all three equations and the complex form hold exactly;
  - e ∈ {0, ±1};
  - q·N_tw + N₋ = 3·116 + 11 = 359.

**Corollary 1 (norm identity) [P].** The coefficient at 0 of (S2 real) gives (p+1)·K − 4·N₋ = 4n′ − 4b.
- At 359: K ∈ {9, …, 12} with N₋ = 3K − 25. The published core has K = 12, N₋ = 11.
- The 116/128 "anti-invariance" that the forensics found (p ≈ 10⁻²¹) is therefore forced by arithmetic, not a design
  choice.

**Corollary 2 (stage structure) [P].**
- S1 involves only (α, w).
- Given y, stage 2 is: ΣPAF(u) = (4vδ − ΣPAF(β))/p, which is quadratic in u, plus X = 0, which is linear in u.
- The Y-equation reads exactly ΣPAF(u) = W − T_s with W := Σ_k [(p−2)PAF(w_k) + CC(α_k, w_k) + CC(w_k, α_k)].
- So **"ΣPAF_u is two-valued off H′" ⟺ W vanishes off H′**. The reduced equations alone do not force this.
- In 359 it holds through a visible cancellation [M]:
  - W_D ≡ 0 identically;
  - W_A = −W_B lag by lag off H′.

**Corollary 3 (compression pinning) [P].**
- S1 at ψ ∈ H′^⊥ gives Σ_k |ŷ_{Q,k}(ξ)|² = 4b on Q̂∖1. So every block's Q-compression is "constant + O(1)",
  with Σ_k Σ_q (y_Q − mean)² = 4b(c−1)/c.
- If w is balanced inside H′-cosets (forced whenever the compression values are below p − 1 in size), then
  y_Q = β_Q = α_Q. This holds in 359 [M].
- So the border fixes α's coset sums: **α carries only 50.7 bits at 359** (`theory/paley_compfix.py`, exact count).

**Level decomposition [P].**
- For H′ = Z₄ = ⟨t⟩, characters split by λ = ψ(t) ∈ {1, −1, ±i}.
- S1, S2 and X at level λ involve only the λ-components of the four values on each H′-coset.
- In 359 every component has a level-constant spectrum [M]:
  - Σ_k |û_k|² = 0, 128, 112, 120 at ψ = 1 and at λ = 1 (ψ ≠ 1), −1, ±i;
  - Σ_k |β̂_k|² = 1436, 28, 204, 116 at the same levels.
- It is still not a Davis–Jedwab building set: the λ-components spread over several Q-cosets
  (scratch analysis, logged in this report only).

## 3. Anatomy of 359 in the new coordinates [M]

These are all computed by `theory/paley_algebra.py` and `theory/paley_compfix.py`.

- **Compressions.** In every block, y_Q = β_Q = α_Q, because w is balanced in every H′-coset.
  - The pattern is constant + one deviation at q* = (1,1,1).
  - Constant parts (2, 2, 4, 0) and deviations (−2, +2, −2, +4); Σd² = 28 = 4b.
  - Row sums (14, 18, 30, 4); every Σu_k = Σw_k = 0.
- **w-gadget** = (bent function on Z₂²) ⊗ (pair δ₀ − δ_{3t}) in block A, plus single pairs at difference 2t in B and D.
  - ΣPAF_w = 12, −4, −4, −4 on H′ and 0 off H′.
  - General form (§5.1): a complementary family on a subset of Q, tensored with balanced patterns on H′.
- **Low-weight α.** α_C = 1 − 2δ (one −1), and α_D has PAF zero outside ⟨(1,1,0), (0,0,2)⟩. Both are forced by the
  compression, not designed: with only 4 cells per coset, coset sums near 2√v/c = 4.7 leave 1–2 minus signs per coset.
- **No hidden block relations.** I searched all 384 automorphisms of Z₂²×Z₈, the translations, signs and α↔u swaps
  (scratch analysis). The best is α_A ↔ α_B at 30/32 agreement, which is expected for two low-weight sequences.
  The u's reach at most 26/32.
- **u is not a b = 1 core with holes.** Filling the 12 holes in all 2¹² ways leaves Σ|ΣPAF_u + 4| ≥ 48 (scratch
  analysis).

## 4. A new exact obstruction: the w-magnitude screen [P, M]

**Lemma 4.**
- For every h′ ≠ 0, S1 reads (p−1)²·s(h′) + a(h′) + (p−1)·c(h′) = T_J(h′), where:
  - s = Σ_k PAF(w_k)(h′);
  - a = Σ_k PAF(α_k)(h′) ∈ [−4n′, 4n′], with a ≡ 0 (mod 4) (the PAF of a full ±1 sequence is ≡ n′ (mod 4));
  - c = Σ_k [CC(α_k,w_k) + CC(w_k,α_k)](h′) ∈ [−2K, 2K], even.
- Summing PAF(w) over G′ gives the exact identity K + Σ_{h′≠0} s(h′) = Σ_k t_k², with t_k = Σ_g w_k(g).
- The row sums are r_k = Σα_k + (p−1)·t_k, with |Σα_k| ≤ n′.

**Theorem 5. 599 with b = 1 (G′ = Z₂₆, p = 23) has no Paley-coset core.**

Proof:
1. The h′ = 0 equation, Σy² = T_J(0) = 2304, gives 5N₋ + 6N₊ = 25, so K = 5 and N₋ = 5 (prompt: "forced K = 5").
2. For b = 1, T_J = −92 off 0. Then |484·s| ≤ 92 + 104 + 22·10 = 416 < 484, so s(h′) = 0 for all h′ ≠ 0.
3. Hence Σ_k t_k² = K = 5. As a sum of at most four squares, that forces some |t_k| = 2.
4. For b = 1 the row sums are (±2, 0, 0, 0). A block with t_k = ±2 then needs |Σα_k| ≥ 44 − 2 > 26 = n′.
   Contradiction. ∎

- This is checked mechanically by `theory/paley_wscreen.py b1` [M].
- The repo's h0 screen passed this option because it uses only h′ = 0 and row sums.

**The screen over all options** [M] (`theory/logs/paley_wscreen_b1.jsonl`, `…_bge3.jsonl`):
- **b = 1, dead:**
  - 523 (p = 29);
  - 537 (p = 67);
  - **599 (p = 23)**;
  - 643 (p = 107);
  - 659 (p = 47);
  - 669 (p = 167);
  - 719 (p = 359);
  - 739 (p = 41).

  Several of these were already dead at h′ = 0.
- **b = 1, alive:**
  - 571 (p = 19, n′ = 30);
  - 573 (p = 11, n′ = 52; p = 13, n′ = 44);
  - 599 (p = 13, n′ = 46);
  - 647 (p = 17, n′ = 38; p = 19, n′ = 34);
  - 659 (p = 7, n′ = 94).

  For 571 and 647 (p = 19) the survivors need ΣPAF_w(h′) = −1 at a few lags, with c(h′) ≥ 8, i.e. extreme α–w
  correlation.
- **b ≥ 3:**
  - Every option with n′ ≤ 30 and a large p dies.
  - The p = 11 options (599 c = 6, 669 c = 10, 719 c = 16, 643 c = 28) survive, because (p−1)² = 100 is too small for
    the bound to bite.
- Stage 1 uses only coset sums, so the screen applies verbatim to p ≡ 1 (mod 4).

## 5. Question 1: iterated Paley

### 5.1 Reduced equations [P, M]

- Take G = G″ × F_{q₁} × … × F_{q_r}, invariant under ∏ squares(F_{qᵢ}^×), acting factorwise.
- A block is one ±1 value per (g″, product orbit): 3^r values per element of G″, i.e. "9-term coset types" for r = 2.
- **Reduced equations.** For each lag (h″, S), with S a class vector:

      R(h″,S) = Σ_{g″} Σ_{O₁,O₂} x(g″,O₁) x(g″+h″,O₂) N_S[O₁,O₂],   N_S = ⊗ᵢ (order-2 cyclotomic numbers of F_{qᵢ}).

  That is |G″|·3^r equations, fewer after h ↔ −h.
- **Field split.** span{δ, χ, ε} ⊗ Q ≅ Q ⊕ Q(√q*), so the GS condition splits into 2^r "sum of norms" conditions over
  G″, one per subset of the factors: rational, Q(√q*), Q(√p*) and the biquadratic field.
- Validated against FFT for prime fields and for F₉, F₂₅, F₂₇ (`theory/paley_multi.py validate`).

### 5.2 Results [M, E]

| option | free bits | result | method |
|---|---|---|---|
| 599 b=1, Z₂ × F₁₃ × F₂₃ | 72 | **empty**: no compatible compression onto G″ = Z₂ | exhaustive join [M] |
| 647 b=1, Z₂ × F₁₇ × F₁₉ | 72 | **empty** (compression) | [M] |
| 659 b=1, Z₂ × F₇ × F₄₇ | 72 | **empty** (compression) | [M] |
| 643 b=1, Z₂ × F₃ × F₁₀₇ | 72 | **empty** (compression) | [M] |
| 571 b=1, Z₆ × F₅ × F₁₉ | 216 | log₂E = −10.8, 202 constraint bits | [E] |
| 573 b=1, Z₄ / Z₂² × F₁₁ × F₁₃ | 144 | log₂E = −24.5 / −45.6 | [E] |
| 599 b=5 (c=6), Z₆ × F₉ × F₁₁ | 216 | log₂E = −11.0 | [E] |

**Controls [M].**
- b = 0 at v = 35 (F₅×F₇) and v = 55 (F₅×F₁₁): 15,360 and 1,280 solutions, each verified by FFT.
- Planted targets over G″ = Z₂ (F₅×F₇): 49,152 solutions found. The join is sound.
- Natural b = 1 slices over Z₂ × F_q × F_p are all empty: v = 71, 111, 131, 155, 171, 183, 191
  (`theory/logs/paley_multi_controls.log`).

**Answer to Q1.**
- A doubly-twisted 599, 647, 659 or 643 slice is exactly empty.
- 571's slice has 216 bits, is predicted empty, and as a 4-list problem would cost about 2^67 even if it were not.
- No iterated-Paley slice falls below 2^45 with solutions.

## 6. Questions 2–3: designing w, α, u; and why the X-killing designs fail

### 6.1 w-gadgets (Q2) [P, M]

- A gadget with PAF_w supported on H′ is a sum over blocks of f_k ⊗ φ_k:
  - f_k is a ±1 function on a subset S_k ⊆ Q whose PAF vanishes off 0 (a bent or perfect function, or a single point);
  - φ_k is a balanced pattern on H′ (Σφ_k = 0, which keeps y_Q = α_Q).
- Then ΣPAF_w = Σ_k |S_k|·PAF(φ_k) on H′, and the H′-values must be about T_J|H′ / (p−1)². At 359 that is −308/100 ≈ −3.
- 359 is the minimal instance: |S| = 4 (bent on Z₂²) with a pair at difference t, and two singletons with pairs at 2t,
  giving (−4, −4, −4).
- At b = 1 with large p, Lemma 4 forces the gadget to be a signed difference family with zero off-peak, which is
  impossible at 599 (Theorem 5).
- Randomly placed balanced gadgets miss S1 on H′ by thousands of bits (`theory/paley_option_cost.py`, [E]). The gadget
  must be designed; it costs few bits once the family is fixed.

### 6.2 The α-stage is a relative SDS (Q2)

- In separable solutions, ΣPAF_α = 4p off H′ (prescribed on H′), i.e. the −1 sets form a 4-{n′; k₁…k₄; λ} family
  **relative to H′**.
- It is accompanied by a linear condition from W = 0.
- At 359 the border pins α to 50.7 bits, so this stage needs no algebra: **exact meet-in-the-middle solves it**
  (`theory/paley_stage1.py`, §7).
- I found no classical family producing the full α-quadruple. α_C and α_D are compression-forced degenerate cases, not
  instances of one.

### 6.3 The u-stage (Q3)

- ΣPAF_u = −T_s + W is a b = 1-type core "relative to H′" with holes at the gadget.
- For b = 1 with large p it is a near-perfect quadruple: ΣPAF_u(h′) = 4 + c(h′) + (p−2)s(h′) ∈ [−4, 4].
- **Freedom for X = 0 [M].**
  - Reflecting one block (g ↦ −g), or applying the NQR multiplier on Z_p (u ↦ −u), preserves every PAF and negates X_k
    (checked on 359).
  - Translations preserve X_k.
  - So the only freedom is the 4 signs in Σ ε_k X_k = 0: 4 bits against 12 equations (about 50 bits) at 359.
- **Algebraic u-families tried, all fail:**
  - **Legendre products on a factor of G′** give too many holes (14 per block at 571, against K ≤ 7 in total).
  - **QR_q-invariant u** is the doubly-twisted slice, empty (§5).
  - **Kronecker bent ⊗ Golay** has ΣPAF_u = 0 off the Z₂² slice instead of −4 [P].
  - **"b = 1 core with holes"** is not satisfied by the published u (§3).

### 6.4 Designs that make X ≡ 0 automatically, and why each fails

| design | mechanism | result |
|---|---|---|
| even slice (α, u, w all even on G′) | Uβ* self-adjoint | log₂E = −21.5 at 359; −49.8 (571); −58.9 (647); −43.2 (573) [E] (`theory/paley_slices.py`) |
| co-symmetric blocks (center t) | same | −21.5 at 359 (center 2); center 1 infeasible [E] |
| swap pair (α_B = −c∘ρ, u_B = a∘ρ from block A = (a, c)) | X_A + X_B = 0, S1 : S2 = 1 : (p+1) | exact search: solutions at v = 55 (b = 0); none at the 13 b ≥ 3 / b = 1 cases tried (v = 29–95; at v = 57, b = 1, the (C, D) system has 2048 solutions but no Golay-type pair completes them) [M] (`theory/paley_swap.py`). The pair's contribution is in fixed ratio, so C, D must absorb S2 = −(p+1)P: arithmetically incompatible with b ≥ 3 targets (§6.5) |
| odd/odd (β, u both odd) | X_k = 0 | needs J-cosets at all self-inverse elements: K ≥ 8 at 571, 647 against K ≤ 7 / N₋ ≥ 7 [P] |
| convolution β_k = u_k ∗ s_k, ŝ real | X_k = 0 | requires (p+1)·(…) row-sum identities. At 359 the published row sums contradict B_k = U_k·S_k. In general about 44 equations on about 50 α-bits [P/E] |
| pure Legendre block (α = ±u) | X_k = 0 | ratio 1 : (p+1) again. All blocks of this type would need n′ = b [P] |

**Why they all fail.**
- X = 0 is about 12 equations (about 50 bits at 359).
- Each design that makes X vanish identically removes far more free bits: about 144 for the even slice at 359.
- The remaining S1/S2 equations keep their per-equation entropy, so log₂E drops by 20–60 bits.

### 6.5 The swap-pair ratio obstruction (sketch) [P]

- A swap pair contributes P := PAF(a) + PAF(c) to S1 and (p+1)P to S2.
- Off H′, S2 needs S2_CD = −(p+1)P and S1 needs S1_CD = 4p − P.
- Summing S1 over h′ gives Σ_{h′} P = (Σa)² + (Σc)² ≥ 0. With n′ = 32 this forces P ≥ −12 on average, while C, D can
  supply S2 values of size about 10 only.
- The exact small-case searches (13 cases with b ≥ 3 or b = 1) found no (C, D) solution except at v = 57.

## 7. Question 4: positive control at 359

**Stage 1 is solved exactly [M]** (`theory/paley_stage1.py check359`):
- Inputs: the published compression and gadget.
- Per-block α lists: 98,304, 16,384, 4 and 279,936 (2^50.7 combinations).
- A hashed meet-in-the-middle finds **all 4096 stage-1 solutions in 25 min** on one core.
- The published α is among them, and every solution has N₋ = 11, as the norm identity requires.
- The solutions are per-block symmetry images: 4, 16, 4 and 64 distinct α per block.

**Stage 2 is the bottleneck [M, E]:**
- **Gaussian count** (`theory/paley_compfix.py`):
  - 116 u-bits against 125.9 constraint bits gives log₂E = −9.9 per y;
  - with compression and gadget fixed, the whole system has log₂E = −7.9.
  - So with the border-forced compression and the published gadget, the published core is essentially the only one.
- **Exact Q-level count** (`theory/paley_stage2_levels.py`):
  - given the published y, **2^30.2 Q-level u-compressions** satisfy the Q-level part of S2 exactly;
  - 2^84.8 of the 2^116 u-configurations pass it.
  - The λ = −1 and λ = ±i levels must then remove about 95 bits through further joins of about 2^27 each.
  - About 2^57 per y.
- **Local search on stage 2 fails [M].** SA on u alone, given the published y, never finds the known solution in
  6 × 3·10⁸ steps (final cost 160–320).
- **The other exact decompositions all cost about 2^(N/2) = 2^55–2^58 per y**: Z₂²-character halves, G′/A₁ phases and
  block meet-in-the-middle.

**Cost of the design route at 359 [E]:**
- about 2^8 (compression, gadget) classes (needed because each class holds about 2^−8 solutions);
- × 2^57 per y in stage 2;
- ≈ **2^65**, against a target of ≤ 2^45.
- The authors' local search evidently did better in the full ansatz, where log₂E = +67.7 and solutions are dense.

**Positive control at v = 183** (the smallest member of 359's family: c = 8, b = 7, p = 11, G′ = Z₂²×Z₄, |H′| = 2) [M]:
- `theory/paley_c8family.py` rebuilds the published H(1436) exactly. It recomputes the corner from
  C = −K_big Core^T Kc_big/(32h) and reuses K and Kc, inflated by h.
- **Exactly 2 of the 64 single-deviation 183-compressions are accepted by the 359 border**: deviations (2, −2, 2, 4) with
  constants (−2, 2, −2, −2), and its negation.
- The norm identity forces K = 4 with N₋ = 3. Two pairs in two different blocks are needed, because a same-block pair of
  pairs leaves ±200 off H′.
- Stage 2 decouples exactly for |H′| = 2:
  - the Q-level involves only u_Q;
  - the λ = −1 level involves only the signs on the u_Q = 0 cosets.
- Measured per stage-1 solution: about 3·10⁶ Q-level solutions.
- RESULT: (pending; filled in below)

## 8. Question 5: the b ≥ 3, p = 11 options and their borders

| option | n′ | \|H′\| | K range (norm identity) | w-screen | border |
|---|---|---|---|---|---|
| 599 c = 6 (b = 5) | 54 | 9 | 17–24 | alive | **reuse 347's border (same b = 5)**: admissible compressions below |
| 669 c = 10 (b = 9) | 60 | 6 | 17–25 | alive | none known |
| 719 c = 16 (b = 15) | 64 | 4 | 17–24 | alive | none known. Sylvester H(64) ⊇ H₄⊗1₁₆ gives candidate K = H₄ ⊗ (rows 2–16 of H₁₆), corner untested |
| 643 c = 28 (b = 27) | 56 | 2 | 10–14 | alive | none known |

**Which design transfers.**
- Stage 1 transfers (exact meet-in-the-middle given compression and gadget).
- Nothing else does. Stage 2 has (4n′ − K) = 192–239 free bits and about n′ = 54–64 equations, a deficit of about −10 to
  −20, so roughly 2^100 or more per y.

**Border–compression matching (computable).**
- The corner is forced, C = −K Cb^T Kc^T/(4c), so a compression is admissible for a given border iff K Cb^T Kc^T ∈ ±4c
  entrywise.
- For c = 8 (359's border) the accepted single-deviation compressions are tabulated in §7 (183) and include 359's own.
- For c = 6 (347's border) at v = 599: see below. (pending)

## 9. What would change the picture

1. **A structured stage 2.** A family of ±1 quadruples u on a 2-group G′ with prescribed level-constant spectra and
   Σ_k Im(β̂ conj(û)) = 0 for a given β, i.e. a "complex building set" over Z[π].
   - None exists in the literature I know [L].
   - The exact level decomposition (§2) is the place to look: each level is a complementary-family problem over Q.
2. **Better local search in the full ansatz** rather than more structure.
   - The published 359 sits in a space with log₂E = +67.7 and about 300 constraint bits; the authors solved it.
   - The open b = 1 options 571 (+40.7, 298 bits) and 647 (+42.6, 343 bits) are in the same class.
   - The repo's tabu/SA fail already at n′ = 12–28, so the solver, not the ansatz, is the gap.
3. **A c = 16 border** (Sylvester-type, untested) would open 719 c = 16 to the same stage-1 machinery. Its stage 2 is
   still out of range.

## 10. Reproduce

```
python3 theory/paley_algebra.py check359            # Theorem 1 on the published core
python3 theory/paley_multi.py validate               # iterated-Paley reduction vs FFT (incl. F9, F25, F27)
python3 theory/paley_multi.py mitm 599 1 2 -- 13 23  # doubly-twisted 599: empty
python3 theory/paley_wscreen.py b1                   # w-magnitude screen (599 b=1 dead)
python3 theory/paley_wscreen.py bge3
python3 theory/paley_slices.py 359 7 11 2 2 8 --Hgen 0.0.2 --slice full,even
python3 theory/paley_compfix.py                      # 359 counts with compression + gadget fixed
python3 theory/paley_stage1.py check359              # exact stage 1 (25 min)
python3 theory/paley_stage2_levels.py                # exact Q-level count of stage 2 (3 min)
python3 theory/paley_c8family.py selftest            # 359 border reused for the c = 8 family
python3 theory/paley_pos183.py                       # two-stage positive control at v = 183
python3 theory/paley_swap.py search 55 0 11 5        # swap-pair design control
```
