# Algebraic design inside the Paley-coset ansatz: exact reformulation, new obstructions, and the irreducible core

Date: 2026-10-09. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_paley_design.md`.
Scripts: `theory/paley_*.py`. Logs: `theory/logs/paley_*`, `theory/logs/stage1_*`.

**Tags.**
- **[M]** verified by exact computation here; the script is named.
- **[P]** proved here; the proof is in the text.
- **[E]** Gaussian estimate (same conventions as `scripts/paley_bge3.py`, recalibrated on 359 below). It matched every
  exact control at p = 7, but overcounts badly at b = 1, p ≥ 11 (§7.2, §7.4).
- **[L]** literature recollection, not re-checked.

---

## 0. Bottom line

1. **(A) No Hadamard matrix of an open order was constructed.**

2. **(B) Not reached, but one route comes close.** No design is shown to bring an open-order option to ≤ 2^45
   operations.
   - There is an exact reason why the natural designs fail (§2, §6).
   - After every exact reduction found here, the problem is a "stage 2": a critically constrained complex-complementary
     ±1 problem on G′ with about 4n′ free bits and about n′ equations.
   - Every decomposition I could make still costs about 2^(2n′) for it: 2^55–2^58 per candidate at 359, more at the open
     orders.
   - **The closest route is the authors' own 419 design, transferred to 573 (§7.5). It is the best (B)-type result
     here: validated end-to-end at 419, costed at the edge of the 2^45 target at 573.**
     - The published 419 core turns out to be an exact Paley(11) core inside an order-3 multiplier slice.
     - The same slice at 573 (p = 11, G′ = Z₅₂) is predicted non-empty (log₂E = +23 [E]).
     - The exact slice pipeline rebuilds a verified 419 core from its gadget in 4 minutes.
     - The same design route, run from scratch, also finds a **new** H(1676) in 272 s.
     - At 573 it already produces stage-1 solutions: 392 y's from 19 of 61 designed gadgets, in about 45 minutes on
       three cores.
     - Ranking gadgets by predicted core yield (calibrated at 419) costs the first 573 core at about 2^48–2^49 [E]
       with plain stage-2 joins, about 2^3–2^4 above target. The stage-2 joins are independent per y, so the work is
       embarrassingly parallel. Stage 2 was not run at 573.

3. **(C) New exact obstructions.**
   - **599 at b = 1 (G′ = Z₂₆, p = 23), the prompt's top b = 1 option, contains no Paley-coset core at all** (§4, [P]).
   - The new "w-magnitude" screen kills it. Its other kills (the b = 1 Paley options of 523, 537, 643, 659 (p = 47), 669,
     719 and 739, and 25 b ≥ 3 options) already fail the norm identity, i.e. the repo's h0 screen [M].
   - The doubly-twisted ("iterated Paley") b = 1 slices of 599, 647, 659 and 643 are exactly empty, killed by their compression
     onto G″ [M].
   - The doubly-twisted slices of 571, 573 and 599 (c = 6) are predicted empty (log₂E = −11 to −46) [E].
   - Every "X-killing" design is predicted empty or arithmetically impossible: even slices (log₂E = −22 at 359, −50 to −59
     at 571/647/573, −15 at 659), swap pairs, convolution designs and odd/odd slices (§6).
   - **Lemma 6 (new, [P, M]):** S2 at the real characters of G′ fixes the row sums of u by stage-1 data, with a parity
     condition. Consequences:
     - for b = 1 the four blocks cannot all contain an odd number of constant cosets;
     - it is the exact reason why every stage-1 solution at v = 111 and v = 133 (p = 11) fails to complete;
     - it prunes gadgets, but kills no open option by itself (§4.1).

4. **New structure [P, M].**
   - The Paley-coset ansatz splits exactly into two independent quadratic systems on G′:
     - stage 1, a rational GS-type system for the Z_p-coset sums y;
     - stage 2, a GS system over the ring Z[(1+√−p)/2] for Z = β + √−p·u.

     They are coupled only through the ±1 entry constraint.
   - **The ansatz is a multiplier orbit slice:** that of the automorphism group 1 × QR_p of G′ × Z_p (§7.5). So
     "Paley + multiplier" designs are orbit slices of product groups.
   - **The published 419 core (b = 1) is an exact Paley(11) core on G′ = Z₃₈, invariant under ×7.** So both published
     Paley-coset cores, 359 and 419, use p = 11 [M].
   - Consequences:
     - the norm identity (forced K);
     - the w-magnitude screen;
     - an exact description of 359's "separately two-valued" components and of its w-gadget;
     - the fact that the border-forced compression pins α to about 50 bits at 359.

5. **Positive controls [M].**
   - The exact two-stage design pipeline (§7.2) constructs and verifies Hadamard matrices H(116), H(172), H(228),
     H(284), H(340) and H(396): b = 1, p = 7, n′ = 4–14 (`theory/paley_b1pipe.py`).
     - Each takes at most 30 s, except H(396), which takes 3 min.
     - H(396) lives on the non-cyclic group Z₁₄ × Z₇.
     - Every core is re-checked independently by `paley_b1pipe.py verify`.
   - At p = 11 it proves v = 67, 89 and 111 exactly empty, and v = 133 empty with at most 4 constant cosets per block,
     although the Gaussian estimates are +23 to +40.
     - At 67 and 89 stage 1 is already empty, because the gadget is rigid.
     - At 111 and 133 every stage-1 solution violates Lemma 6.
     - So the Gaussian model is not trustworthy for p ≥ 11 at b = 1. 599 at b = 1 shows the same pattern (prompt: +28,
       exactly empty here).
   - Exact stage 1 at 359, given the published compression and gadget, enumerates all 4096 stage-1 solutions in 25 min;
     the published α is among them.
   - The published 359 border, inflated, serves the whole c = 8, p = 11 family. For v = 183 exactly 2 single-deviation
     compressions are admissible. The c = 4 builder rebuilds a verified H(1916) from the 479 core with a border found
     automatically.
   - None of this reaches ≤ 2^45 at 359 or at an open order (§7.3).

6. **Where the difficulty sits at the open orders [E] (§7.4).**
   - The most slack-rich open option is 659 at b = 1 (G′ = Z₉₄, p = 7): log₂E up to +236 in the full ansatz.
   - But stage 2 given y has a deficit of 99–167 bits, so no staged method can reach it. It needs a joint solver, and the
     X-killing even slice is dead (−15).
   - For 571 and 647 (p = 19) stage 2 per y is near break-even, and all the cost sits in the rigid stage 1. That is the
     regime where the Gaussian model failed at every p ≥ 11 control.

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
  - S1: Σ_k |ŷ_k(ψ)|² = 4v off H′^⊥, 4b on H′^⊥∖1 (b ≥ 3), and Σ r_k² = T̂_J(1) at ψ = 1 (4v for b ≥ 3, 4 for b = 1).
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

**The screen over all options** [M] (`theory/logs/paley_wscreen_b1.jsonl`, `…_bge3.jsonl`).
**The only new kill is 599 (p = 23).** Every other dead option has no (K, N₋, N₊) at all, i.e. it fails the norm identity,
which is part of the repo's h0 screen.
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

### 4.1 Lemma 6: S2 at the real characters fixes the u row sums [P, M]

**Lemma 6.**
- Let ψ be a character of G′ with ψ² = 1, and write f̂(ψ) = Σ_g ψ(g) f(g).
- On every core of the ansatz (any b):

      Σ_k û_k(ψ)² = (4v − T̂_J(ψ))/p + 2 Σ_k ŷ_k(ψ) ŵ_k(ψ) − p Σ_k ŵ_k(ψ)²,

  with û_k(ψ) ≡ |T_k| (mod 2) and |û_k(ψ)| ≤ |T_k|, where T_k is the set of twisted cosets of block k (|T_k| = n′ − K_k).
- At ψ = 1, writing ρ_k = Σu_k, this reads Σ_k (r_k − p t_k)² + p Σ_k ρ_k² = 4v.

Proof:
1. S2 (real) at ψ gives Σ_k β̂_k² + p Σ_k û_k² = 4v.
2. Substitute β = α − w = y − p·w, and use S1 at ψ: Σ_k ŷ_k² = T̂_J(ψ).
3. The parity holds because û_k(ψ) is a ±1-sum over T_k. ∎

The proof uses only S1 and the 1-component of XX*, which is ββ* + pUU* for either sign of χ* = ±χ. So Lemma 6 also
holds verbatim for p ≡ 1 (mod 4).

**Corollary 6.1 (b = 1).**
- At ψ = 1, Σ_k ρ_k² = 4n′ + 4t₀ − p Σ_k t_k², where block 0 is the block with row sum 2.
- Since t_k ≡ K_k ≡ ρ_k (mod 2), reducing mod 8 shows that **the four blocks cannot all have an odd number of constant
  cosets**:
  - if all four were odd, both sides would be ≡ 4 (mod 8) only if n′ + t₀ − p is odd;
  - but n′ is even and t₀, p are odd, so n′ + t₀ − p is even.

**Verification and use [M].**
- Verified exactly at every real character on the published 359 core and on all six b = 1 controls of §7.2
  (`theory/paley_lemma6.py`, `theory/logs/paley_lemma6.log`).
  - At 359, K_k = (8, 2, 0, 2) and t = ρ = 0.
- **It explains the p = 11 failures.**
  - At v = 111 the exact pipeline finds 8 stage-1 solutions. All have t = (1, ±1, ±1, ±1) and K_k = (1, 1, 1, 1), so
    Σρ² = 40 + 4 − 44 = 0 while every ρ_k is odd. Impossible.
  - At v = 133 the same holds: 8 stage-1 solutions, all with K_k = (1, 1, 1, 1), which Corollary 6.1 forbids.
  - The lattice test `theory/paley_s2lattice.py` independently finds the stage-2 target outside the congruence lattice
    for every one of these 16 stage-1 solutions (`theory/logs/paley_s2lattice_111.log`, `…_133.log`). The Gaussian
    model predicted 2^13.5 (111) and 2^15.9 (133) completions per y.
- **Screen.** Added to `theory/paley_wscreen.py` (with K_k ≡ t_k and |ρ_k| ≤ n′ − |t_k|). It changes no verdict for
  the open-order options screened (`theory/logs/paley_wscreen_*_rho.jsonl`): it removes witnesses, not options.
- **Pipeline filter.** In `theory/paley_b1pipe.py`, ψ = 1 filters gadgets and ψ = (−1)^g filters stage-1 solutions.
  - At 111 it removes the 8 gadget shapes behind those stage-1 solutions, and at 133 the 16 shapes behind its 8, so
    stage 2 is never entered.
  - At p = 7 it cuts nothing; every control still completes.

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

## 7. Question 4: positive controls, and why 359 stays out of reach

### 7.1 Stage 1 at 359 is solved exactly [M]

`theory/paley_stage1.py check359`:
- Inputs: the published compression and gadget.
- Per-block α lists: 98,304, 16,384, 4 and 279,936 (2^50.7 combinations).
- A hashed meet-in-the-middle finds **all 4096 stage-1 solutions in 25 min** on one core.
- The published α is among them, and every solution has N₋ = 11.
- The solutions are per-block symmetry images: 4, 16, 4 and 64 distinct α per block. So up to symmetry, the published
  compression and gadget admit essentially one y.

### 7.2 The exact two-stage pipeline and its positive controls [M]

**Pipeline.**
1. The norm identity gives (K, N₋, N₊).
2. Enumerate w-gadgets up to per-block translation and reflection, filtered by Lemma 4.
   - Translation preserves every PAF and X_k.
   - Reflection is combined with u ↦ −u (the NQR multiplier), which also preserves X_k (checked on 359).
3. **Stage 1.** Take every α with the required row sums (b = 1) or the border-forced compression (b ≥ 3), and run an
   exact 2+2 join on S1.
4. **Stage 2.** Exact 2+2 join on (ΣPAF(β) + pΣPAF(u), X) = (4vδ, 0).
   - For |H′| = 2 the levels decouple exactly: the Q-level sees only the coset sums u_Q, and the λ = −1 level sees only
     the signs on the cosets with u_Q = 0 (`theory/paley_pos183.py`).
5. **Completion.**
   - b = 1: the bordered GS array with border rows and columns S, T ∈ H(4), constant on blocks. The corner is forced,
     Cn = −S Σᵀ T/4, with Σ the block row-sum matrix. Checked by rebuilding H(668) from the published 167 core
     (`paley_b1pipe.py verifycyclic data/core_v167_b1.txt`, `theory/logs/paley_b1pipe_verifycyclic_167.log`).
   - b ≥ 3, c = 4: `border_c4`. Checked by rebuilding H(1916) from the published 479 core with an automatically found
     border (`paley_pipeline.py check479`, 10 s, `theory/logs/paley_pipeline_check479.log`).
   - b ≥ 3, c = 8: the 359 border, inflated. Checked by rebuilding H(1436).

**Results** (`theory/paley_b1pipe.py`, `theory/logs/paley_b1pipe_*.log`; every H checked as H Hᵀ = N I in exact integers):

| v | b | p | n′ | result | time |
|---|---|---|---|---|---|
| 29 | 1 | 7 | 4 | **H(116)** | 0.7 s |
| 43 | 1 | 7 | 6 | **H(172)** | 0.8 s |
| 57 | 1 | 7 | 8 | **H(228)** | 1.4 s |
| 71 | 1 | 7 | 10 | **H(284)** | 2.6 s |
| 85 | 1 | 7 | 12 | **H(340)** | 27 s |
| 99 | 1 | 7 | 14 | **H(396)**, on the non-cyclic G = Z₁₄ × Z₇ | 191 s |
| 45 | 1 | 11 | 4 | exactly empty: no gadget has feasible row sums | 0.0 s |
| 67 | 1 | 11 | 6 | exactly empty (16 gadget classes, 0 stage-1 solutions); Gaussian log₂E = +23 (split estimator) | 0.0 s |
| 89 | 1 | 11 | 8 | exactly empty (33 classes, 0 stage-1 solutions); Gaussian log₂E = +28.7 | 0.1 s |
| 111 | 1 | 11 | 10 | exactly empty (1,700 classes). 8 stage-1 solutions, all killed by Lemma 6; Gaussian log₂E = +28 (split estimator) | 6.5 s |
| 133 | 1 | 11 | 12 | exactly empty with ≤ 4 constant cosets per block. Without Lemma 6: 24,856 classes, 8 stage-1 solutions, none completes. With Lemma 6: 24,792 classes, 0 stage-1 solutions (the 16 gadget shapes it removes carry all 8). Gaussian log₂E = +33 (split) / +40 (`paley_slices.py`) | 358 s |
| 177 | 1 | 11 | 16 | not run: the class count grows from 1,700 (111) to 24,856 (133), so 177 should have ≳ 10^5 classes [E], each with a 2^27-pair stage-1 join | — |
| 179 | 3 | 11 | 16 | no control possible with a balanced gadget (\|y_Q\| ≤ 4): Corollary 3 forces three single-bump blocks and one flat block (Q = Z₂²), with row sums in {±2, ±6, ±10, ±14} and {0, ±8, ±16}, and no four such squares sum to 4v = 716 [P, M: `paley_pipeline.py run179`]. The b ≥ 3 control was run at 183 instead | 0 s |
| 183 | 7 | 11 | 16 | 17 stage-1 orbit classes with the 359-border compression swept exactly, no core. About 2^21.5 Q-level solutions and 30,912 zero patterns each. The sweep then stopped (memory) | 9 min |

**Reading.**
- At p = 7 the norm identity gives 2K − N₋ = n′ − 1, so K is at least about n′/2. Stage 2 keeps slack, and the first
  stage-1 solution always completed.
- At p = 11 every order through n′ = 12 is empty, even though the Gaussian models predict 2^23–2^40 solutions
  (`paley_split_b1.py`: 67 +23, 111 +28, 133 +33; `paley_slices.py`: 89 +28.7, 133 +40). The models miss two exact
  effects:
  - the rigid gadget (Lemma 4), which empties stage 1 at 67 and 89;
  - Lemma 6, which blocks every stage-1 solution at 111 and 133.
- This is the same failure as at 599 (b = 1), where the Gaussian estimate is +28 and the ansatz is exactly empty
  (Theorem 5).
- No p = 11, b = 1 core has been seen at any order. The next orders, 155 and 177 (n′ = 14, 16), pass both screens but
  are beyond the present exact pipeline.

### 7.3 Why the pipeline does not reach 359 or an open order at ≤ 2^45

**Stage 2 is the wall.**
- The exact 2+2 join holds 2^(n′ − K/4) per list, i.e. about 2^(2n′) pairs. That is feasible here up to n′ ≈ 16.
- At 359 the Q-level count is exact: **2^30.2 Q-level solutions** given the published y, with 2^84.8 of 2^116
  u-configurations passing.
- Solving the remaining levels per Q-solution:
  - naively: about 2^27 each, so about 2^57 per y;
  - with the pair-sum / zero-pattern decoupling of §7.2 extended to Z₄: estimated 2^44–2^48 per y [E].
- Stage 2 has log₂E = −9.9 per y [E]. So about 2^10 distinct y-classes are needed, i.e. about 2^8 (compression, gadget)
  classes.
- Total: **2^54–2^58 at 359** [E].
- **SA on stage 2 fails [M]:** given the published y, SA never finds the known u in 6 × 3·10⁸ steps (final cost
  160–320).

**At the open orders it is worse.**
- b = 1 needs n′ ≥ 26: 571 (30), 647 (34, 38), 599 (46), 573 (44, 52), 659 (94).
- b ≥ 3, p = 11 needs n′ = 54–64.
- So the stage-2 joins are 2^60 or more, and Gaussian log₂E per y is negative wherever p ≥ 11.

### 7.4 Where the difficulty sits at the open b = 1 options [E]

**Method.** `theory/paley_split_b1.py` splits the Gaussian estimate (same density estimator as `paley_compfix.py`) for
each (K, N₋, N₊) allowed by the norm identity, with the exact row sums imposed:
- stage-1 log₂E: (gadget, J/J′ choice, w, α) placements minus the S1 bits;
- stage-2 log₂E per y: u bits minus the S2 and X bits, given one sampled y.

Logs: `theory/logs/paley_split_b1_*.log`.

| option | n′ | stage 1 log₂E | stage 2 per y | core log₂E (best K) | exact pipeline |
|---|---|---|---|---|---|
| 85, p = 7 (control) | 12 | +35 … +44 | +2 … +8 | +48 | core at the first stage-1 solution |
| 67, p = 11 (control) | 6 | +14 | +9 | +23 | stage 1 empty |
| 111, p = 11 (control) | 10 | +18 | +9 | +28 | stage 2 empty (Lemma 6) |
| 133, p = 11 (control) | 12 | +18 … +23 | +10 … +11 | +33 | stage 2 empty (Lemma 6) |
| **659, p = 7, Z₉₄** | 94 | +255 … +372 | **−99 … −167** | **+236** (K = 72) | — |
| 573, p = 11, Z₅₂ | 52 | +86 … +121 | −23 … −37 | +89 (K = 22) | — |
| 571, p = 19, Z₃₀ | 30 | +25 … +32 | +1 … +4 | +33 (K = 7) | — |
| 647, p = 19, Z₃₄ | 34 | +35 … +38 | −4 … 0 | +35 (K = 7) | — |

**Reading.**
- **659 at p = 7 is the most slack-rich open option** in the ansatz: +236 here, and +288.5 from the full-ansatz
  estimator `paley_slices.py` (`theory/logs/paley_slices_659_p7.log`). p = 7 is also the regime where the Gaussian
  model was right at every control (n′ = 4–14).
  - Stage 2 per y is 99–167 bits short. No staged method (fix y, then solve for u) can work: it needs about 2^100
    stage-1 solutions per core.
  - The X-killing even slice is dead (log₂E = −15.3) and the center slice is infeasible (`theory/paley_slices.py`).
  - What is needed is a joint solver on the 1128-bit, 140-equation system. That is the single most promising open-order
    target I found for a non-staged search. Nothing here reaches it.
- **571 and 647 at p = 19:** stage 2 per y is near break-even, so the whole cost is in stage 1.
  - Stage 1 at p ≥ 11 is the rigid-gadget regime, where the Gaussian model overcounted at every control (67, 89, 111,
    133, and 599 at p = 23).
  - So their +33 and +35 should be read as upper bounds.
- **573 at p = 11** has both problems: a rigid stage 1, and a stage 2 that is short by 23–37 bits.

**Scaling at p = 7 (b = 1, cyclic G′) [E].** Stage 2 per y, over the K range, falls steadily:

| n′ | 12 (v = 85) | 20 (141) | 30 (211) | 40 (281) | 60 (421) | 94 (659) |
|---|---|---|---|---|---|---|
| stage 2 per y | +2 … +8 | −8 … +5 | −26 … −3 | −39 … −15 | −88 … −43 | −167 … −99 |
| core log₂E (best K) | +48 | +73 | +100 | +126 | +171 | +236 |

- So the staged route is open only up to n′ ≈ 20–25.
- The exact 2+2 join for stage 2 costs about 2^(1.75 n′), since u has about 3.5 n′ bits. That is 2^24.5 at n′ = 14
  (H(396) took 3 min) and 2^35 at n′ = 20.
- So **the exact pipeline tops out near n′ ≈ 20 even where solutions are plentiful.**
- The joint annealer of the repo (`scripts/paley_ls.py --sa 1`, 10⁸ proposals per restart) is weaker still: it solves
  v = 57 (n′ = 8) but not v = 71 or 85 (n′ = 10, 12), where the exact pipeline succeeds (`theory/logs/paley_sa_p7_controls.log`).
- **659 needs n′ = 94.**

### 7.5 The published 419 core is a Paley(11) core in a multiplier slice, and 573 has the same slice [M, E]

**The Paley-coset ansatz is a multiplier orbit slice [P].**
- On G = G′ × Z_p, consider the automorphism group {id_G′ × (z ↦ az) : a a nonzero square mod p}. When G is cyclic
  this is the multiplier group {μ ≡ square (mod p), μ ≡ 1 (mod n′)}.
- A sequence is constant on its orbits iff on each Z_p-coset it takes one value at z = 0, one on the squares and one on
  the non-squares.
- Writing these three values as α, u + w and −u + w, that is exactly x = αδ + uχ + wε, with exactly one of u, w nonzero
  because the values are ±1.
- So the Paley-coset ansatz **is** the orbit slice of that order-(p−1)/2 group. Adding a multiplier μ on G′ gives the
  orbit slice of the product group: the labels (α, type, sign) become constant on the orbits of g ↦ μg.

**419 [M]** (`theory/paley_anatomy419.py`, `theory/logs/paley_anatomy419.log`):
- The published b = 1 core `data/core_v419_b1.txt` becomes, after translating blocks B, C, D by 208 (which preserves
  every PAF), an **exact Paley-coset core with p = 11 on G′ = Z₃₈**.
- Its labels are invariant under g ↦ 7g (order 3), so it lies in the "Paley(11) + mult:7" slice. That is the order-15
  multiplier slice the forensics identified (`report/alpoge_methods_forensics.md` §5).
- The forensics' "8 of 56 orbit pairs deviate from the Legendre twist" are its constant cosets.
- Data:
  - K = 16 (N₋ = 11, N₊ = 5) with K_k = (1, 13, 2, 0); the norm identity holds.
  - t = (−1, −1, 0, 0) and ρ = (1, 5, −6, 8).
  - Lemma 4 values s(h) ∈ {−2, −1, 0}.
  - S1, S2 and X = 0 hold exactly, and Lemma 6 holds at both real characters.
  - The gadget is a single J coset in A; in B, the J coset {0}, three J′ triples and one J triple; the J′ pair {0, 19}
    in C; and nothing in D.
- **So the two published Paley-coset cores, 359 (b = 7) and 419 (b = 1), both use p = 11.**

**Slice estimates [E]** (`theory/paley_slices.py`, `theory/logs/paley_slices_paleymult.log`):

| option | slice | orbits / block | free bits | constraint bits | log₂E |
|---|---|---|---|---|---|
| 419 (control, solved by its authors) | Paley(11) + ×7 (order 3) | 14 | 168 | 124 | +21.1 |
| 419 | full Paley(11) | 38 | 456 | 340 | +96.1 |
| **573**, p = 11, G′ = Z₅₂ | Paley(11) + ×29 (order 3 on Z₁₃) = the repo's d = 15 slice | 20 | 240 | 193 | **+23.2** |
| 573 | + ×5 (order 4) / + ×49 (order 6) | 16 / 12 | 192 / 144 | 181 / 144 | −13.0 / −23.7 |
| 571, p = 19, G′ = Z₃₀ | + ×7 / ×11 / ×19 | 12 / 20 / 18 | 144 / 240 / 216 | 261 / 232 / 202 | −142 / −15.9 / −10.3 |
| 647, p = 19, G′ = Z₃₄ | + ×13 (order 4) | 10 | 120 | 135 | −37.5 |

- **573 + ×29 is the only open-order analogue of 419's slice that is predicted non-empty.** Its log₂E is about the same
  as 419's (+23 against +21), but it has 69 more constraint bits.
- The repo had costed this slice as "d = 15: 193 bits, never run" (`report/unattempted_2000_3000.md` §3).

**Stage split inside the slices [E]** (`theory/paley_split_slice.py`, `theory/logs/paley_split_slice_*.log`):

| slice | stage 1 log₂E | u bits (stage 2) | stage 2 per y | core log₂E |
|---|---|---|---|---|
| 419, ×7 (K = 13–17) | +11 … +18 | 49–50 | −0.5 … −3.5 | +10 … +15 |
| 573, ×29 (K = 17–25) | +19 … +34 | 68–73 | −8.6 … −15.7 | +10 … +23 (best K = 21) |

**Exact slice pipeline at 419 [M]** (`theory/paley_slicepipe.py`):
- **The published gadget, held fixed (`given419`).** α and u are searched exactly in the slice: per-block α lists
  (495, 56, 924, 1848), an exact S1 join, Lemma 6 at both real characters, then an exact stage-2 join over u on the
  twisted orbits.
  - It produced **a verified H(1676) in 145–241 s** (two runs; the first needed 3 stage-1 solutions, the second 1).
  - The core is re-verified by `paley_b1pipe.py verify`. The command writes it to `data/core_v419_b1_paley.txt`; it is
    not committed, because the committed 419 core is the from-scratch one below.
  - This is the stage-1 + stage-2 analogue of the 359 control of §7.1, at a b = 1, p = 11 order.
- **Gadgets drawn at random in the slice:** 400 gadgets each at K = 15 and K = 16, all passing Lemma 4 and Lemma 6,
  give **0 stage-1 solutions** (`theory/logs/paley_slicepipe_v419_mu7_K1*_s*.log`).
- **"Designed" gadgets:** signed unions of ×7-orbits with |PAF_w| ≤ 2 in every block and s(h) ∈ {−2, …, 1}. The
  published gadget is one of them: its 13-coset block has PAF_w ∈ {−1, 0} off 0, a cyclotomic near-difference set of
  Z₃₈. But:
  - most such gadgets have |t_k| ≥ 2, which makes α infeasible;
  - the feasible ones have Gaussian stage-1 yields of 2^−3 to 2^−23, against **2^+6.9 for the published gadget**.
- **Yield-guided search, from scratch [M]** (`paley_slicepipe.py guided 419 11 38 7 15 8 400 4`,
  `theory/logs/paley_slicepipe_guided_v419_mu7_K15_s4.log`). Designed gadgets as above, restricted to |t_k| ≤ 1 and
  α-feasible, with the exact stages run only when the Gaussian stage-1 yield of that gadget exceeds 1:
  - 1,069 candidate gadgets were screened: 568 failed the s(h) filter, 436 were α-infeasible, 60 had low yield.
  - Their predicted yields have median 2^−26.9, 90th percentile 2^−2.2 and 99th percentile 2^+5.6.
  - **The 5th gadget run exactly (predicted 2^+6.6) gave a core: a new verified H(1676) in 272 s** on one core,
    gadget screening included.
  - It has K = 15 (the published core has K = 16) and is ×7-invariant, so it is a different core. It was not tested
    for Hadamard equivalence with the published matrix. It is in `data/core_v419_b1_paley_guided_K15.txt` and
    re-verified by `paley_b1pipe.py verify`.
  - For comparison, the repo's k-tree tool for this slice size (`kt_b1.py`) had only passed a planted 419 test, where
    the solution's bucket is given (`logs/kt_b1_planted_v419.log`).
  - **This is a complete positive control of a design route at a b = 1, p = 11 order:** gadget designed, α and u
    exact, border by closed form, no input from the published core.
- **Reading.**
  - Given a good gadget, the slice pipeline is fast: about 2^25 work per gadget at 419.
  - Gadget design was the bottleneck, as at 359 (§6.1). It is solved at 419 by three steps:
    - **design:** cyclotomic near-difference sets on the multiplier orbits, with |PAF_w| ≤ 2 and |t_k| ≤ 1;
    - **exact filters:** Lemma 4, Lemma 6 and α-feasibility;
    - **ranking:** each gadget's Gaussian stage-1 yield.
  - Random gadgets give nothing (0 in 800).

**Cost at 573 in the same design [E].**
- Per block the slice has 20 orbits (4 fixed points and 16 triples of ×29 on Z₅₂). With 0–4 constant orbits per
  block, α lists hold 2^11–2^16 (exact DP counts) and u lists 2^16–2^20.
- Given a good gadget:
  - stage 1 costs about 2^28–2^32 per gadget (exact 2+2 join);
  - stage 2 per y costs about 2^34–2^38 (2+2 join), less after splitting by the Lemma 6 classes.
- A priori, the slice split predicts a stage-2 shortfall of 9–16 bits per y, which gives 2^43–2^52 in stage 2. The
  measured gadgets below sit near the upper end.
- **Gadget search, measured at 573 [M]** (`paley_slicepipe.py guided1`, stage 1 only,
  `theory/logs/paley_slicepipe_guided_v573_mu29_K2*_s*.log`). The designed-gadget generator works unchanged on Z₅₂.
  - It takes about 1–2 minutes per exact stage-1 join, with α lists of 700–77,000 per block.
  - The first K = 22 gadget run (predicted 2^1.4) gave **8 stage-1 solutions** in 129 s, all passing Lemma 6 at both
    real characters (`theory/logs/paley_slicepipe_stage1_v573_mu29_K22_s1.npz`).
  - Over the three K = 21–23 runs: **61 gadgets were run exactly in about 45 minutes on three cores.**
    - The K = 21 run was stopped after 22 minutes by the container memory limit, with three large joins running at
      once.
    - 19 gadgets produced stage-1 solutions: **392 y's in total**, all passing Lemma 6 at both real characters
      (`theory/logs/paley_slicepipe_stage1_v573_mu29_K2*_s*.npz`).
    - Per-y stage-2 predictions are 2^−11.9 to 2^−12.6 at K = 21, 2^−12.7 to 2^−15.8 at K = 22 and 2^−14.0 to
      2^−17.5 at K = 23.
  - The first two K = 21 gadgets gave none.
- **Stage 2 for those y's [E]** (`theory/paley_s2predict.py`, `theory/logs/paley_s2predict_573_K22.log`):
  - u-lists of 2^15, 2^16, 2^19 and 2^18 per block, so the exact 2+2 join is about 2^34 per y (not run here);
  - predicted yield 2^−15.5 per y, against the slice average of 2^−9 to 2^−16.
  - On this instance, stage 2 alone needs about 2^15 y's × 2^34 ≈ 2^49.
- **Ranking gadgets by predicted core yield [E]** (`theory/paley_coreyield.py`).
  - Score per gadget = Gaussian stage-1 yield + Gaussian stage-2 yield at fixed y, with α drawn from the gadget's exact
    lists.
  - Calibrated on the two 419 gadgets that gave cores: predicted core yield 2^+2.7 and 2^+2.9, with stage 2 at
    2^−3.4 and 2^−4.2 per y. At their actual y's the stage-2 prediction is 2^−4.3 and 2^−2.8, and both completed.
  - With about 2^34 per y for stage 2, the cost to a first core is about 2^(34 − s₂) when stage 1 is plentiful
    (s₂ = per-y stage-2 log-yield).
  - **At 573 [E]** (`theory/logs/paley_coreyield_573.log`, 60 designed gadgets per K):
    - The best predicted core yields are 2^−9.5 (K = 21), 2^−12.7 (K = 22) and 2^−13.8 (K = 23). The median
      per-y stage-2 yields are 2^−14 to 2^−16.
    - The best K = 21 gadget has stage 1 at 2^+2.8, stage 2 at 2^−12.3 per y, and u-lists of 2^15–2^20, so its stage-2
      join is about 2^36.6.
    - Hence about 2^9.5 gadgets × 2^2.8 y × 2^36.6 ≈ **2^49 to a first core**.
  - The exact stage-1 runs agree: 8, 22 and 16 y's from three gadgets, with per-y stage-2 predictions of 2^−15.5,
    2^−15.5 and 2^−17.5 (`theory/logs/paley_slicepipe_guided_v573_mu29_K2*_s*.log`).
- **Practical cost [E].**
  - Stage 1 is cheap: stage-1 solutions came at about one per 10–20 core-seconds (392 in about 2.3 core-hours).
  - A first core needs about 2^12 of them (per-y stage-2 yield about 2^−12).
  - Stage 2 is a 2+2 join of u-lists of 2^15–2^20 per block, about 2^36 per y when bucketed by a residue of the
    signature (2^25 memory per bucket).
  - Total: **about 2^48 operations, embarrassingly parallel over y.** Roughly days on one core in C, or hours on a
    small cluster. Not run here.
- **So for 573 the route is costed at about 2^48–2^49 [E] with plain 2+2 stage-2 joins, not validated.** That is
  about 2^3–2^4 above the target. It is the most concrete route I found:
  - it is the authors' own 419 design;
  - it is validated end-to-end at 419;
  - it is transferred to the only open order where its slice is predicted non-empty.
- **The gap to ≤ 2^45 is the stage-2 join (about 2^35–2^37 per y).** A level (Z₄ × Z₁₃) decomposition of the u-join,
  of the kind that gave the 359 Q-level count (§7.3), needs to cut it by about 2^4. Alternatively, a richer gadget
  design is needed whose per-y stage-2 yield is nearer the split average of 2^−9. Neither is built here.

## 8. Question 5: the b ≥ 3, p = 11 options and their borders

| option | n′ | \|H′\| | K range (norm identity) | w-screen | border |
|---|---|---|---|---|---|
| 599 c = 6 (b = 5) | 54 | 9 | 17–24 | alive (also with Lemma 6) | **347's border (same b = 5) admits 384 explicit two-valued compressions** [M], below |
| 669 c = 10 (b = 9) | 60 | 6 | 17–25 | alive (also with Lemma 6) | none known |
| 719 c = 16 (b = 15) | 64 | 4 | 17–24 | alive (also with Lemma 6) | none known. Sylvester H(64) ⊇ H₄⊗1₁₆ gives candidate K = H₄ ⊗ (rows 2–16 of H₁₆), corner untested |
| 643 c = 28 (b = 27) | 56 | 2 | 10–14 | alive (also with Lemma 6) | none known |

**Which design transfers.**
- In principle, all of the 359 machinery transfers: border-forced compression, designed gadget, exact stage 1, then
  stage 2. **None of it is cheap at these sizes.**
- **Stage 1 is no longer small.** At 359, the compression pins α to 50.7 bits because |H′| = 4.
  - At 599 c = 6 we have |H′| = 9. With the compression fixed and a balanced gadget, α still carries
    Σ log₂ C(9, #minus) = 136.9 bits, the same for every one of the 384 admissible compressions below
    (`theory/logs/paley_border_c6_599_admissible.json`).
  - So the exact stage-1 join costs about 2^68, not the 2^25 of 359.
- **Stage 2** has 4n′ − K = 192–239 free bits and about n′ = 54–64 equations.
  - Its yield per y is short by about 10 to 20 bits.
  - An exact 2+2 join over u costs about 2^96–2^120 per y.
  - So without a multiplier slice (§7.5) these options are far out of reach.

**Border–compression matching (computable).**
- The corner is forced, C = −K Cb^T Kc^T/(4c), so a compression is admissible for a given border iff K Cb^T Kc^T ∈ ±4c
  entrywise.
- For c = 8 (359's border) the accepted single-deviation compressions are tabulated in §7 (183) and include 359's own.
- **For c = 6 (347's border) at v = 599 [M]** (`theory/paley_border_c6.py`, `theory/logs/paley_border_c6_599.log`):
  - The 347 border reproduces 347's own corner.
  - It annihilates block-constant compressions [P, checked on 20 random shifts]. So admissibility depends only on the
    deviation pattern.
  - Of the 1,244,160 two-valued deviation quadruples (each block's compression takes two values 2 apart) with the
    compressed spectrum required by Corollary 3 (Σ_k |d̂_k(ξ)|² = 4b = 20 for ξ ≠ 1), **384 are accepted by the
    border**.
  - All 384 extend to row sums with Σr_k² = 4v, e.g. (−42, −2, 22, −12).
  - So **599 c = 6 has explicit, border-admissible compressions.** The border is not the obstruction there; the
    stage costs above are.
  - Not tested: compressions with a deviation of ±4 in some block (compatible with Corollary 3, but outside the
    two-valued family).
- For c = 10, 16 and 28 no border is known, and none was constructed here.

## 9. What would change the picture

0. **Finish the 573 slice route (§7.5).** It is the one route here that is validated end-to-end on a known order (419,
   new core in 272 s) and transferred to an open order (573, slice log₂E = +23, stage-1 solutions found).
   - What remains is stage 2: about 2^35–2^37 per y with plain 2+2 joins, against a predicted yield of 2^−12 to
     2^−17 per y for the gadgets found so far. The ranked estimate to a first core is about 2^49.
   - Two concrete steps:
     - decompose the u-join by the Z₄-level characters of Z₅₂ = Z₄ × Z₁₃ (Lemma 6 gives the real ones);
     - widen the gadget design beyond |PAF_w| ≤ 2 and |t_k| ≤ 1, ranked by `theory/paley_coreyield.py`.
   - A 2^4 gain in either brings the route to about 2^45 [E].
1. **A structured stage 2.** A family of ±1 quadruples u on a 2-group G′ with prescribed level-constant spectra and
   Σ_k Im(β̂ conj(û)) = 0 for a given β, i.e. a "complex building set" over Z[π].
   - None exists in the literature I know [L].
   - The exact level decomposition (§2) is the place to look: each level is a complementary-family problem over Q.
2. **Better joint (non-staged) search in the full ansatz**, rather than more structure.
   - The published 359 sits in a space with log₂E = +67.7 and about 300 constraint bits; the authors solved it.
   - **The best-placed open target is 659 at b = 1 (G′ = Z₉₄, p = 7), not 571 or 647.**
     - It has log₂E ≈ +236 to +288 (two estimators, §7.4).
     - It is in the p = 7 regime, where every exact control (n′ = 4–14) behaved as the Gaussian model predicts.
     - Its stage 2 per y is 99–167 bits short, so the search must be joint.
   - 571 and 647 (p = 19; +29 to +43 by either estimator) are in the rigid-gadget regime, where the Gaussian model
     overcounted at every p ≥ 11 control. Lemma 6 and Lemma 4 should be built into any search there.
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
python3 theory/paley_wscreen.py b1 --norho           # the screen without Lemma 6 (same verdicts)
python3 theory/paley_lemma6.py                       # Lemma 6 at every real character: 359 + all b = 1 controls
bash theory/run_paley_b1pipe.sh                      # exact b = 1 pipeline on all small controls (H(116)..H(396), p = 11 empties)
python3 theory/paley_b1pipe.py 99 7 14               # H(396) over the non-cyclic Z14 x Z7 (3 min)
python3 theory/paley_b1pipe.py 111 11 10             # p = 11: exactly empty (Lemma 6 cuts the 8 gadgets)
python3 theory/paley_b1pipe.py verify data/core_v*_b1_paley.txt   # independent check of every control core
python3 theory/paley_s2lattice.py 111 11 10          # stage-2 congruence-lattice test (without Lemma 6)
python3 theory/paley_split_b1.py 659 7 94            # stage split [E]; also 571 19 30, 647 19 34, 573 11 52, 85 7 12, ...
python3 theory/paley_slices.py 659 1 7 94 --slice even    # 659 even slice: -15.3
python3 theory/paley_border_c6.py 599                # 347 border vs the 599 c = 6 compressions (about 25 min)
python3 theory/paley_pipeline.py check479            # c = 4 completion: rebuilds and verifies H(1916) from the 479 core
python3 theory/paley_pipeline.py run179              # 179: no single-bump compression
python3 theory/paley_anatomy419.py                   # published 419 = exact Paley(11) core, x7-invariant (Theorem 1, Lemmas 4, 6)
python3 theory/paley_slices.py 573 1 11 52 --slice mult:29   # Paley + multiplier slices (also 419 ... 38 --slice mult:7)
python3 theory/paley_split_slice.py 573 11 52 29     # stage split inside the slice (also 419 11 38 7)
python3 theory/paley_slicepipe.py given419           # 419 with the published gadget fixed: H(1676) in about 3 min
python3 theory/paley_slicepipe.py guided 419 11 38 7 15 8 400 4   # 419 from scratch (designed gadgets): H(1676) in 272 s
python3 theory/paley_slicepipe.py random 419 11 38 7 16 11 400 1  # random gadgets: no stage-1 solution
python3 theory/paley_slicepipe.py guided1 573 11 52 29 22 15 60 1 13 2700   # 573, stage 1 only (records y's + stage-2 prediction)
python3 theory/paley_s2predict.py theory/logs/paley_slicepipe_stage1_v573_mu29_K22_s1.npz 573 11 52 29   # stage-2 yield per y
python3 theory/paley_b1pipe.py verifycyclic data/core_v167_b1.txt   # bordered-GS completion check: H(668)
```
