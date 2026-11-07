# H(2396) by a Paley plug-in Legendre quadruple, and two combinations of Kumar's methods with ours

Date: 2026-11-07. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_599_plugin.md` (commit d36066f).
Scripts: `theory/h2396_*`. Logs: `logs/h2396/`.

**Tags.** [P] proved; [M] exact computation (script named); [E] estimate. A cap or a timeout is never a zero.

**Credit.**
- The Paley plug-in decomposition, the K-level equations (I), (J), (Q), the defect identity and the index-r plug-ins
  are K. Kumar's (`vatskishan03/hadamard`, Apache-2.0: `search/NOTES.md` §81, §82, §85, §86, §92); so are the
  character-bordered family of §18 and the 599 border block of §103, and the conventions of `verify/verify.py`.
  No code of his was copied; every program here is new.
- The Paley-coset ansatz of `report/paley_design_answer.md` (Theorem 1, Lemma 4, Lemma 6) is this repository's form of
  the same plug-in for q ≡ 3 (mod 4); its Lemma 6 is the ℤ₁₃ level used below.
- The estimator (`theory/open_screen_deficit.py`, `theory/ab_gauss.py`, `theory/q523_gauss.py`) and the exact stage-1
  tool (`theory/open_screen_stage1_exact.py`) are the repository's earlier work, reused unchanged.

---

## 0. Bottom line

<!-- BOTTOM -->
1. **Part A (599 by a Paley plug-in over ℤ₄₆ × F₁₃): killed at gate A2, the go/no-go. No H(2396) was built; A4 was
   not run.**
   - **A1 passes [P, M].** The plug-in system for Legendre quadruples (b = 1, q ≡ 1 mod 4) is derived and checked
     against direct G-level PAF (66/66 random cases, q = 3 … 17, cyclic and non-cyclic K). The defect identity
     14e₊ + 10e₋ = 4(|K| − 1) has three allocations at |K| = 46: (e₊, e₋) = (0, 18), (5, 11), (10, 4).
   - **No analogue has a plug-in core where it was decided [M].** At |K| = 6 and 8 stage 1 is empty (three methods).
     At |K| = 12 stage 1 has exactly 1,327,104 ordered solutions, all "spikes" (each block a constant plus one big
     letter), and every one dies at the ℤ₁₃ level. <!-- N14 --> The Gaussian predicted 2^25 to 2^40 cores there.
   - **The one non-spike stage-1 solution found (|K| = 18, CP-SAT) dies at the coarsest stage-2 level**, the ℤ₁₃
     compression, where the estimator predicted 2^10.2 completions and 2^14.2 S-level solutions for that very y.
     <!-- SAMPLER -->
   - So the plug-in class is **not rich** here. Kumar's undercounting (§81, Szöllősi targets) does not carry over to
     q = 13 Legendre quadruples; they behave like this repository's q = 11 Paley-coset controls (exact zeros
     against +23 to +40 bits).
   - **The cost prong fails too, even under the Gaussian [E].** At |K| = 46 the staged route needs 2^35.3 S-level
     solutions per core: 2^72.2 operations with plain joins, and at least 2^54 even with free S-solutions (each needs
     a 2^19 D-check). Both exceed 2^46 by more than one bit; the calibration says the Gaussian is optimistic.
   - At |K| = 46 the ℤ₁₃ level alone does not decide the slice: 153,392 of the 1,617,624 stage-1 profiles survive it
     [M].
2. **Part B (599 over ℤ₄² × ℤ₃₇, t = 7): killed at B1/B2. The structure has no core at all [M].**
   - The structure is pinned down uniquely: G = GR(4, 2) × ℤ₃₇, Q = GR(4, 2) ≅ ℤ₄², H = ⟨(ξ, ×7)⟩ of order 9 (the
     only multiplier group with 70 orbits), T = {1} ∪ {χ_a : a ∈ ±𝒯*} (unique up to an automorphism commuting
     with H) [P, M].
   - **Its Q-level system has no solution**, by two independent methods (a Fourier-class join and a spatial join on all
     16 lags), for both choices of T. Both methods give the same non-zero count on 3/3 planted targets.
   - So the border question never arises. Kumar's "no border in any symmetric class tested" is true but not the
     obstruction: the first quotient is empty. Gaussian estimate here: 10^29.3 (his 10^24.6) against an exact zero.
   - Side result [P, M]: only 14 ±1 functions on GR(4, 2) have Fourier support in T. So every frame row is one signed
     basis function per block, and the border problem reduces to a ±1 condition on Z M_Tᵀ Z′/4 (§2.4).
3. **Part C (1291, index-4 Szöllősi plug-in on ℤ₁₉ × F₁₇): stopped at C3 on cost.** <!-- C2 -->
   - The repository's estimator, written for this family, reproduces Kumar's numbers: K = ℤ₅ gives 10^5.97 against his
     exact 742,400 = 10^5.87; K = ℤ₁₉ gives 10^9.1 against his 10^8.8 [E].
   - In the staged model the deficit is 26.3 bits per stage-1 solution, and ℤ₁₉ has no subgroup for an S/D split.
     So the floor is 2^26.3 stage-1 solutions per core, each needing a join over per-block fibres of 2^34.8:
     **2^61 per core**, against the 2^44 bar.
   - "The obstruction is search" is therefore confirmed in our cost model. It is our kind of slice, but not a cheap
     one. 2083 was not run (1291 was not decided early).
4. **Set aside** (no compute; one line each in §4): 1187, 1987, 1223, 719 at index 16, 1879 and 2215.
<!-- /BOTTOM -->

---

## 1. Part A: 599 by a Paley plug-in Legendre quadruple

### 1.1 A1: the plug-in system [P, M]

**Set-up.**
- G = K × F_q, q = 13, H = {1} × S (S the squares of F_q*). A block constant on H-orbits is three functions on K:
  x₀ on K × {0}, x₁ on K × S, x₂ on K × N.
- With e = (x₁ + x₂)/2, t = (x₁ − x₂)/2 (exactly one of them is non-zero at each k) and s = x₀ − e, the block is the
  group-ring element X = s ⊗ δ₀ + t ⊗ χ + e ⊗ 1_F, with χ the Legendre symbol.

**The K-level system [P].** In ℤ[F_q]: δ₀·f = f, 1·1 = q·1, 1·χ = 0, χ⁽⁻¹⁾ = εχ (ε = χ(−1)) and
χ·χ⁽⁻¹⁾ = qδ₀ − 1. Hence, with ~ the inversion on K,

    X X⁽⁻¹⁾ = (s s~ + q t t~) ⊗ δ₀ + (t s~ + ε s t~) ⊗ χ + (s e~ + e s~ + q e e~ − t t~) ⊗ 1.

δ₀, χ and 1 are independent, so the Legendre-quadruple target (4m + 4)δ_G − 4·1_G (m = q|K|) is equivalent to

- (I) Σᵢ (sᵢ sᵢ~ + q tᵢ tᵢ~) = (4m + 4) δ_K;
- (Q) Σᵢ (tᵢ sᵢ~ + ε sᵢ tᵢ~) = 0, with ε = +1 for q = 13;
- (J) Σᵢ (sᵢ eᵢ~ + eᵢ sᵢ~ + q eᵢ eᵢ~ − tᵢ tᵢ~) = −4·1_K.

These are Kumar's (I), (J), (Q) of §81 with the b = 1 right-hand sides.
- Stage 1 is (I) + q(J): u = s + qe = x₀ + (q − 1)e, and Σᵢ PAF_{uᵢ} = (4m + 4)δ − 4q on K. The letters are ±1
  (clean), ±13 (x₀ = e, type e₊) and ±11 (x₀ = −e, type e₋).
- Given u, stage 2 is: Σᵢ PAF_{tᵢ}(h) = −(1/13) Σᵢ PAF_{sᵢ}(h) for h ≠ 0, together with (Q), which is linear in t.
- **Check [M]** (`theory/h2396_plugin_system.py check`, `logs/h2396/a1_system_check.json`): the G-level PAF equals
  the K-level formula on 66 of 66 random inputs, for q ∈ {3, 5, 7, 11, 13, 17} and K ∈ {ℤ₂, ℤ₃, ℤ₄, ℤ₅, ℤ₆, ℤ₂²,
  ℤ₃ × ℤ₂, ℤ₂ × ℤ₃}, together with the stage-1 compression.

**The defect identity [P].** (I) at 0 gives (q + 1)e₊ + (q − 3)e₋ = 4(|K| − 1), i.e. 7e₊ + 5e₋ = 2(|K| − 1) at
q = 13 (`h2396_plugin_system.py identity`, `logs/h2396/a1_defect_identity_q13.jsonl`).
- At |K| = 46: (0, 18), (5, 11), (10, 4), i.e. 18, 16 or 14 defects. **A1's kill condition is not met.**
- No allocation exists at |K| = 2, 3, 4, 5, 9 or 10. The first even analogues are |K| = 6 (0, 2), 8 (2, 0),
  12 (1, 3), 14 (3, 1), 16 (0, 6), 18 (2, 4), 20 (4, 2), 22 (1, 7) and (6, 0).

**Row sums and the ℤ₁₃ level [P].**
- r_i = Σ_K u_i = ŝ_i + 13ê_i (hats are sums over K), and Σ r_i² = 4. So for even |K| the rows are (±2, 0, 0, 0).
- The trivial character of K applied to (I), (Q), (J) gives the ℤ₁₃ level: Σ(ŝ_i² + 13t̂_i²) = 4m + 4 and
  Σ ŝ_i t̂_i = 0. Equivalently, with block 0 the block of row sum ±2:

      −2 r₀ ê₀ + 13 Σ ê_i² + Σ t̂_i² = 4|K|,        r₀ t̂₀ = 13 Σ ê_i t̂_i,

  with t̂_i ≡ ê_i ≡ D_i (mod 2) (D_i = the defects of block i), |t̂_i| ≤ |K| − D_i, and
  |r_i − 12ê_i| = |Σ x₀| ≤ |K|.
- It is the compression along K to F₁₃, so it is the energy identity of `report/h4892_answer.md` §5:
  Σ[w(0)² + 6(w_S² + w_N²)] = 4 + 48|K| is (I) + (J) at the trivial character.
- **Consequences:**
  - 13 | t̂₀;
  - |ê_i| ≤ (|K| + 2)/12, i.e. |ê_i| ≤ 1 for |K| ≤ 21 and |ê_i| ≤ 4 at |K| = 46;
  - **parity lemma [P]:** for even |K| no core has an odd number of defects in all four blocks. Mod 8 the three left-hand
    terms are each ≡ 4 (−2r₀ê₀ with ê₀ odd; 13 times four odd squares; four odd squares), so the left side is ≡ 4,
    while 4|K| ≡ 0. This is
    the q = 13 form of Corollary 6.1 of `report/paley_design_answer.md`.
- The order-2 character ψ of K gives the same two equations with ψ-sums (Σ_i A_i² = 4m + 4 is the ℤ₂ screen of
  `theory/gp_z2screen.py`).

**Congruences on stage 1 [P].** With u = x₀ + 12e, stage 1 at d ≠ 0 reads a(d) + 12c(d) + 144w(d) = −52. Here
a = Σ PAF_{x₀}, c = Σ [CCF(x₀, e) + CCF(e, x₀)] and w = Σ PAF_e. So:
- **Σᵢ PAF_{x₀,i}(d) ≡ 8 (mod 12) for every d ≠ 0** (it is ≡ 0 mod 4 automatically for even |K|, so this is a mod-3
  condition on the x₀-quadruple: the analogue of Kumar's mod-3 level-1 filter at 643);
- the w-magnitude bound of Lemma 4 of `report/paley_design_answer.md`: at |K| = 46, w(d) ∈ [−4, 3];
- Σᵢ PAF_{sᵢ}(h) ≡ 0 (mod 13) at every h ≠ 0, so the stage-2 target is integral; this is implied by stage 1;
- the per-block PSD bound |ûᵢ(φ)|² ≤ 4m + 4 at every nontrivial character, and |ŝᵢ(φ)|² ≤ 4m + 4 for stage 2.

**Controls [M].**
- Stage-1 counts: an independent brute force (`theory/h2396_s1brute.py`: every letter sequence per block, a
  dictionary join) equals the exact enumerator (`theory/h2396_s1enum.c`): q = 7, |K| = 4 gives 12,288 for (0, 3) and
  6,144 for (1, 1); q = 13, |K| = 6 gives 0.
- Full system: a CP-SAT model of the whole plug-in system (`theory/h2396_fullcp.py`) finds cores at q = 7, |K| = 4
  (the shape of the repository's H(116)) and q = 5, |K| = 4 and 6, each checked by direct G-level PAF; it proves
  q = 13, |K| = 6 and 8 infeasible (`logs/h2396/fullcp_controls.jsonl`).
- **End to end [M]:** a G-level brute force sharing no code with the K-level tools (`theory/h2396_gbrute.py`: every
  plug-in block, direct PAF on ℤ₄ × F₅, dictionary join) counts 1,572,864 ordered cores at q = 5, |K| = 4.
  - The pipeline (`h2396_s1enum` stage 1 for every allocation, then `h2396_stage2` on every class) gives the same
    1,572,864, with or without the ℤ_q-level prefilter (`logs/h2396/gbrute_n4_q5.log`).
  - With the prefilter the stage-1 set shrinks from 557,056 to 282,624 solutions, but the core count is unchanged.

### 1.2 A2: the existence trend [M, E]

**Analogues.** K = ℤ_n, q = 13, slice `["b1", 13n + 1, [], H1, 13, 2, []]`.
- Estimates: the repository estimator (`theory/h2396_gauss_analogues.py`, `logs/h2396/gauss_analogues.jsonl`).
- Exact stage 1: three independent methods.
  - `theory/h2396_s1enum.c`: canonical defect patterns per block under translation, reflection and sign, then an
    exact 2+2 join over x₀ with the row-sum restriction; orbit-weighted ordered counts.
  - The repository tool `theory/open_screen_stage1_exact.py` (methods A and B).
  - CP-SAT (`theory/h2396_s1cp.py`).
- Exact stage 2 for every stage-1 class: `theory/h2396_stage2.c` (§1.3).
- Table: `theory/h2396_a2_table.py --md`.

<!-- A2TABLE -->

**Notes.**
- **|K| = 6, 8 [M].** Stage 1 is empty: repository tool (A and B, planted 2/2 each), `h2396_s1enum` and CP-SAT
  (INFEASIBLE in under a second).
- **|K| = 12 [M].** Stage 1 has 1,327,104 ordered solutions in 4,368 classes up to block order and sign; the
  repository tool's method B finds the same 4,368. Every one is a spike quadruple (each block a constant ±1 sequence with one big letter: one ±13 block, three
  ±11 blocks), so all four blocks have one defect. The parity lemma kills them, and the exact stage-2
  counter gives 0 for every class (`logs/h2396/s1enum/stage2_n12.log`). **No core.**
- **The spike family [P].** If every block is c_i·1 + β_i δ_{p_i} (one defect each), stage 1 holds iff
  Σ c_i b_i = −2|K| − 22, with b_i the big letters. Since |c_i b_i| ≤ 13, this needs |K| ≤ 15. The family fills
  stage 1 at |K| = 12 and is all CP-SAT finds in seconds at |K| = 12 and 14; it is always dead by the parity lemma.
- <!-- N14NOTE -->
- **|K| = 18.**
  - CP-SAT found one non-spike stage-1 solution after 40 minutes (`logs/h2396/s1cp/n18_y0.txt`), with defects
    D = (1, 2, 2, 1).
  - The exact counter gives 0 completions and 0 S-level solutions. The quotient ladder (`theory/h2396_s2ladder.py`)
    shows it is already dead at ℤ₁, the ℤ₁₃ level.
  - By hand: ŝ = r − 13ê = (−11, 0, 0, 13), so Σ ŝ_i t̂_i = 0 forces (t̂₀, t̂₃) = j(13, 11). j = 0 violates the
    parity of t̂₀, and j = ±1 violates Σ t̂² = 50.
  - The estimator, read at this real y, predicted E per y = +10.18 bits and 14.21 bits of S-level solutions
    (`logs/h2396/realy_n18.json`).
  - At |K| = 14 the estimator's lattice test does see the spike's parity obstruction (E per y = −∞,
    `logs/h2396/realy_n14.json`). At 18 the obstruction is an inequality, which no Gaussian sees.
- **The ℤ₁₃ level over all stage-1 profiles** (`theory/h2396_zq_profiles.py`). A profile is (D_i, ê_i, r_i); the
  table counts alive / all:

| \|K\| | 12 | 14 | 16 | 18 | 46 |
|---|---|---|---|---|---|
| alive / all | 464 / 720 | 80 / 976 | 1,120 / 2,592 | 160 / 2,592 | 153,392 / 1,617,624 |
| odd-D blocks in alive profiles | 0, 2 | 0 | 0, 2 | 0 | 0, 2 |

  At |K| = 14 and 18 a core needs an even number of defects in every block. At |K| = 46 the level leaves 9.5% of the
  profiles and does not decide the slice (`logs/h2396/zq_profiles_n46.json`).
- <!-- SAMPLERNOTE -->

**The trend and its extrapolation.**
- No non-zero count was found at any analogue: exact zeros at |K| = 6, 8, 12 <!-- N14SHORT -->, and zero completions
  for every sampled instance. So log(count) cannot be fitted.
- What can be said honestly:
  - the estimator over-predicts by at least 25 to 40 bits where the count is exact;
  - it over-predicts by at least 10 bits per y at the one real y beyond the spike regime;
  - in every case the instances die at the coarsest stage-2 level, as at 4892 and 7572.
- Kumar's concentration near the Paley point does not appear: the defects are forced (14 to 18 at |K| = 46), and
  the small analogues are empty or rigid.
- This agrees with Kumar's own real plug-in censuses.
  - §86: "the estimator gave up to 10^6.9 solutions on exactly empty small shapes".
  - §92 (Spence shapes of odd |G|): "it gives up to 10⁷ on exactly empty small shapes".
  - The undercount of §81 was a property of the Szöllősi (complex) classes. For Legendre quadruples (Spence's
    array) both his data and ours show over-prediction.
  - His §92 also kills the other U = 1 Paley slice of 599 (ℤ₂ × ℤ₁₃ × F₂₃), as `report/paley_design_answer.md`
    Theorem 5 does.
- **Kill on the first prong: no solution appears at any |K| decided, well below 46.**

### 1.3 A3: the sampler [M]

**Level 1.** CP-SAT on stage 1 (`theory/h2396_s1cp.py`).
- With `--zq`, the ℤ₁₃ and ψ levels of stage 2 are added: integers t̂_i and t′_i with the exact parity and range
  links, including t̂ ± t′ = 2·(sum of t over the even / odd clean positions). Only instances that survive both
  real-character levels are produced.
- Controls: `--zq` finds q = 7, |K| = 8 at once, and that instance has 224 completions. At q = 13, |K| = 12, where
  every stage-1 solution is known to die, it never returns FEASIBLE (UNKNOWN at 600 s).
- A stochastic local search (`theory/h2396_s1sa.c`) solves q = 7 at once, but fails on planted q = 13 targets for
  |K| ≥ 22. Its failures on the real target are therefore not evidence of anything, and it is not used.

**Level 2.** The exact stage-2 counter `theory/h2396_stage2.c`.
- It splits along L = {0, |K|/2}: an S-level 2+2 join on ℤ_{|K|/2} (PSD-filtered lists), then for every S-solution an
  exact D-check (a 2+2 join on the full K-level key). Every completion is re-checked from scratch.
- Validation:
  - equal to an independent reference (`theory/h2396_stage2_ref.py`, direct enumeration of t, no split) on all real
    and planted instances it can run (q = 7: 422 and 320; q = 13: 0; planted |K| = 8, 10: 44 and 6);
  - planted recovery 7/7 at |K| = 8 to 20 (`logs/h2396/stage2_planted.log`).
- The stage-1 enumerator is validated by the brute force and the repository tool (§1.1, |K| = 12).

**Time per instance.**
- At |K| = 18: level 1 about 40 CPU-minutes (CP-SAT), level 2 14 s.
- At |K| = 46: level 2 is out of reach by this design. A block has about 19 clean-clean pairs and 4 clean-defect
  pairs, so 3¹⁹·2⁴ ≈ 1.9·10¹⁰ ≈ 2^34 S-configurations, 1–2 TB per list (keys of 24 lag values) before the join.
- The deficit model's own per-instance cost at |K| = 46 is ≥ 2^47 (2^28 S-solutions, each a 2^19 D-check).
- No level-1 instance was produced at |K| ≥ 22 within the budget (<!-- N22 -->). A cap is not a zero.

### 1.4 The cost at |K| = 46 [E]

From `logs/open_screen/deficit_s*.jsonl` (synthetic y; Gaussian, so optimistic by the calibration above):

| split | C_S | H_S | h_D | D_def | S-step | cost (log₂ per core) |
|---|---|---|---|---|---|---|
| L = ℤ₂ | 101.9 | 130.0 | 38.0 | 35.3 | 36.9 | **72.2** |
| L = ℤ₂₃ | 22.6 | 25.9 | 142.1 | 10.5 | 9.6 | 81.6 |
| L = G′ | 12.4 | 14.9 | 153.1 | 9.7 | 5.0 | 86.2 |

- E per y = −7.2 bits (H₂ = 168, C₂ = 175.2).
- Even if S-solutions were free, a core needs 2^35.3 D-checks of at least 2^19 each, so ≥ 2^54.
- The 523-style two-level Wagner tree needs a subgroup inside ℤ₂₃, and there is none.
- The Kumar-style per-instance join is worse: per-block fibres of 2^42, against 2^20 at 643.
- **Kill on the second prong as well:** ≥ 2^54 against 2^46 + 1 bit.

### 1.5 Verdict

**Part A stops at A2.** It is killed on both prongs: no solutions at any decided analogue, and a Gaussian cost at
|K| = 46 that is 8 to 26 bits over budget. A4 (production) was not run, and no matrix is claimed.

---

## 2. Part B: 599 over ℤ₄² × ℤ₃₇ (Kumar §103)

### 2.1 B1: the structure, reconstructed uniquely [P, M]

**The structure.**
- Q = GR(4, 2) = ℤ₄[ξ]/(ξ² + ξ + 1), additively ℤ₄². G = Q × ℤ₃₇ with |G| = 592 = 599 − 7, K = 0 × ℤ₃₇, t = 7,
  border 4t = 28.
- Characters: χ_a(x) = i^{Tr(ax)}, with Tr(x₀ + x₁ξ) = 2x₀ − x₁ (mod 4).
- Target: Σᵢ PAFᵢ(g) = −Γ(g + K) for g ≠ 0, with Γ = 4 Σ_{χ∈T} χ.

**H is forced** (`theory/h2396_b_alternatives.py`).
- Among all ⟨(A, ×m)⟩ ≤ GL(2, ℤ₄) × ℤ₃₇*, exactly those with ord A = 3 and ord m = 9 have 70 orbits
  (48 pairs (A, m), all conjugate to (ξ, ×7)). The orbits are 1 + 5 + 4 + 60.
- The only non-cyclic subgroup of order 9 has 78 orbits.
- So H = ⟨(ξ, ×7)⟩ (`theory/h2396_b_qlevel.py orbits`).

**T is forced up to an automorphism.**
- T must be H-invariant and closed under conjugation, so its six order-4 characters form one ⟨ξ, −1⟩-orbit. There
  are two such orbits: a ∈ ±𝒯* = ±{1, ξ, ξ²} and a ∈ ±(1 + 2ξ)𝒯*.
- The automorphism x ↦ (1 + 2ξ)x commutes with ξ and swaps them. Both were run anyway.

**The Q-level [P].**
- y_i(x) = Σ_k x_i(x, k) is ξ-invariant. y(0) = ε₀ + 9·(4 signs) and y(x ≠ 0) = ε₀ + 3·(12 signs), from the H-orbits
  on the fibres.
- Σᵢ |ŷᵢ(χ)|² = 28 on T and 2396 off T. The 16 characters fall into four classes (trivial; T \ 1; order 2; the other
  units), so the Q-level is four equations: Σ r² = 28 (r = (±4, ±2, ±2, ±2)), Σ |ŷ(χ_1)|² = 28,
  Σ ŷ(χ_2)² = 2396 and Σ |ŷ(χ_{1+2ξ})|² = 2396.

### 2.2 B2: the Q-level classes: there are none [M]

- **Fourier-class join** (`theory/h2396_b_qlevel.py classes`, `logs/h2396/partB/qlevel_classes.json`): 428 admissible
  blocks (energy, rows, all four class bounds), 66 key classes, **0 ordered quadruples**. The other T gives
  428 blocks and 0 as well.
- **Spatial join** (`theory/h2396_b_qcheck.py`, shares no code with the above): 31,404 blocks with 3,109 distinct
  PAF vectors over all 16 lags of Q, **0 ordered quadruples**.
- **Planted controls** (`theory/h2396_b_qplant.py`, `logs/h2396/partB/qlevel_planted_rerun.log`): random H-invariant
  quadruples give Q-level targets, and the two methods count the same non-zero numbers on 3/3 (799,426,560;
  1,046,019,840; 1,759,427,328).
- **So the structure has no core [M].** There is no compression class to border, and B3 is not reached.

### 2.3 Estimate [E]

- The lattice Gaussian of the full structure (`theory/h2396_b_gauss.py`: 35 lag classes, 4·70 free bits) gives
  2^97.3 = 10^29.3. Kumar's figure is 10^24.6, presumably with another normalisation.
- Either way: positive estimate, exact zero, at the coarsest quotient.

### 2.4 Side result: the border problem in this family [P, M]

- Only **14** ±1 functions on GR(4, 2) have Fourier support in T. They are ±1 and ±(six functions), each supported on
  one conjugate pair {χ_a, χ_{−a}}, forming an orthogonal basis of V_T per block up to sign.
- So a frame row is one signed basis function per block. The frame is R̂ = Z·B, with Z a weighing matrix W(28, 4)
  carrying one ±1 in each block's group of 7 columns.
- With bt_border's reduction, a border exists iff X = Z M_Tᵀ Z′/4 is a ±1 matrix for two such Z, Z′. Here M_T is Mc
  restricted to V_T: block-diagonal by frequency, depending only on the row sums and on ŷ_i(χ_1).
- This is the setting of Kumar's §102 counting argument against single characters per block. It is recorded for the
  other ℤ₄² members of the family; no search was needed here.
- **A normal form for B2's question [P].** M_T depends only on (r_i) and (ŷ_i(χ₁)), because ξ-invariance makes the
  three frequency pairs carry the same values.
  - The data are finite: r is a signed permutation of (4, 2, 2, 2) (from Σ r² = 28 with r even), and
    ŷ_i(χ₁) = 2w_i with w_i ∈ ℤ[i] and Σ|w_i|² = 7 (y is odd-valued).
  - So the border problem for every core of every ξ-invariant structure with this T reduces to a finite list of
    (r, w) classes. One certificate per class would settle it, as Theorem B does for k4.

    That list was not worked out here, because this block has no Q-level class to border.

### 2.5 Verdict

**Part B is closed at its first gate: the block is real, but it is not a border block.** The structure is empty at
Q = ℤ₄² before any border question. This agrees with Kumar's rule in §103 ("test the border and the level-1 join
over Q before any ladder"): the level-1 join over Q is what fails.

---

## 3. Part C: the index-4 Szöllősi plug-in for 1291

### 3.1 C1: the formulation [P]

**Core.** A zero-diagonal Szöllősi array on G = ℤ₁₉ × F₁₇, N = 323: A and B symmetric, vanishing exactly at 0; C and
D ±1; Σ PAF = 1290δ. These are Kumar's `verify.py` conventions.

**Plug-in.** H = {1} × S₄ with S₄ = {1, 4, 13, 16}, the fourth powers; −1 ∈ S₄. A block is x₀ on K × {0} and
x_j on K × C_j, j = 0..3 (C_j = 3^j S₄).

**Stage 1** (compression along F₁₇).
- y = x₀ + 4σ with σ = Σ_j x_j, letters {±1, ±7, ±9, ±15, ±17}; at A(0) and B(0), {0, ±8, ±16}.
- Σᵢ PAF_{yᵢ} = 1290δ on ℤ₁₉. Its lag 0 is Kumar's balance identity Σ σ(4σ + 2x₀) = 304.

**ℤ₁₇ level** (compression along ℤ₁₉).
- w = aδ₀ + Σ_j W_j 1_{C_j}, with a = Σx₀ and W_j = Σ_k x_j(k).
- Σᵢ [aᵢ² + 4Σ_j W_ij²] = 1290, and for l = 0..3, Σᵢ [2aᵢW_il + Σ_{j,j′} W_ij W_ij′ (j − l, j′ − l)] = 0, with the
  cyclotomic numbers of order 4 mod 17 (`theory/h2396_c_zq.py`).

**Generalised-Paley slice.** Stage 1 is the K-level with fixed per-k class sums σ(k). The ℤ_q level is the coarse
F-side level. Stage 2 is the cyclotomic remainder (49 nonzero G-lag classes in all).

### 3.2 C2: stage 1 and the ℤ₁₇ level

**Exact enumeration is out of reach.** A C-block has about 10^13 letter sequences inside the PSD and energy bounds.

**Instead, CP-SAT on stage 1 and the ℤ₁₇ level jointly** (`theory/h2396_c_s1cp.py`). Each (k, class) value is a
variable, so both levels describe the same blocks. Every answer is re-checked by direct PAF on ℤ₁₉ and on F₁₇.
- Controls: FEASIBLE at K = ℤ₃ and ℤ₅ (q = 17) in 1 s and 20 s.
- <!-- C2RESULT -->

### 3.3 C3: cost [E]

The estimator (`theory/h2396_c_estimate.py`, `logs/h2396/partC/estimate_z19.json`) follows the model of
`report/open_screen_answer.md` §2.3:
- stage 1: an exact y-level energy DP, then the lattice Gaussian of lags 1..9;
- stage 2: H₂ is the exact fibre (1, 4 or 6 class splittings per position), and C₂ is the lattice Gaussian of all
  G-lag classes at the sample mean.

| K | stage 1 log₂E | H₂ | C₂ | E per y | total | Kumar |
|---|---|---|---|---|---|---|
| ℤ₅ (control) | 16.4 | 38.4 | 35.0 | +3.4 | 10^5.97 | 742,400 = 10^5.87 (exact) |
| ℤ₁₉ (1291) | 56.6 | 139.1 | 165.4 | −26.3 | 10^9.1 | 10^8.8 (tilted Gaussian) |

- A staged search needs 2^26.3 stage-1 solutions per core.
- Each needs an exact stage-2 solve over per-block fibres of 2^(H₂/4) = 2^34.8, and ℤ₁₉ (prime) offers no S/D split.
- **Floor: 2^61 per core, against 2^44. Stop.**
- Kumar's other figures (95 bits per block, joins only to |K| ≤ 5) say the same from the other side.

### 3.4 2083

Not run: 1291 was not decided early.
- 2083 also needs a different model: q = 13 has −1 ∉ S₄, so A and B pair the classes C_j ↔ C_{j+2}; and the
  bordered form is needed.

---

## 4. Set aside after review (no compute)

- **1187** (k4 over ℤ₄² × ℤ₂ × ℤ₃₇): Kumar's join work 10^26 ≈ 2^86, far above 2^46.
- **1987** (ℤ₈² × ℤ₃₁): estimate 10^1.4 with work 10^21 ≈ 2^70.
- **1223** over GR(4, 3) × ℤ₁₉: heavy-leaf fibres of 10^7.2–7.6 per block and a yield of 10^−4.1 per leaf, about 2^63 per
  solution by his figures. The S/D splits measured here (2^61 floor for 1291, §3.3) do not bridge it.
- **719** at index 16: work 10^21.8 ≈ 2^72, and exact quotient tests killed every index-16 structure he took to the
  battery.
- **1879 and 2215** (t = 1 over ℤ₇ × ℤ₆₇, ℤ₇ × ℤ₇₉): about 0.03 inequivalent solutions expected.

---

## 5. Next computation

<!-- NEXT -->
1. **Gate every staged slice at its real-character levels while stage 1 is generated, not after.**
   - At 599 (this report), 4892 and 7572, every stage-1 solution met died at the ℤ_q level.
   - `theory/h2396_zq_profiles.py` decides that level over all stage-1 profiles in minutes. `h2396_s1cp.py --zq`
     builds it into level 1.
   - For a slice with |ê| unrestricted (larger |K|) the profile table is the cheap first gate. It does not decide 599
     (9.5% of profiles survive).
2. **599 b = 1 plug-in at |K| = 46:** only if a level-1 method produces ℤ₁₃-surviving stage-1 solutions at |K| ≥ 22.
   Then read the estimator at those real y's. None was produced here.
3. **The ℤ₄² × ℤ_p border family (Part B):** before any border search, run the four-equation Q-level join
   (`theory/h2396_b_qlevel.py`) on the other members (669 ℤ₄₁, 719 and 739 ℤ₄³ × ℤ₁₁). It takes seconds.
<!-- /NEXT -->

---

## 6. Reproduce

```
python3 theory/h2396_plugin_system.py check
python3 theory/h2396_plugin_system.py identity 13
gcc -O3 -march=native -o h2396_s1enum theory/h2396_s1enum.c -lm && ./h2396_s1enum 12 13 1 3 sols_n12.txt
gcc -O3 -march=native -o h2396_stage2 theory/h2396_stage2.c -lm && theory/h2396_stage2_all.sh 12 13 sols_n12.txt ./h2396_stage2
python3 theory/h2396_s1brute.py 4 7
python3 theory/h2396_s1cp.py 18 13 2 4 --secs 3600 --workers 1 [--zq]
python3 theory/h2396_s2ladder.py 18 logs/h2396/s1cp/n18_y0.txt 1,2,3,6
python3 theory/h2396_zq_profiles.py 46
python3 theory/h2396_realy_estimate.py 18 logs/h2396/s1cp/n18_y0.txt
python3 theory/h2396_a2_table.py --md
python3 theory/h2396_b_qlevel.py orbits && python3 theory/h2396_b_qlevel.py classes
python3 theory/h2396_b_qcheck.py && python3 theory/h2396_b_qplant.py 3 7
python3 theory/h2396_b_alternatives.py && python3 theory/h2396_b_gauss.py
python3 theory/h2396_c_s1cp.py 19 17 --secs 8400 --workers 1
python3 theory/h2396_c_estimate.py 19 17 2000 2
```
