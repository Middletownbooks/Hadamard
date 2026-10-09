# Open orders 2000–3000: review of unattempted ideas

Date: 2026-10-06. Sources: every report, pre-registration and README in the repo, plus the screen and enumeration
logs (`logs/b1_mult_open_*`, `logs/b1_rerank.log`, `logs/bge3_screen_*`, `logs/kt523_*`), and today's T-matrix work.

**Open list:** v = 523, 537, 571, 573, 599, 643, 647, 659, 669, 719, 739.
- 515 (order 2060) was done by Schneider.
- 721 (order 2884) was done this session (`cw/RESULT_2884.md`).

**Context:** both of those are 2026 results from other groups, Alpöge et al. below 2000 and Schneider + Codex at 2060.
This range is being worked on actively.

**Reference points for difficulty** (constraint bits n; lower is easier):
- T(103): solved blind by the support-split k-tree in 12 min at n ≈ 63 [M].
- v = 487, b = 0, d = 9: solved by the flat-pool k-tree at n′ ≈ 83 (flat-restricted) [M].
- v = 523, b = 0, d = 9: n = 142 (n′ ≈ 91 flat). About 0.3 expected hits accumulated, none found [M/E].

## 1. Best known route per open order (cheapest feasible configuration on record)

| v (order) | b = 0 | b = 1 | b ≥ 3 | T-matrix / Cooper–Wallis | best | status |
|---|---|---|---|---|---|---|
| 523 (2092) | **d = 9: 142 (flat ≈ 91)** | d = 14: 172 | c = 4: 226 | = b = 0 | 142 | **partially run** (local + Kaggle run 1, about 0.3 expected hits) |
| 537 (2148) | none | ruled out (d = 22, 11 infeasible); d = 6: 290 | no border (v+1 = 2·269) | T(179): no subgroup | 290 | effectively closed |
| 571 (2284) | d = 5: 281; **d = 15 twisted: see §2** | d = 12: 203 | c = 22: 192 (screen undecided) | = b = 0 | 192–203, or about 95 twisted [E] | never run |
| 573 (2292) | none | **d = 15: 193** | c = 14: 269 | T(191) dead | 193 | never run |
| 599 (2396) | none | d = 11: 267 | c = 20 (Z20): 200 (survives, undecided) | dead | 200 | never run |
| 643 (2572) | d = 3: 530 | none | c = 14: 232 (screen empty); c = 4: 425 | = b = 0 | 232 | never run |
| 647 (2588) | none | d = 12: 213 | c = 18: 235 | dead | 213 | never run |
| 659 (2636) | d = 7: 237 | d = 23 infeasible | c = 6: test 1 failed; c = 10: 237–242 | = b = 0 | 237 | never run |
| 669 (2676) | none | none | – | T(223, k = 3): 157 | 157 | twisted T screen: no gain (`tmat_routes.md` §8) |
| 719 (2876) | none | d = 179 infeasible | **c = 20: 227** (mod-5 check undecided) | dead | 227 | never run |
| 739 (2956) | d = 9: 210 | d = 60: −203 log2E | c = 10: 229–235, not screened | = b = 0 | 210 | never run |

## 2. Unattempted ideas, ranked

### 1. Twisted multipliers for prime-v b = 0 cores: new, and opens 571 at d = 15

**Ansatz.**
- The multiplier group U has order d, with 3 | d and d odd.
- The generator u fixes block A and cycles the blocks B → C → D → B, allowing signs whose product is +1.
- Equivalently, X(u·x) = g(X)(x) with g a 3-cycle in the signed block-permutation group.
- ΣPAF is invariant under this action, so there is one equation per class of ⟨U, −1⟩, i.e. (v − 1)/(2d) classes.

**Why it matters.**
- With no twist, d = 15 at v = 571 has 39 orbits but no realisable row sums (`nbits.py` gives ∞).
- With the twist, block A's sum is ≡ ±1 (mod 30), and B, C, D share one sum ≡ ±1 (mod 10), so the row-sum
  condition becomes a² + 3s² = 2284.
- **a = 31, s = 21 works** [M, arithmetic].
- That gives a 39-orbit, 19-class slice. The row-sum-blind estimate for this shape was n ≈ 95 bits [E, needs a proper
  twisted estimate]. That is far below every other route for 571 (192–281), close to the demonstrated k-tree range
  (63–91), and cheaper than 523's 142.

**Other twisted slices** (exact arithmetic):
- 571/57, 523/87 and 739/123: infeasible.
- 523/9 and 739/9: feasible, but with the same orbit counts as untwisted. 523/9 forces profile (35, 17, 17, 17),
  which Kaggle run 1 already searched untwisted.
- 599, 647, 659 and 719 have no d with 3 | d among the useful divisors.

**Precedent:** the published v = 479 core has a twisted automorphism of order 12 (`report/bge3_construction.md` §3).
`REPORT.md` §6 lists "twisted automorphisms … not searched".

**Next steps:**
1. Write a twisted-GS estimator, adapted from `tmat_twist.py` with 16 sign labels per orbit and no product condition.
   Validate it against brute force on small primes p ≡ 1 (mod 3).
2. Decide feasibility at 571/15.
3. If n is at most about 100: run a blind positive control with the same twist at a smaller prime, then search 571/15.
   The flat-pool k-tree applies (four lists A, B, C, D), with C and D determined by B on each orbit, which needs a
   2-list or 3-list variant. Alternatively, the support-split idea from `tmat_kt.py`.

### 2. Continue v = 523, b = 0, d = 9 on Kaggle (partially run)

- The v2 notebooks with the INSTANCE fix are ready (`kaggle/kt523_A2_*.py`, `kt523_B2_*.py`).
- Accumulated so far: about 0.3 expected hits.
- At about 3·10⁻³ hits per bucket, one expected hit takes about 300 buckets, roughly 50–75 Kaggle session-hours [E].
- This is the most mature validated route in the repo, but it competes with Zarankiewicz/Hill for Kaggle time.

### 3. b = 1 searches for 523, 573, 571 and 647 (never run)

| v | d | n |
|---|---|---|
| 523 | 14 | 172 |
| 573 | 15 | 193 |
| 571 | 12 | 203 |
| 647 | 12 | 213 |

- The tools exist and were validated: `kt_b1.py` passed a planted v = 419 test.
- These are 80–120 bits beyond anything demonstrated. Low priority.

### 4. Undecided b ≥ 3 screens (never decided)

| v | configuration | n |
|---|---|---|
| 571 | c = 22 | 192 (screen returned no verdict) |
| 599 | c = 20, Q = Z20 | 200 (survives) |
| 719 | c = 20 | 227 (mod-5 check undecided) |
| 643 | c = 14 | 232 |
| 739 | c = 10 | not screened |

The c = 6 pipeline failed its blind positive control (`posctrl/RESULTS.md`), so even a pass would feed an
unvalidated lift. Low priority.

### 5. Literature routes considered and not recommended

- **Williamson(-type) quadruples of order 191, 179 or 223**, combined with T(3) via Cooper–Wallis, would give 573,
  537 or 669.
  - W-type of order 179 or 223 would also have given 716 or 892, which were open until 2026.
  - Circulant symmetric quadruples at about 190 sit near the Gaussian non-existence edge. W(65) is already empty.
- **Yang / base-sequence multiplication.** It needs BS(n+1, n) with n ≥ 89 (537 = 3·179, 573 = 3·191,
  669 = 3·223). Known BS only reach n ≈ 40, and no open v factors as (Yang number)·(2n+1) with a known BS.
- **D-optimal transfer** (`report/dopt_transfer_audit.md`): a different problem. Only the arithmetic overlaps
  (537, 573 and 643 satisfy the two-squares condition).
- **The b = 1 half-shift structure** (`logs/b1_halfshift_cost.log`): about 2000 bits at 515 and 719. Hopeless.
  Update (2026-10-09): the structure is now characterised as a two-arc template (`report/alpoge_methods_forensics.md` §4). With the published-style arc supports, each of the two independent halves has log₂E ≈ +31 at 167/179, +6 at 419, +2 at 443 and −13 at 523 (about 535 constraint bits per half). The verdict stands.

## 3. Closed today (for the record)

- T-matrix / Cooper–Wallis routes: 573 (T(191) dead by T-form row sums), 537 (no subgroup), 669 (157 bits; twisted
  T gives no gain).
- Untwisted small T slices: 191/19, 223/37 (SAT UNSAT); 571/57 and 523/29 (row sums infeasible).

## Recommendation

1. Do idea 1 first. It is the only unattempted idea with a plausible n ≤ 100. A twisted estimator plus small-prime
   validation is about an hour of work. If the estimate holds, 571 becomes the best target in the range, and the
   search is k-tree style, so it would suit Kaggle.
2. Keep idea 2 as the Kaggle-time fallback.
3. Ideas 3–5 are not worth compute now.

## Update: twisted estimator for 571 (2026-10-06, later)

**Estimator** (`scripts/gs_twist.py`): multiplier group U of order d, with a signed block permutation g ∈ B₄ satisfying
g^d = 1. It uses 16 sign labels per orbit, no product condition, exact block-sum feasibility (ILP) and a Gaussian
log₂E.

**Validation [M]:**
- **Consistency with `nbits.py`** when untwisted: 523/9 gives 141.9 vs 141.9, and 191/5 gives 79.1 vs 80.2.
- **Brute force** (`logs/tmat/gs_twist_validate_final.txt`): every "infeasible" verdict has a brute-force count of 0.

  | (t, d, twist) | brute force | estimate |
  |---|---|---|
  | (37, 9, 3-cycle) | 2^11.8 | 2^10.9 |
  | (61, 15, 3-cycle) | 2^11.4 | 2^9.6 |

  The estimate runs 1–2 bits low.
- **A bug caught by this validation:** block-sum targets were generated over even integers whenever ⌊√(4t)⌋ was even,
  giving false "infeasible" at 37/9. It is fixed. 571 was unaffected (⌊√2284⌋ = 47 is odd).

**Result:** 571, d = 15, g = 3-cycle on (B, C, D), A fixed.
- **Feasible**, with block sums (31, 21, 21, 21), 39 orbits and 19 classes.
- log₂E = 59.4 (2^50.5 per symmetry class), **n = 94.6 bits** [E].
- It is the only feasible g at d = 15.

**Why it is still not a practical target.**
- The twist gives C(x) = B(u⁻¹x) and D(x) = B(u⁻²x). The problem becomes a **2-list** equation,
  PAF_A(c) = −Σⱼ PAF_B(uʲ c), with lopsided lists:
  - A: U-invariant (order 15) with sum 31, about 2^35 choices;
  - B: U³-invariant (order 5) with sum 21, about 2^110 choices.
- Each B matches some A with probability about 2^−54.6.
- Even a complete A table needs about 2^50 B evaluations. Flatness can recover at most about 14 bits, which still
  leaves about 2^40 lookups into a table of roughly 300 GB.
- The repo's measured local-search scaling (about 2^(0.55·n)) gives about 2^52 flips.
- The 4-list independence that made v = 487 (k-tree) and T(103) (support split) tractable is gone.

**Verdict:** the lead is documented as feasible but without an efficient algorithm. It does not beat 523, b = 0, d = 9,
where the four blocks stay independent and the Kaggle k-tree is ready.

## Update: b ≥ 3 screens finished, up to n = 300 bits (2026-10-06, later)

**Tools**
- **`bge3_screen2.py`:** the exact DP screen with three changes:
  - per-block row sums restricted to feasible profiles;
  - chunked Fourier pruning;
  - a verdict printed after each test, so an OOM kill on a later t does not lose earlier results.
- **`bge3_normscreen.py`:** a new exact necessary condition based on the identity
  Σ_k ‖x̄_k‖² = (4v(ct − c + 1) + 4b(c − 1))/(ct) for the compression onto Q × Z_t.
  - Verified exactly on the 479 and 347 cores at four values of t each.
  - It decides cases the DP cannot. Example: for 599 with c = 20 and Q = Z20, group 1, the coset values are
    restricted to {±1, ±27, ±29}, which makes the norm budget of 192 impossible.
  - Positive control: the real 347 configuration passes at t = 1, 3 and 19.
- **Two related facts:**
  - **Q-character decomposition** (verified on 479 and 347). The ψ = 1 component is a complementary integer
    quadruple over Z_h, with Σ PSD = 4v everywhere. Every ψ ≠ 1 has Σ PSD = 4v at ω ≠ 0 and 4b at ω = 0.
  - **The t = 1 norm identity** gives Σ_k Σ_q x̄_k(q)² = (4v + 4b(c−1))/c.

**Results** (every log₂E > 0 configuration with n ≤ 300; raw records in `logs/bge3_normscreen_part{1,2}.jsonl` and
`logs/bge3_screen2_*`)

| v | surviving b ≥ 3 configurations (n bits) | dead (screen that killed it) |
|---|---|---|
| 523 | c = 4, Q = Z4, d = 12 (226): passes DP t = 1, 5 and norm at every t (group 1) | – |
| 571 | **c = 22, Q = Z22, d = 10 (192)**: DP t = 1 passes; norm t = 5 passes (DP t = 5 is memory-bound) | c = 4, all variants (norm t = 1) |
| 573 | c = 14, d = 12 (269, 285): DP t = 1, 2 pass; norm passes every t | – |
| 599 | none at or below 300 | c = 20, Z2×Z10 (DP t = 1); c = 20, Z20, both groups (norm); c = 8, all; c = 12, Z12 |
| 643 | c = 14, d = 12 (232, 234): DP t = 1 passes; norm t = 3, 5, 9, 15 pass | – |
| 647 | none at or below 300 | c = 12 (DP and norm); c = 18, all three variants (DP t = 1); c = 4 (norm) |
| 659 | c = 10, d = 12, r = 70 (242; 2 of 6 groups); c = 4, r = 90 (288) | c = 6, all (DP and norm); c = 10, r = 66; c = 4, r = 76 and 80 |
| 719 | c = 20, Z2×Z10 (227): DP t = 1; norm t = 5, 7. c = 6 (289) | c = 8, all; c = 10, all; c = 18, Z3×Z6 (DP t = 1) |
| 739 | none at or below 300 | c = 10, all three (DP t = 1) |

**Reading.**
- **No b ≥ 3 route at or below 300 bits remains for 599, 647 or 739** (nor for 537 or 669, which had none).
- **The cheapest survivors** are 571/c = 22 (192), then 523/c = 4 (226), 719/c = 20 (227), 643/c = 14 (232),
  659/c = 10 (242) and 573/c = 14 (269).
- **All of these are far beyond demonstrated search capability** (about 63–91 bits). A b ≥ 3 construction for any of
  them needs a structural method, not more screening. The open questions for that are in
  `report/prompt_bge3_opus.md`.

## Check: are the 2026 matrices H(716) and H(892) Williamson-type? (2026-10-09)

**Motivation.** 537 = 3·179 and 669 = 3·223. A Williamson-type quadruple of order 179 or 223, combined with trivial
length-3 T-matrices via Cooper–Wallis, would give H(2148) or H(2676).

**H(892), from the b = 0 core of length 223 [M].**
- For circulant blocks, Williamson type means every pairwise cross-correlation is symmetric about some common shift.
- Tested over all 16 reversal patterns and all per-block shifts: **no pair of blocks has a symmetric cross-correlation
  about any centre.** Negations and a common multiplier cannot change this.
- So the core is not Williamson-type in its circulant form.
- A Williamson-type form of order 223 would need a Z₂₂₃ action on H(892). Sylow conjugacy (assuming 223² ∤ |Aut H|)
  makes any such action conjugate to the core's own, so this is near-conclusive.

**H(716), from the b = 1 bordered core of length 178.**
- It does not plug in directly.
- The decisive test is whether Aut(H(716)) has an element of order 179 with four row cycles. That needs nauty or
  Traces, which could not be built in this container.
- **Untested.**

**Update (2026-10-08, theory-arm session): resolved [M].** Traces was built from the nauty 2.8.8 source bundled in the
pynauty sdist (`theory/hadaut_traces.c`).
- |Aut(H(716))| = 4 = Z₂², with no element of order 179.
- |Aut(H(892))| = 6 = Z₆, with no element of order 223. This confirms the Sylow argument above.
- Neither group contains Q₈, which every Williamson-type (amicable, not necessarily circulant) form requires.
- So neither 2026 matrix is a Williamson-type ingredient in any equivalent form. Details:
  `report/theory_arm_2000_3000.md` §2.

## Paley-coset ansatz (from the 359/419 forensics), 2026-10-09 [M unless marked]

**Forensics, verified independently:**
- The published v = 359 core is *exactly* of this form on its Z₁₁ factor. 116 of its 128 Z₁₁-cosets are ±Legendre(11)
  on Z₁₁* and the other 12 are constant on Z₁₁*; random signs give about 0.2 Legendre cosets.
- 419 is approximately of this form: 99 twisted and 16 constant cosets out of 152.
- The half-shift template of 167 and 179 is also confirmed (combined log₁₀ P = −71.4 and −78.0 on the repo's own cores).

**Ansatz** (`scripts/paley_ansatz.py`). Take G = G′ × Z_p. On each coset, x is α at z = 0 and, on z ≠ 0, either ±χ(z)
(twisted) or ±1 (constant). That gives 3 bits per coset.
- The core condition reduces **exactly** to 3 equations per h ∈ G′: lags (h, 0), (h, square) and (h, non-square).
  Validated against FFT on random labellings.
- **Exhaustive search** finds solutions at v = 15, 21, 33, 35 (b = 0) and 29 (b = 1); these are saved as positive
  controls.

**Exact screen** (row sums plus the within-coset equations, which are additive per coset). It agrees with the exact
counts on all 11 calibration cases. For example, it shows why v = 45 (b = 1, p = 11) is empty: the equations force
exactly one constant coset, which makes that block's |row sum| ≥ 6.

**Open orders** (`logs/paley_h0_open.jsonl`):
- **Empty (exact):**
  - b = 0 with G′ = Z₁, and with G′ = Z₃ for 537, 573 and 669;
  - b = 1 for 719 (Z₂), 669 (Z₄), 643 (Z₆), 537 (Z₈), 659 (Z₁₄), 523 (Z₁₈) and 739 (Z₁₈).
- **Pass the screen:**
  - b = 1: 599 (Z₂₆×Z₂₃, Z₄₆×Z₁₃), 571 (Z₃₀×Z₁₉, Z₁₁₄×Z₅, Z₁₉₀×Z₃), 647 (Z₃₄×Z₁₉, Z₃₈×Z₁₇), 573 (Z₄₄×Z₁₃,
    Z₅₂×Z₁₁), 659 (Z₉₄×Z₇) and 643 (Z₂₁₄×Z₃);
  - b = 0 with p = 3: 537, 573 and 669.

**Cost estimates [E]** (row-sum-conditioned Gaussian, within about 2 bits of the exact counts where the screen
passes):
- **Full slices:** expected solutions are positive, log₂E = +28 (599 over Z₂₆×Z₂₃) to +47 (647 over Z₃₈×Z₁₇).
  But the constraint counts are 264–386 bits, so the list-method cost is roughly 2^118 or more.
- **Multiplier-symmetric slices** (u a square mod p, labels constant on u-orbits of G′): the best so far is 571 over
  Z₃₀×Z₁₉ with u = 19 (18 orbits) at log₂E ≈ −10 (202 constraint bits); 599 with u = 3 gives about −3 (99 bits).
- **Anti-symmetric slices** (u a non-square mod p, as in 359) **are exactly impossible for b = 1** [P/M]:
  - twisted signs alternate around each orbit, so Σ_g x(g, z) = W (the constant-coset sum) for every z ≠ 0;
  - hence Σ_h PAF_k(h, 0) = A_k² + (p − 1)·W_k² exactly;
  - with b = 1 row sums (±2, 0, 0, 0), W_k ∈ {−1, 0, 1}, and the total falls short of 4m − 4(n′ − 1). For example,
    at 599 over Z₂₆×Z₂₃ the maximum is 2116 against the 2292 needed.

**Symmetric slices, completed runs (b = 1)** [E]: 571 over Z₃₀×Z₁₉ with u = 19: −9.9 (202 constraint bits); 599 with u = 3:
−3.2 (99); **647 over Z₃₈×Z₁₇ with u = 7: −4.8 (148)**. Every other completed slice is at −15 or below. No slice has
log₂E > 0.

## Paley-coset ansatz for the bordered b ≥ 3 architecture (the form 359 actually has), 2026-10-09 [M unless marked]

**Tool** (`scripts/paley_bge3.py`):
- G = G′ × Z_p with p | h and p ∤ c, so Z_p ≤ H and H = H′ × Z_p with H′ of index c in G′.
- Lag (h′, s) lies in H iff h′ ∈ H′, so the reduced target is −4b on H′ × {0, QR, NQR} ∖ (0, 0) and +4 elsewhere.
- G′ may be any abelian group.
- Two exact screens:
  - **h0**: within-coset equations plus row sums, with Σr² = 4v for b ≥ 3;
  - **qnorm** (new): compression onto Q = G′/H′. Per block, (r, N = Σ_q y_q², s, d) are additive over the cosets, and
    the conditions are Σr² = 4v, ΣN = 4(h + b), Σs = −4b and Σd = 0.

**Positive control: the published 359 core.**
- On G′ = Z₂²×Z₈ with H′ = ⟨(0, 0, 2)⟩, p = 11, it lies exactly in the ansatz (116 twisted and 12 constant cosets).
- Its reduced equations equal the b = 7 target exactly, and it passes both screens.

**Exact decomposition for p ≡ 3 mod 4** (proved by expanding the correlation table; verified on 359 and on solutions
found here):
- Write x = α·δ + u·χ + w·ε per coset, and let Y and X be the symmetric and antisymmetric parts of the (QR, NQR)
  equations. Then:
  - R₀ = PAF_α + (p−1)(PAF_u + PAF_w);
  - Y = −PAF_u + (p−2)·PAF_w + CC(α,w) + CC(w,α);
  - X = CC(α−w, u) − CC(u, α−w).
- So R₀ + (p−1)·Y = PAF(y), where y = α + (p−1)·w is the vector of Z_p-coset sums. The core condition therefore splits
  into two stages:
  - **(stage 1)** a G′-level "compressed core" y ∈ {±1, ±(p−2), ±p} with Σ PAF_y = 4p off H′ and −4bp on H′ ∖ 0;
  - **(stage 2)** the twisted signs u: a ±1 GS-type system Σ PAF_u = (p−2)·PAF_w + CC(α,w) + CC(w,α) − T_s, together
    with the linear condition X = 0.
- In 359 every component is separately two-valued off H′:
  - Σ PAF_u = −4;
  - Σ PAF_α = 44 = 4p;
  - Σ PAF_w = 0;
  - CC(α,w) + CC(w,α) = 0.
- 359's w is **algebraic**:
  - block A's w is a bent function on Z₂², placed on the two Z₈-levels 0 and 6 with opposite signs;
  - blocks B and D have w on a single H′-coset; block C has none.
  - So the authors designed w and searched α and u.

**Cost of 359 in this slice** [E]:
- The full Paley slice has 384 free bits, **303.6 constraint bits** and log₂E = +67.7.
- For the published y, stage 2 alone has 116 free bits against 119.6 constraint bits (log₂E = −3.6). The published y is
  a lucky one, so staging moves the difficulty around but does not remove it.
- **The authors therefore solved a system of about 300 constraint bits.** That is far beyond the 63–91 bits our list
  methods have reached.
- Competitive reading: any open-order Paley slice at or below about 300 bits is within their demonstrated reach. That
  includes the b = 1 slices 599/Z₂₆×Z₂₃ (264 bits, +27.9) and 571/Z₃₀×Z₁₉ (298, +40.9) above.

**Our solvers on the reduced system** (`scripts/paley_ls.py`: tabu, or numba simulated annealing with `--sa 1`;
`scripts/paley_staged.py`: two-stage):
- They find verified cores quickly for bordered v = 15, 27, 31, 39, 47, 59 (n′ ≤ 10) and for b = 0 v = 21, 33, 35.
- They fail at n′ = 12–28: v = 51, 63, 65, 87, which have log₂E ≈ +60 and 70–76 constraint bits [E], even with
  3·10⁸-step anneals.
- Exact compressions onto Q and onto Z_p do not explain those failures; both pass.
- On 359 itself (n′ = 32), plain tabu reaches cost 4800 and staged search stalls in stage 1.
- **The SA wall for bordered Paley systems is about n′ = 12. Every open-order option has n′ ≥ 28.**

**Open orders** (59 options with p ∥ m, p ∤ c, c | n′; `logs/paley_bge3_screen_snapshot.jsonl`, 23 decided so far,
the remaining ones still running):
- **Dead (exact):**
  - every option with p ≥ 29: 599 c = 4, 8, 10, 20; 647 c = 6, 12; 659 c = 4, 6, 22; 719 c = 4, 8, 10, 12, 24;
    739 c = 10; 571 c = 4;
  - p = 23: 719 c = 30, killed by h0; **739 c = 4 and 647 c = 4 (the c = 4 options, whose borders always exist), killed
    by qnorm**;
  - p = 19 and 17: 599 c = 30; 719 c = 36, 40; 647 c = 36.
- Still to screen: the small-p options (n′ ≥ 40), including 599 c = 6 and 669 c = 10 with p = 11 (n′ = 54 and 60). Even
  if they pass, n′ is far beyond our SA wall, and each would need a border matched to its compression.

**Verdict.**
- The 359 route is now understood: a Paley(p) twist with a designed sparse w, then search on a reduced system of about
  300 bits.
- It does not give us a cheap route to an open order. Every open-order option either dies on an exact screen or sits
  at n′ ≥ 28, where our annealer fails even on instances with abundant solutions.
- What would change this is a better solver for the reduced system: a large GPU anneal, or the designed-w-then-α-then-u
  pipeline with a much stronger stage solver. Kaggle GPU time is the natural place to test it, with 359 (n′ = 32) as
  the calibration target before any open order.
