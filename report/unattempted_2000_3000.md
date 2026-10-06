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
