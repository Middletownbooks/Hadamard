# Review of `report/bge3_answer.md` (higher-reasoning session), 2026-10-06

## Independently verified here [M]

- **Border–compression coupling** (answer §4.1). On the verified v = 47 rung core with 347's border, corner and layout:
  - the original core gives a Hadamard matrix;
  - every H-translation of any block preserves the core condition and still works, because the compression is
    unchanged;
  - Q-translations of blocks B, C or D (6 of 9 tested), and negations of B, C or D, preserve the core condition
    (Σ PAF = −20 on H∖{0} and +4 off H) **but break H·Hᵀ = 4v·I**;
  - block A tolerated every move tested (all three shifts and negation).

  So a border must be matched to the compression placement. This corrects `bge3_construction.md` §2.
- **Screen verdicts:** where the two analyses overlap they agree. The answer's `bge3_cost2.py` and this session's
  `bge3_screen2.py` / `bge3_normscreen.py` agree on 647/c12, 659/c6, 599/c20 (both Q), 719/c8, 719/c10, 739/c10 and
  571/c22.

## Items where this session has newer results than the answer

| configuration | answer says | this session (DP t = 1 screen, `bge3_screen2.py`) |
|---|---|---|
| 647 c18 Z18 d12 r74 | "DP overflow" | **dead**: t = 1 fails in both groups |
| 647 c18 Z18 r72 and Z3×Z6 r72 | survivor of the norm screen | **dead**: t = 1 fails in all groups |
| 719 c18 Z3×Z6 r78 | survivor of the norm screen | **dead**: t = 1 fails in all 3 groups |
| 659 c10 r70 | survivor of the norm screen | **alive**: t = 1 passes in 2 of 6 groups (n ≈ 242) |
| 659 c4 r90, 719 c6 r84 | survivors of the norm screen | alive at t = 1 (n ≈ 288–289) |

## Accepted as reported (not independently re-run)

These rest on the answer's own exact-count scripts and positive controls; I did not re-run them:
- the calibration counts (§2.2);
- the 182 exact slice decisions (§3.3);
- the refined-model numbers (571: n_lift 168, log₂E 23.5; 523 c4: log₂E −32; 643: log₂E 2);
- the forensics (§5).

The list-method barrier (§2.1) and the decoupling cost argument (§2.4) are derivations that I find consistent with
this repo's measurements: T(103) took 12 min, while the T(191) and twisted-571 analyses hit the same ≈ 2r-bit wall.

## Net status of b ≥ 3 for open orders

- **No b ≥ 3 configuration at or below 300 bits survives for 537, 599, 647, 669 or 739.**
- **Cheapest survivor:** 571/c22, which needs ≳ 2^80 operations and also a border matched to its compression (none
  known for c ≥ 10).
- **Everything else is at ≳ 2^93 or worse:** 523/c4, 719/c20, 643/c14, 659/c10 and 573/c14.
- **Conclusion (agreed by both analyses):** no Kaggle or local compute on b ≥ 3 for open v. 523/b = 0/d = 9 stays the
  primary route.

## Leads from the answer worth a cheap follow-up

1. **The generalized-Ψ architecture (answer §6)** has b < c − 1 with a border that always exists. For example c = 8,
   b = 3, Ψ a basis of Q = Z2³, which covers 523 (h = 65), 571 (71), 643 (80), 659 (82) and 739 (92). It is
   unconfirmed: there is no core at h = 1, 3 or 5. The next step is the Gaussian cost estimate with this target.
2. **Reproduce a solution in 479's own slice** (d = 12, r = 46, log₂E ≈ 11, T_min ≈ 2^71 by list methods). This is the
   test for any method that claims to beat the barrier.
