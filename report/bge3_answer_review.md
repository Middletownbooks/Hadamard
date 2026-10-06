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

## Follow-up: cost of the generalized-Ψ architecture (answer §6), 2026-10-06

**Tool:** `scripts/bge3_cost_psi.py`, which is `bge3_cost.py` with three changes:
- the target is Σ_k PAF_k(q, z) = −4·Σ_{ψ∈Ψ} ψ(q) off 0, with b = |Ψ|;
- only symmetry groups that preserve the target are kept;
- the Gaussian model is unchanged.

**Validation [M]:** with Ψ = all nontrivial characters it reproduces `bge3_cost.py` exactly on 347 (n = 200.9,
log₂E = 55.5).

**Runs** (Q = Z2³, Ψ = coordinate basis, b = 3, c = 8, v = 8h + 3, orbit counts r = 20–200; logs in
`logs/genpsi/*_final.txt`):

| v | h | configurations costed | log₂E > 0 | cheapest | list-method bound 2^((n − log₂E)/2) |
|---|---|---|---|---|---|
| 571 | 71 | 28 | 5 | d = 7, r = 88, L = 47: **n = 253, log₂E = 75.6** | ≈ 2^89 |
| 523 | 65 | 315 (stopped; only larger r left) | 32 | d = 6, r = 108, L = 59: n = 314, log₂E = 93 | ≈ 2^111 |
| 659 | 82 | 45 | 1 | d = 5, r = 144: n = 440, log₂E = 110 | ≈ 2^165 |
| 643 | 80 | 378 | 0 | – | – |
| 739 | 92 | 23 | 0 | – | – |

**Verdict:** the generalized architecture does not open a cheaper route.
- Its best case (571 at about 2^89) is worse than the original architecture's 571/c22 (about 2^80, answer §2.1).
- 643 and 739 have no configuration with expected solutions in range.
- The unscreened numbers are upper bounds on promise: exact screens can only remove configurations.
- Its advantage (a border that always exists) does not offset the higher n, because b = 3 leaves the four blocks with
  a less favourable target than b = c − 1.

**Lead closed for these orders.**
