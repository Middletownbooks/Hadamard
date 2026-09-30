# Williamson controls (W0–W3, W2-T63): results

Closed 2026-09-30. Pre-registration and amendments A1–A3: `PREREGISTRATION.md`.

| Bar | Verdict | Evidence |
|---|---|---|
| W0: pipeline = independent brute force (n = 15, 21) | **PASS** | Identical solution sets: 24 and 96. Multiplier reduction on and off; stage-A inclusion; product condition holds. `W0_result.md` |
| W1: n = 35 (known nonexistence) | **PASS** | 0 solutions, every stage exact, ledger complete. `logs/wil_n35.json` |
| W2 (amended A3): blind hits at n = 45, 51, 55, 57 | **PASS** | Orbit reps / quadruples: 45: 2 / 12, 51: 4 / 32, 55: 2 / 20, 57: 6 / 54. All independently verified; every ledger complete. `logs/wil_n{45,51,55,57}.json` |
| W2-T63: targeted check on Turyn's order-63 quadruple | **PASS** | (a) compressions in stage A; (b) canonical combination kept, blocks liftable (252, 252, 5,088, 3,858); (c) stage D returns 6 verified solutions, including a multiplier image of the constructed one. `W2_T63_result.json`, `turyn63.json` |
| W3: cost gate for 65 (≤ 500 core-hours) | **PASS** | Exact work 8.17×10^13 pair operations over 14.25M combinations. Final kernel, loaded: 7.2×10^7/s at n = 63 list sizes, giving about 316 core-hours; 4.7×10^7/s at n = 55, giving about 479 core-hours. `logs/wil_n65_cost.json` |

**W4 preconditions already met:**
- The independent second stage-A implementation matches the main one at n = 65, on 5 cases × 2 moduli
  (`W4_stageA_crosscheck_n65.jsonl`).
- No caps exist in any stage.

Still to come from the n = 65 run: the completeness ledger (`scripts/wil_collect65.py`).

**Not claimed:**
- n = 63 was not searched blind. A partial run covered about 2% of its first case with no solution
  (`logs/wil/n63_partial_A3.log`); that says nothing either way.
- The Turyn construction was recovered by searching a small pool of trace-character sequences. The independent
  verifier certifies the quadruple itself, whatever its provenance.
