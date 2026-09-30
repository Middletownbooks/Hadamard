# W0 (build validation): PASS

Run 2026-09-30. Raw output: `W0_checks.jsonl`.

| n | brute force N(n) | pipeline, no multiplier reduction | pipeline, reduced then orbit-expanded | same set | stage-A inclusion | product condition |
|---|---|---|---|---|---|---|
| 15 = 3·5 | 24 | 24 | 24 | yes | yes | holds |
| 21 = 3·7 | 96 | 96 | 96 | yes | yes | holds |

(Brute force at n = 25, which is outside W0 and has no coprime split: N = 480, product condition holds.)

**Extra build checks, not pre-registered:**
- **Stage-C completeness.** Summed over every margin pair, the lift counts equal 2^((n+1)/2) exactly: 256, 2,048
  and 262,144 at n = 15, 21 and 35. Every lift has the requested margins, and none is repeated.
- **Sensitivity.** Flipping the sign of ε in the product condition at n = 21 gives 0 solutions (96 with the
  correct sign), so the product-condition key is doing real work.

**Code:** `scripts/wil_common.py`, `scripts/wil_brute.py`, `scripts/wil_search.py`, `scripts/wil_w0.py`,
`scripts/wil_liftcheck.py`.

**Known limitation before 65.** Stage D materialises every (C,D) pair of a combination at once. That is fine for
the controls, but it has to be chunked if the lift lists at 65 are large. W3 measures that first.
