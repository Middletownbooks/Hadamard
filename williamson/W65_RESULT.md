# Williamson matrices of order 65: result

Completed 2026-10-05. Pre-registration: `PREREGISTRATION.md` (W0–W4, amendments A1–A3). Controls:
`CONTROLS_RESULT.md`.

## Result

**The exhaustive search finds no Williamson quadruple of order 65.** That is, no four symmetric ±1 circulants
A, B, C, D of order 65 have AAᵀ + BBᵀ + CCᵀ + DDᵀ = 260·I.

| Row-sum case | Combinations (exact) | Pair operations (exact) | Shards | Solutions |
|---|---|---|---|---|
| (9,9,7,7) | 2,880,500 | 1.723×10¹³ | 32/32 | 0 |
| (11,9,7,3) | 2,756,880 | 1.569×10¹³ | 32/32 | 0 |
| (11,11,3,3) | 2,626,820 | 1.438×10¹³ | 32/32 | 0 |
| (13,9,3,1) | 2,980,740 | 1.822×10¹³ | 32/32 | 0 |
| (15,5,3,1) | 3,008,536 | 1.622×10¹³ | 32/32 | 0 |
| **Total** | **14,253,476** | **8.174×10¹³** | **160/160** | **0** |

Every case's ledger and work equal, exactly, the totals computed independently by the cost pass
(`logs/wil_n65_cost.json`). Collector output: `n65_collected.json`.

## W4 conditions for a nonexistence claim (pre-registered)

| Condition | Status |
|---|---|
| W0, W1, W2 (amended) and W2-T63 all PASS | yes (`CONTROLS_RESULT.md`) |
| No stage hit a cap or truncation | yes. The code has no caps; all 107 local shard records show `caps_hit: false`. The lift cache is bounded but only recomputes, and it was checked to be result-neutral. |
| Stage-A lists at 65 match a second, independent implementation | yes (`W4_stageA_crosscheck_n65.jsonl`, 10/10) |
| Completeness ledger: a stage-D record for every canonical combination | yes (`n65_collected.json`, `all_complete: true`) |

**Verdict under the pre-registration: nonexistence claimed.**

## Where the evidence comes from

- **107 shards ran in this container.** Their full JSON records are in `logs/wil65/`.
- **53 shards ran on Kaggle**, recorded from the notebooks' per-shard summary lines (`kaggle_out/wil65/kaggle_logs.txt`).
  - Each line's ledger was checked against the exact expected size for its (case, shard).
  - The per-case work sums equal the independent totals exactly.
  - One Kaggle coincidence (equal work in case 0 shards 9 and 10) was recomputed here and confirmed.

## Caveats (what this result rests on)

- **"Williamson" means four symmetric ±1 *circulant* matrices**, the standard computational definition.
- **The reductions used are exact, but their correctness rests on this pipeline.** They are:
  - normalisation (positive, sorted row sums);
  - (Z₅)*×(Z₁₃)* multiplier reduction;
  - the product condition aⱼbⱼcⱼdⱼ = ε;
  - pair-PSD pruning.

  Each was validated against brute force at n = 15 and 21, by blind recovery of known solutions at
  n = 45, 51, 55 and 57, by the targeted Turyn check at 63, and by 0 solutions at 35.
- **An independent reimplementation would strengthen the claim before publication.** That could be a SAT/CAS
  encoding of the 5 row-sum cases, or a different search order. It was not pre-registered as a condition and has
  not been done.
