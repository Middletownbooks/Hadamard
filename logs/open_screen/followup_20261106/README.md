# Follow-up to the open-order screen (2026-11-06)

The screen's step 5 took candidates only from tier 0. Tier 0 required at least 2^10 expected cores per stage-1
solution y. That rule left out cheap tier-2 slices with log2 E_tot >= 0 and estimated cost <= 2^30 per core. These
runs apply `theory/open_screen_stage1_exact.py KEY --methods A,B --strong --planted 2` to them, with 900 s per slice
and two at a time (keys in `keys.txt`).

| 4v | fam | slice key | result |
|---|---|---|---|
| 6268 | k4 | ["k4", 1567, [4], [17], 23, 2, [[1, 16]]] | stage 1 **empty** by A and B, planted 3/3 by both [M] |
| 4412 | b1 | ["b1", 1103, [], [19, 2], 29, 4, [[0, 9]]] | stage 1 **empty** by A and B, planted 2/2 by both [M] |
| 4892 | b1 | ["b1", 1223, [], [47, 2], 13, 4, [[0, 23]]] | stage 1 **non-empty**: method B finds 384 solutions in 55 s (20 stored); A hit MemoryError; the planted phase did not finish within 900 s |
| 7180 | b1 | ["b1", 1795, [], [13, 3, 2], 23, 2, [[0, 12]]] | undecided (900 s timeout) |
| 3004 | k4 | ["k4", 751, [4], [17], 11, 2, [[1, 8]]] | undecided (900 s timeout) |
| 5164 | k4 | ["k4", 1291, [4], [7, 2], 23, 2, [[1, 6]]] | undecided (900 s timeout) |
| 7916 | k4 | ["k4", 1979, [4], [13, 2], 19, 2, [[1, 12]]] | undecided (900 s timeout) |
| 4792 | b1 | ["b1", 1198, [], [9, 7], 19, 2, [[1, 6]]] | undecided (900 s timeout) |
| 5892 | b1 | ["b1", 1473, [], [32, 2], 23, 2, [[0, 8, 3]]] | tool error: multi-generator U on one axis (`b_setup`, `5892_b1_error.log`) |

**Check of the 4892 solutions [M].** Five of the stored solutions were checked by a separate script (index a*2 + b on
Z47 x Z2). It used the stage-1 target of the b = 1 family compressed along Z13: 4*1222 - 4*12 = 4840 at 0 and
-4*13 = -52 at every other lag. All 94 lag equations hold, the row sums are (-2, 0, 0, 0), and every entry lies in the
alphabet +-1 + 3*{-4, -2, 0, 2, 4}. U-invariance was not checked.

**4892 screen numbers [E]** (`table.jsonl`): stage-1 log2E 14.79; E per y = -12.46; C2 = 58.31; best split G'
with D_def 17.3, h_D 16.84; cost 2^27.5 per core; log2 E_tot = 2.33.
