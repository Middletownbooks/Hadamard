# Positive control for the c=6 compression-and-lift pipeline: results

Closed 2026-09-30. Pre-registration: `PREREGISTRATION.md` (PC1–PC4, addendum PC1b, addendum 2 PC0). Per-slice
records: `slice_*.json` in this folder. Logs: `logs/posctrl/`.

## Verdicts

| Bar | Verdict | Why |
|---|---|---|
| PC0: known v=47 slice found blind | **FAIL (unevaluable)** | The exact join did not finish (see `PC0_result.md`). Targeted checks on the known solution all pass. |
| PC1: ≥ 1 matrix of order 380 (v=95) | **FAIL** | 0 found. 2 slices exhaustively empty, 1 slice with 837,168 compressions of which 206,481 were lifted exactly with none liftable, and the rest too large. |
| PC1b: ≥ 1 matrix at v=131 or 215 | **FAIL** | 0 found in any slice (table below). |
| PC2: no false rejections by the residue filters | **FAIL (unevaluable)** | No blind slice produced a liftable compression, so the bar was never tested. The only evidence is the PC0 targeted check (the known solution passes both filters). |
| PC3: fast lift = exact join on liftable compressions | **FAIL (unevaluable)** | No liftable compressions. Outside the registration, wc_lift matched the exact join on the real 347 core and on 515 sol 3. |
| PC4: negative test rejected | **FAIL (unevaluable)** | It only runs on a found solution, and none was found. |

Reading: **the blind positive control did not validate the pipeline.** No slice produced a solution, and every
bar that depends on one could not be evaluated. What does hold is component-level: stage-A inclusion, the filters
and assembly on the one known solution (PC0 targeted checks), and exact agreement of the fast lift with the exact
join on the real 347 core. This is weaker than a blind hit, so the 515 null stays an "incomplete search", not
evidence of absence.

## Slice table

| v (h) | p | u | compressions | lifted (exact) | liftable | solutions | status |
|---|---|---|---|---|---|---|---|
| 95 (15) | 3 | 7 | 0 | – | – | 0 | exhaustive, empty |
| 95 (15) | 3 | 17 | 0 | – | – | 0 | exhaustive, empty |
| 95 (15) | 3 | 11 | – | – | – | – | stage A too large |
| 95 (15) | 5 | 17 | 837,168 | 206,481 | 0 | 0 | lift partial (container restart) |
| 95 (15) | 5 | 7, 11, 19 | – | – | – | – | too large (stage A or lift) |
| 131 (21) | 3 | 5 | 0 | – | – | 0 | exhaustive, empty |
| 131 (21) | 3 | 19 | 0 | – | – | 0 | exhaustive, empty |
| 131 (21) | 3 | 11 | 2,283,136 | 402 | 0 | 0 | sampled; all 402 pass the filters, none lifts |
| 131 (21) | 7 | 5 | 1,798,464 | 20,000 | 0 | 0 | sampled |
| 131 (21) | 7 | 11 | 1,798,464 | 20,000 | 0 | 0 | sampled |
| 131 (21) | 7 | 19 | – | – | – | – | too large |
| 215 (35) | 5 | 17 | 0 | – | – | 0 | exhaustive, empty |
| 215 (35) | 5 | 73 | 0 | – | – | 0 | exhaustive, empty |
| 215 (35) | 5 | 37 | – | – | – | – | stage-A join too large |
| 215 (35) | 5 | 23 | – | – | – | – | OOM-killed twice |
| 215 (35) | 7 | 17 | 18,295,008 | 39,107 | 0 | 0 | sampled |
| 215 (35) | 7 | 23 | 18,295,008 | 38,297 | 0 | 0 | sampled |
| 215 (35) | 7 | 37 | 0 | – | – | 0 | exhaustive, empty |
| 215 (35) | 7 | 73 | 0 | – | – | 0 | exhaustive, empty |
| 215 (35) | 7 | 11, 19, 59, 79, 101, 149 | – | – | – | – | OOM-killed in stage A (18–30 GB; container limit 15 GB) |

Other attempts: annealing at h = 15 with d = 2 did not reach a solution (best energy about 672).

## Contrast with the Williamson pipeline

The Williamson search (`williamson/`) has what this control lacked: orders where solutions are known to exist
and are small enough to enumerate exhaustively. Its W0 and W2 controls reproduce known solution sets exactly.
