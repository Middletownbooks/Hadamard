# Gate A of the OpenAI fresh-method proposal: run 2026-11-07, failed

This is the blind calibration gate of `report/openai_fresh_answer.md` §4 (branch `codex/order-3452-compression`,
commit d152756), run exactly as specified.

- **Problem:** b = 0 at p = 61, d = 3, rows (11, 11, 1, 1).
- **Method:** parametric cycle oracle with spectral concurrence, beta = 1.
- **Budget:** seeds 1001–1008, 4 trials each, 10,000 iterations per trial, check interval 10.

Each trial took about 19 s, four processes at a time.

**Result [M]: 0 of 32 trials found a core.** The pass criterion was at least 16 of 32.
`theory/openai_fresh_verify.py` verified 0 matrices, and no core outside the requested row sums was recorded.
The best integer PAF energy reached per trial was 96 in 2 trials, 192 in 27 and 288 in 3. Zero is a core.

**Decision:** by the proposal's own rule ("stop at the first failed gate"), the outer iteration is not scaled.
Gates B and C, and the 4412 pilot, were not run. The exact per-lag cycle oracle (§1.2 of the answer) remains a
valid, verified projection; only the search built on it fails.
