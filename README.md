# Hadamard matrices of previously unknown orders

This repository holds constructions of Hadamard matrices of eight orders between 2000 and 4000. As far as we know,
none of these orders had a construction before. It also holds the search code, logs and write-ups behind them.

Every matrix is given as short generating data: four ±1 sequences, plus a border where one is used. A standalone
verifier checks H Hᵀ = N·I exactly in integer arithmetic.

## Results

| order N = 4v | v | construction | first committed (UTC) | data and verifier |
|---|---|---|---|---|
| 2092 | 523 | bordered Goethals–Seidel (b = 1), cyclic Z₅₂₂; quartic-residue multiplier slice (×181) | 2026-10-14, `5ef41fb` | `pub/H2092/` |
| 2292 | 573 | bordered Goethals–Seidel (b = 1), Z₅₂ × Z₁₁; Paley(11) cosets with an order-3 multiplier | 2026-10-09, `39569f5` | `pub/H2292/` |
| 2588 | 647 | bordered Goethals–Seidel (b = 1), cyclic Z₆₄₆; multiplier ×353 acting on Z₁₇ and Z₁₉ | 2026-10-14, `556c3ca` | `pub/H2588/` |
| 2884 | 721 = 7·103 | Cooper–Wallis: T-matrices of order 103 (from L. Schneider's order-2060 construction) with Williamson matrices of order 7 | 2026-10-06, `866e105` | `cw/` (`python3 scripts/cw_build.py 7`) |
| 3532 | 883 | Goethals–Seidel (b = 0), cyclic Z₈₈₃; order-21 multiplier slice; flat-pool k-tree search | 2026-10-16, `52f6aa1` | `pub/H3532/` |
| 3668 | 917 = 7·131 | Cooper–Wallis: new cyclotomic T-matrices of order 131 with Williamson matrices of order 7 | 2026-10-15, `363a9cb` | `pub/H3668/` |
| 3676 | 919 | Goethals–Seidel (b = 0), cyclic Z₉₁₉; order-51 multiplier with a 3-cycle twist on blocks B, C, D | 2026-10-16, `d3ef318` | `pub/H3676/` |
| 3732 | 933 = 3·311 | Goethals–Seidel (b = 0), Z₃ × Z₃₁₁; index-10 cyclotomic classes on Z₃₁₁ | 2026-10-14, `a5b63b5` | `pub/H3732/` |

**Verify** (Python 3 and numpy only):

```
cd pub/H2092 && python3 verify.py      # likewise in pub/H2292, H2588, H3532, H3668, H3676, H3732
python3 scripts/cw_build.py 7          # H(2884)
```

Each verifier:
1. checks that the four sequences have the required periodic autocorrelation sums;
2. builds the full matrix (the Goethals–Seidel array, plus the recorded border for b = 1);
3. checks H Hᵀ = N·I exactly;
4. checks the published normalized matrix the same way.

Every result was also checked by at least one further script that shares no code with the search
(`scripts/verify_h*.py`, `scripts/verify_b0_bundle_rebuild.py`, `scripts/verify_h2092_rebuild.py`; logs in
`logs/verify/`).

**Which orders were open.** The orders below 3000 were open on the list in `report/unattempted_2000_3000.md`. The
3000–4000 orders were open on a list compiled from public sources; SageMath's `hadamard_matrix(n, existence=True)`
also reports them as unknown (`data/sage_unknown_orders_le4000.json`). We would welcome corrections.

## How they were found

**Authorship.** The searches, the mathematics behind them and the code were done by Claude (Anthropic's Opus 5.5
models) in Claude Code sessions, directed and checked by the repository owner.

**Background.** Two 2026 constructions set the context:
- Alpöge, Voinov, Reynolds-Haertle and Claude, every open order below 2000;
- L. Schneider, order 2060.

Our analysis of their published matrices is in `report/alpoge_methods_forensics.md`.

**Methods, in brief:**
- **Staged multiplier slices (b = 1).** Choose a multiplier group that fixes one factor of the group pointwise and
  acts on a prime factor by cyclotomic classes. Compressing along the prime gives a small, exactly solvable first
  stage; the rest becomes norm equations over Gauss periods, solved by a hierarchical join. Records:
  - `report/RESULT_2292.md`;
  - `report/paley_523_answer.md`;
  - `report/paley_647_answer.md`.
- **b = 0 (Goethals–Seidel):**
  - the same staging at composite v (`report/alpoge_design_answer.md` §4);
  - a twisted multiplier slice (`report/b0_919_answer.md`);
  - pools of locally optimal "flat" blocks joined by a 4-list k-tree (`report/b0_883_answer.md`,
    `report/b0_flatpool_ladder.md`).
- **Cooper–Wallis** with T-matrices: `cw/RESULT_2884.md`, `report/b0_917_answer.md`, `report/tmat_routes.md`.
- **Exact screens and obstructions** for many other slices, so negative results are recorded too:
  - `report/slice_screens_20261013.md`;
  - `report/next_design_answer.md`;
  - `report/alpoge_design_answer.md`;
  - `report/bge3_stage1_answer.md`.

## Repository map

| path | contents |
|---|---|
| `pub/` | publication bundles (sequences, border, normalized matrix, standalone verifier) |
| `report/` | write-ups: results, methods, estimates and their calibration, obstructions, and the prompts given to each session |
| `theory/`, `scripts/` | search and analysis code |
| `data/` | cores and packed matrices |
| `logs/` | run logs, including the logs of the runs that found each matrix |
| `forensics/` | analysis of the published Alpöge et al. matrices |
| `cw/` | the H(2884) Cooper–Wallis construction |

Work in progress (open orders 537, 571, 599, 643, 659, 669, 719, 739, 751, 789, 823, 839, 859, 863, 907, 947, 955,
971 and 991) is described in the most recent files of `report/`.
