# Hadamard matrices from Claude Code searches

This repository holds Hadamard matrices of eight orders between 2000 and 4000, the search code that found them, and
the write-ups, logs and negative results behind them.

**Status, following Egan, Kotsireas, London and Ó Catháin [EKLO].** Three of the orders, 2092, 2292 and 2588, were
open in the published literature and were first constructed here. They are credited to this repository in [EKLO].
The other five were constructed here independently of earlier or concurrent work, listed below.

Every matrix is given as short generating data: four ±1 sequences, plus a border where one is used. A standalone
verifier checks H Hᵀ = N·I exactly in integer arithmetic.

## Results

| order N = 4v | v | construction here | first committed here (UTC) | status | data and verifier |
|---|---|---|---|---|---|
| 2092 | 523 | bordered Goethals–Seidel (b = 1, a Legendre quadruple of length 522); quartic-residue multiplier slice (×181) | 2026-10-14, `5ef41fb` | **first construction** [EKLO]; an independent matrix of this order, by a different construction (an SDS in ℤ₅₂₃ with a multiplier group of order 9), is in [K] ([hadamard/2092](https://github.com/vatskishan03/Hadamard/tree/main/hadamard/2092)) | `pub/H2092/` |
| 2292 | 573 | bordered Goethals–Seidel (b = 1), ℤ₅₂ × ℤ₁₁; Paley(11) cosets with an order-3 multiplier | 2026-10-09, `39569f5` | **first construction** [EKLO]; [K] has a later independent matrix | `pub/H2292/` |
| 2588 | 647 | bordered Goethals–Seidel (b = 1), cyclic ℤ₆₄₆; multiplier ×353 acting on ℤ₁₇ and ℤ₁₉ | 2026-10-14, `556c3ca` | **first construction** [EKLO] | `pub/H2588/` |
| 3676 | 919 | Goethals–Seidel (b = 0), cyclic ℤ₉₁₉; order-51 multiplier with a 3-cycle twist on blocks B, C, D | 2026-10-16, `d3ef318` | independent; Schumer gave a skew Hadamard matrix of this order earlier [S] (cited with this repository in [EKLO]) | `pub/H3676/` |
| 3732 | 933 = 3·311 | Goethals–Seidel (b = 0), ℤ₃ × ℤ₃₁₁; index-10 cyclotomic classes on ℤ₃₁₁ | 2026-10-14, `a5b63b5` | independent; found concurrently by the authors of [EKLO] and in [K] | `pub/H3732/` |
| 3532 | 883 | Goethals–Seidel (b = 0), cyclic ℤ₈₈₃; order-21 multiplier slice; flat-pool k-tree search | 2026-10-16, `52f6aa1` | independent; the order was already settled (Đoković, arXiv:1008.2043, 2010) | `pub/H3532/` |
| 3668 | 917 = 7·131 | Cooper–Wallis: cyclotomic T-matrices of order 131 found here, with Williamson matrices of order 7 | 2026-10-15, `363a9cb` | independent; already obtainable classically (T-matrices of order 131 from London and Kotsireas's TT(44), 2025, with Cooper–Wallis) | `pub/H3668/` |
| 2884 | 721 = 7·103 | Cooper–Wallis: L. Schneider's T-matrices of order 103 with Williamson matrices of order 7 | 2026-10-06, `866e105` | an immediate consequence of Schneider's published T-matrices; not open in [EKLO] | `cw/` (`python3 scripts/cw_build.py 7`) |

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

**Open orders.** By Table 5 of [EKLO], the smallest order with no known Hadamard matrix is now 2396 = 4·599. The
open orders below 4000 are 4v for v = 599, 643, 669, 719, 739, 751, 839, 859, 863, 907, 947, 955 and 971. Our earlier
working lists (`report/unattempted_2000_3000.md`, `data/sage_unknown_orders_le4000.json`) predate the 2026 results
and should not be read as current.

**References.**
- [EKLO] R. Egan, I. S. Kotsireas, S. London and P. Ó Catháin, *Recent progress in the construction of Hadamard
  matrices* (2026). Their matrices are at https://www.math74.com/Hadamard/.
- [K] K. Kumar, *New Hadamard matrices*, https://github.com/vatskishan03/Hadamard (2026).
- [S] A. Schumer, *Periodic T-matrices of order 151 and a skew spin Goethals–Seidel family over ℤ/919ℤ*, Zenodo,
  doi:10.5281/zenodo.22820997 (17 September 2026).

## How they were found

**Authorship.** The searches, the mathematics behind them and the code were done by Claude (Anthropic's Opus 5.5
models) in Claude Code sessions, directed and checked by the repository owner.

**Background.** Two 2026 constructions set the context:
- Alpöge, Voinov and Reynolds-Haertle (announced 12 August 2026), every open order below 2000;
- L. Schneider, order 2060.

Related work of 2026 on the next orders includes [EKLO] and [K]; in particular [K] constructs 2148, 2636 and 3964 with
bordered Goethals–Seidel arrays of border width 12 over non-cyclic groups.

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
  - `report/bge3_stage1_answer.md`;
  - `report/b599x4_answer.md`: a coset-Parseval obstruction. A multiplier slice has no core when the orbits that
    are not unions of cosets of some subgroup N of H cover fewer than v/(1 + 1/|N|) points (|N| odd; (1 − 1/|N|)·v
    for |N| even). It empties the 599 ×4 slice and 259 of the generated configurations in all.

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

Work in progress on the open orders is described in the most recent files of `report/`.
