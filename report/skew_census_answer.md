# New skew-Hadamard orders 4v from the Đoković–Seberry census: coverage, exact screen, calibrated search

Date: 2026-11-07/08. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_skew_census.md`.
Scripts: `theory/skew_*.py`, `theory/skew_*.c`, `scripts/verify_skew.py`. Logs: `logs/skew/`.

**Tags.** [P] proved; [M] exact computation (script named); [E] estimate (estimator named). A cap or a timeout is never
a zero.

**Credit.** The low-energy k-tree is K. Kumar's method (`search/cp_lowe.py`, `search/cp_enum.c`, `search/kt4.c`,
`search/listgen`, Apache-2.0, github.com/vatskishan03/hadamard at d885022). `kt4.c` is compiled and run unchanged;
the list generators `theory/skew_enum_lowe*.c` and `theory/skew_lg*.c` are re-implementations for arbitrary abelian
groups and say so in their headers. `cp_lowe.py` and `cp_hunt.py` were run unchanged from the read-only clone.

<!-- BOTTOM -->
<!-- /BOTTOM -->
