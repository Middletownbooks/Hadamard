# Fresh-method proposals for open Hadamard orders

Status: early checkpoint, 2026-11-06 (user-local date). Work in progress.
No new Hadamard matrix is claimed. Rankings and cost estimates below are pending
the repository audit and short experiments.

## Scope and evidence labels

This answers `report/prompt_openai_fresh_idea.md`. [P] denotes a proof, [M] an
exact computation with a named script, and [E] an estimate or testable hypothesis.
Only `codex/order-3452-compression` is a push target; no pull request is requested.
The current Claude branch was fetched and merged cleanly before this checkpoint.

## Ranked proposals (unfinished)

At most three proposals will survive review. Each will specify:

1. The precise array, group, ansatz, and algorithm.
2. The assumption of earlier failed methods that changes.
3. Work per solution, assumptions, and comparison with a 2^46-operation budget.
4. Blind calibration orders and predetermined pass/kill criteria.
5. A complete first-experiment protocol for Claude Code.
6. Short experiments, with `theory/openai_fresh_*.py` and
   `logs/openai_fresh/` artifacts where applicable.

Candidate directions under audit: an exact reformulation of the half-shift
template; spectral constraints with discrete recovery; and ingredients or rich
hosts with a substantively different search state. These are not yet proposals.

## Reconstruction of the August methods (unfinished)

The originals and existing forensic scripts will be examined. Statistical
evidence will be distinguished from recovery of an actual historical algorithm.

## Ideas considered and rejected (unfinished)

The audit will explicitly exclude reruns of closed exact multiplier slices,
closed algebraic families, and parameter-only variations of existing solvers.

## Reproduction and verification (unfinished)

All short runs must remain below ten minutes each. A claimed new matrix would
require a separate exact-integer check of H H^T = n I, independent of construction.

## Credit and sources

- Egan, Kotsireas, London and Ó Catháin (EKLO): survey, Theorem 7, flatness analysis.
- Alpöge, Voinov and Reynolds-Haertle: the twelve announced matrices.
- Kishan Kumar, https://github.com/vatskishan03/hadamard (Apache-2.0): constructions
  and ingredient routes.
- This repository: previous constructions, exact screens, obstructions, and
  forensic data. Source provenance will be expanded in the final report.
