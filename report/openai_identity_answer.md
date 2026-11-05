# Algebraic identities for the odd-|K| k4 bordered family

Session: 2026-11-05. Authorized push branch: `codex/order-3452-compression`.
Claim tags: [P] proof, [M] exact computation with named script, [E] estimate.

## Bottom line

Initial checkpoint. No new family or new Hadamard matrix is claimed. The objective is an algebraic construction or a rigorous obstruction for a precisely defined natural class. Coverage of the 38 supplied targets is currently unknown.

## Analysis of the known cores

Unfinished: inspect the 1916 and 3964 cores for two-prime multiplicative structure; compare their partial compressions with small examples.

## Proposal 1: quadratic-character and two-prime identities

Unfinished: define the class, derive exact character/augmentation constraints, and implement `theory/openai_identity_build.py`. An unsuccessful finite search is not a nonexistence proof outside its explicitly enumerated class.

## Proposal 2: product identities and relative difference sets

Unfinished: derive the required group-ring identity and test whether standard norm products preserve it. A proposal must eliminate the correlation conditions algebraically.

## Ideas considered and rejected

The previously reviewed small cyclotomic-repair and elliptic search families remain closed, as specified in the prompt.

## What a Claude instance should run first

Unfinished. Exact scripts, commands, input formats, and pass/kill criteria will be added with each completed result. Any successful small core must also pass independent integer bordered-matrix verification.

## Unfinished items

- All new mathematical conclusions and computations.
- Known-core structural analysis, small-instance tables, and 38-target coverage.
- Independent matrix verification for any constructed control.

## Provenance and setup

The bordered construction is Theorem 7 of Egan–Kotsireas–London–Ó Catháin, in Kumar's form. Theorem B and the compression normal form are the prior OpenAI results reviewed in `report/openai_4p3_review.md`.

The two requested branches were fetched. A normal merge preserving both branch histories is being used. No pull request will be opened; only the authorized Codex branch will be pushed.
