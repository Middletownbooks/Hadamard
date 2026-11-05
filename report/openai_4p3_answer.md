# Algebraic proposals for the k4 bordered family at n = 4p + 3

Status: initial checkpoint, 2026-11-05. Work in progress; no new core or Hadamard matrix is claimed.

Tags: [P] proof supplied; [M] exact computation with the script and scope named; [E] estimate with method named. Untested proposals are explicitly conditional.

Credit: the bordered construction is Theorem 7 of Egan–Kotsireas–London–Ó Catháin (EKLO), in Kumar's form (https://github.com/vatskishan03/hadamard). The census, pre-existing obstructions, calibration and border tooling are by Claude Code in this repository.

Setup: GitHub's comparison found codex/order-3452-compression at 5b26dae60414599ccccf4c692220f616fc2aae24 to be an ancestor of claude/focused-planck-j1645f at a87392b4d91fcaafe0bc61d85435b4ff0a06c878 (nine commits ahead). The target branch was fast-forwarded to that reference commit with an expected-head check. The shell execution environment was unavailable, so this used the GitHub API, not a local git checkout. No PR is opened and no other branch is written.

## Bottom line

Pending. The session will prioritize rigorous obstructions and explicit reduced constructions. Coverage of all 16 target orders is currently unknown.

## The explicit coset system

[P] Identify G = Q × F_p and put
C_i(q,r;t) = sum over k in F_p of a_(i,q)(k) a_(i,r)(k+t).
Then the necessary and sufficient equations are

- sum_(i,q) C_i(q,q;0) = 16p;
- sum_(i,q) C_i(q,q;t) = -12 for every t != 0;
- sum_(i,q) C_i(q,q+d;t) = 4 for every d != 0 and every t.

Proof: partition the defining PAF sum over G into its four K-cosets. This is an identity, not an additional ansatz.

For Q = Z_2^2 label its elements 0,1,2,3 with addition XOR. The three cross equations, summed over i, are respectively
C(0,1)+C(1,0)+C(2,3)+C(3,2)=4,
C(0,2)+C(2,0)+C(1,3)+C(3,1)=4,
C(0,3)+C(3,0)+C(1,2)+C(2,1)=4,
at each t. For Q = Z_4 they are
C(0,1)+C(1,2)+C(2,3)+C(3,0)=4,
C(0,2)+C(1,3)+C(2,0)+C(3,1)=4,
C(0,3)+C(1,0)+C(2,1)+C(3,2)=4.
Here C means C_i at the common lag t.

[P] Writing z_iq = sum_k a_iq(k), summing these equations over t gives
sum_iq z_iq^2 = 4p+12 and sum_iq z_iq z_i,q+d = 4p for d != 0.
In particular, an ansatz with every |z_iq|=1 is impossible for every odd prime p: its energy is 16 instead of 4p+12. This immediately excludes using only translated/signed zero-completed quadratic characters in all 16 cells, without repairs.

## Ranked proposals

Pending: definitions, proofs or reduced conditions, exact checks, 16-order coverage and fixed pass/kill criteria will be recorded here.

## Ideas considered and rejected

- Pure translated or signed zero-completed quadratic characters in every coset: compression-energy obstruction above [P].

## What a Claude instance should run first

Pending scripts. Do not interpret this checkpoint as a validated construction.

## Unfinished items

- Read calibration and border interfaces.
- Complete and rank up to three proposals.
- Supply the builder, independent integer PAF verifier and bounded tests.
- Check all 16 target primes.
- Obtain a small core and independent bordered matrix check if feasible.
- Python execution is not currently available in this session; any checks run in the JavaScript orchestration runtime will be labelled as such and their reproducible source committed.
