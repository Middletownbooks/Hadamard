# Prompt (high effort, bounded): verify the OpenAI model's n = 4p + 3 results, test its universal border on real cores, and close or continue the elliptic-seed ansatz

*Self-contained, for a fresh Claude Code session on `middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`.*

**What this reviews.** `report/openai_4p3_answer.md` on branch `codex/order-3452-compression` (commit 27cfc29), the
answer to `report/prompt_openai_4p3.md`. It claims:
1. **A compression normal form [P]** with an obstruction to small repairs of cyclotomic seeds [P, M].
2. **A border for every valid core of the k4 family with odd |K|** over Q = ℤ₂² or ℤ₄, from one finite
   certificate [P, M].
3. **An elliptic-character seed ansatz** that gives cores at p = 3, 5, 7. It is exhaustively empty at p = 11, 17 and
   19, and unknown at the 16 target primes.

No open-order matrix is claimed.

**Your job.**
- **Check** claims 1 and 2 independently. If they hold, record the border theorem as a repository result.
- **Test** claim 2 end to end on real cores of other orders.
- **Decide** whether the elliptic ansatz deserves any more compute.

The answer file is `report/openai_4p3_review.md`. Stop after writing it. **Time box:** about 6 hours of wall-clock
time.

**Read first:**
1. this file;
2. `report/openai_4p3_answer.md` (all of it);
3. `report/prompt_openai_4p3.md` §1, the problem;
4. `report/b3_push_answer.md` §1 and §7, the controls and the border problem.

## Rules

- **Setup:** run `git fetch origin claude/focused-planck-j1645f codex/order-3452-compression`, then
  `git merge --no-edit origin/codex/order-3452-compression` into `claude/focused-planck-j1645f`. A trial merge was
  clean. This brings in `theory/openai_4p3_*.py` and `logs/openai_4p3/`.
- **Pushing:** before every push run
  `git fetch -q origin claude/focused-planck-j1645f && git merge -q --no-edit origin/claude/focused-planck-j1645f`,
  then `git push -u origin HEAD:claude/focused-planck-j1645f`.
- **Branches:** never push to `main` or to the codex branch. Open no pull requests.
- **Commits:** trailers follow the existing commits. Commit after each task, with logs in `logs/4p3_review/`.
- **Processes:** never kill processes you did not start. Never `pkill -f` with a pattern that matches your own
  shell; kill by PID.
- **Credit:** the theorems, the ansatz and the scripts `theory/openai_4p3_*.py` are by the OpenAI model; the bordered
  construction is Kumar's form of Theorem 7 of EKLO; the frame-search border reduction is the repository's
  (`theory/bt_border.py`). Say so in the headers of any file you derive from them.
- **Tags:** mark claims [P] proved, [M] exact computation (script named), [E] estimate. A timeout or a cap is never a
  zero.
- **Independence:** a verification counts only if your code shares none with the code it checks. That means no
  imports from `theory/openai_4p3_*.py`, and for matrix checks none from the construction or border code either.
- **Compute:** 4 cores and about 15 GB of RAM. Install `numpy`, `numba` and `ortools` with pip if they are missing.

## 1. The claims in brief

**The family.**
- k4: G of order 4m, K of order m odd, Q = G/K ∈ {ℤ₂², ℤ₄}.
- Target: Σᵢ PAFᵢ = 4|G| at 0, −12 on K \ 0, +4 off K. With n = |G| + 3 this gives H(4n).
- Q-compressions: zᵢ(q) = Σ_{k∈K} xᵢ(q + k), odd. Their targets are Σᵢ PAF_{zᵢ} = 4m + 12 at 0 and 4m off 0, on Q.

**Claim 1, the normal form [P].** After permuting whole rows and negating whole rows, every odd solution has:
- one constant row (a, a, a, a);
- three rows with three entries tⱼ and one entry tⱼ + 2;
- with m = a² + Σⱼ tⱼ(tⱼ + 1).

The converse holds too.
- *Already checked here:* a brute force sharing no code (all odd 4-vectors with energy ≤ 4m + 12, joined 2 + 2 on
  their PAF vectors) found every solution in normal form for every odd m from 3 to 39, on both ℤ₂² and ℤ₄, with no
  violations. The script is not in the repository; redo it yourself.
- *Example:* the 863 Codex compression z = (1,1,1,1), (1,1,1,3), (1,1,1,3), (13,15,15,15) has
  m = 1 + 2 + 2 + 210 = 215 with t = −15 for the last row (after negation).

**Claim 1′, the repair obstruction [P, M].**
- Seeds are cyclotomic cells of fixed index e plus at most c sign changes per cell.
- Minimum total changes at the 16 targets: 29–61. Minimum common per-cell cap: 4–7. Table in §4 of the answer.
- Fixed e and fixed c cannot give an unbounded prime family (explicit inequalities in §2 of the answer).

**Claim 2, the universal border [P, M].**
- Normalise: permute, negate and translate rows on Q so that the exceptional entries sit at q = 0. Then
  z₀ = c₀·1 and zⱼ = cⱼ·1 + 2δ₀.
- One certificate for z* = (0, 2δ₀, 2δ₀, 2δ₀) then works for every core: X (12 × 12), Rh (12 × 16), Ch (16 × 12), at
  frame index 75, per Q (`logs/openai_4p3/border_certificates.json`).
- The proof: Mc(z) − Mc(z*) has constant blocks and annihilates Rh^T, so the linking identity is unchanged.

**Claim 3, the elliptic seeds [M].** The definition is in Proposal 2 of the answer:
- **Positive controls:** cores and verified matrices at orders 60 (p = 3, both Q), 92 (p = 5, both Q) and 124 (p = 7,
  ℤ₂²).
- **Seed gate:** exhaustively empty at p = 11, 17 and 19. At p = 13 the seed gate holds, but the phase join failed
  for the one seed quadruple tried.
- **Targets:** the 16 targets pass the compression condition (C) only.

## 2. Tasks, in order

**Task 1. Claim 1: proof and exhaustive check (≤ 1.5 h).**
- Read the proof in the answer's Proposal 1 §2 line by line, for both Q, and write any gap and its fix.
- Implement your own brute force and run it to the largest odd m that fits in about 20 minutes per Q (aim for m ≤ 99).
  Report the counts per m.
- Cross-check against the repository's Q-level solution lists for a few census structures with odd |K|, e.g. the 863
  S7 host (`theory/bt_gates.py`, `Gates.q_level`).
- Reproduce the answer's minimum-change table for the 16 targets with your own code. Mismatch kills the table until
  resolved.

**Task 2. Claim 2: certificates, then real cores of other orders (≤ 2 h). The main deliverable.**
- **Certificates.** Verify both certificates in `logs/openai_4p3/border_certificates.json` with independent integer
  code:
  - all entries ±1;
  - X Xᵀ = 12I, Rh Rhᵀ = 16I, Ch Chᵀ = I₄ ⊗ Γ_Q (Γ_Q has diagonal 12 and off-diagonal −4);
  - X Chᵀ + Rh Mc(z*)ᵀ = 0.

  Build Mc from the definition of the Goethals–Seidel array over Q.
- **Real cores.** Assemble and verify Hadamard matrices from cores the OpenAI model did not use, with **only** the
  universal certificate:
  - **1916** (v = 479; G = ℤ₂² × ℤ₁₁₉, |K| = 119 odd; core in `data/core_v479_b3.txt`, loader in
    `theory/bt_controls.py` `c1916`);
  - **3964** (v = 991; G = ℤ₂² × ℤ₂₄₇, |K| = 247 odd).
    - Kumar's certificate is `hadamard/3964/certificate.json` in `https://github.com/vatskishan03/hadamard`
      (Apache-2.0).
    - To get it: call the `add_repo` tool for `vatskishan03/hadamard` with read access, then run
      `GIT_LFS_SKIP_SMUDGE=1 git clone --depth 50 https://github.com/vatskishan03/hadamard /home/user/vatskishan03/hadamard`.
    - Loader: `theory/bt_controls.py` `c3964`.
  - one random valid core from the elliptic controls at p = 5 and at p = 7.
- **How.** For each core:
  1. check the core by integer PAF on all lags;
  2. normalise it as in Claim 2 (record the row permutation, signs and Q-translations);
  3. assemble the bordered GS matrix with the certificate lifted constantly along K;
  4. verify H Hᵀ = 4n·I exactly, with a verifier that shares no code with the assembly.

  Also run a negative control: corrupt one sign and require rejection.
- **Pass:** 1916 and 3964 both verify. Then the border problem of the k4 family with odd |K| is closed.
  - Add one line to `README.md` crediting the OpenAI model.
  - Add a note to `report/b3_push_answer.md` §7 pointing to your review. Don't rewrite that section.
- **Kill:** any real core fails after correct normalisation. Find out whether the certificate, the normalisation or
  the proof is at fault, and say which.

**Task 3. Claim 3: reproduce, then decide (≤ 2.5 h).**
- **Reproduce independently:**
  - the five small cores and matrices (orders 60, 60, 92, 92, 124);
  - the exhaustive empty seed gates at p = 11, 17 and 19, with your own seed enumeration and join, sharing no code
    with `theory/openai_4p3_search.py`.

  A positive at 11, 17 or 19 refutes the answer's exclusion and needs an independently checked counterexample.
- **Run the answer's next experiment:**
  `python3 theory/openai_4p3_search.py --primes 23,29,31 --phase-max-p 0 --seconds 540 --out logs/openai_4p3/next_seeds.json`.
  Confirm any zero with your own implementation.
- **Estimate the scaling [E].** For each p in {11, 13, …, 31} and at the 16 target primes, compare:
  - log₂(number of distinct seed quadruples): λ, ε, u, the root masks and the switched-root choices per row;
  - the bits of the (p − 1)/2 seed-gate equations (K): lattice Gaussian `theory/q523_gauss.gauss_bits` on sampled
    seeds, read at the true target, calibrated against the exact counts at p ≤ 31.

  Give the expected number of seed-gate solutions at each target prime.
- **Decision.**
  - If the gate is exactly empty at p = 23–31 and the calibrated expectation at every target prime is below 2^−10,
    **close the ansatz** [M + E].
  - Otherwise, name the target primes worth an exact seed gate, with its cost.

**Task 4. Write `report/openai_4p3_review.md` (≤ 0.5 h)** with:
- the bottom line: verdicts on claims 1, 1′, 2 and 3, and what changes in the repository (the border theorem, if
  verified);
- per task: methods, scripts (`theory/4p3_review_*.py`), logs and results with tags;
- the next step for the n = 4p + 3 orders, given these results.

Commit, push, and stop.

## 3. Calibration warnings (from this repository)

- Zero counts are common below about 2^30. At the open orders, every exactly decided positive estimate has so far
  been zero, up to 2^85. So a positive expectation alone never justifies a long run.
- A selected-seed failure in a join is not a zero for the class. Equal PAF vectors need not have equal
  cross-correlations.
- A planted or small-order control checks correctness, not existence at the targets.
