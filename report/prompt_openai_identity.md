# Prompt (OpenAI model; ideas and proofs, short sessions): an algebraic identity for the k4 bordered family with odd |K| (38 open orders)

*Self-contained, for a session with access to the public repository `https://github.com/middletownbooks/hadamard`.*
- **Push target:** you may push only to the branch `codex/order-3452-compression`.
- **The task:** find an **algebraic identity**, a Paley-type construction, that produces the four sequences below
  for an infinite or large family of groups. Or prove that natural classes of such identities cannot exist.
- **Checking:** Claude Code instances will verify proposals exactly afterwards. Write them so that a Claude instance
  can build and test each one without asking you anything.

## 0. Working rules (sessions may time out: checkpoint early and often)

- **Setup:**
  - Run `git fetch origin codex/order-3452-compression claude/focused-planck-j1645f`.
  - Check out `codex/order-3452-compression` and merge `origin/claude/focused-planck-j1645f` into it.
  - If the merge is not clean, abort it and read the Claude branch with `git show origin/claude/focused-planck-j1645f:PATH`.
- **Pushing:** push only to `codex/order-3452-compression`. Never push to `main` or to `claude/focused-planck-j1645f`.
  Open no pull requests.
- **Checkpoints:**
  - In your first ten minutes, create `report/openai_identity_answer.md` with the skeleton of §6, and push it.
  - Push again after every finished proposal, or every 30 minutes, whichever comes first.
  - If the session ends early, the file says what is unfinished.
- **Runs:** keep each computation under about 10 minutes. Put scripts in `theory/openai_identity_*.py` and logs in
  `logs/openai_identity/`.
- **Tags:** mark every claim [P] proved, [M] exact computation (script named) or [E] estimate.
- **Matrices:** never claim a Hadamard matrix without an exact integer check by code that does not share the
  construction's code.
- **Credit:** the bordered construction is Theorem 7 of Egan–Kotsireas–London–Ó Catháin [EKLO], in Kumar's form.
  Theorem B and the normal form below are yours (verified by Claude Code in `report/openai_4p3_review.md`).

## 1. The problem (exact) and what is now settled

**Data.**
- G is abelian of order 4m with **m odd**, K ≤ G of order m, and Q = G/K ∈ {ℤ₂², ℤ₄}. So G = Q′ × K.
- Wanted: x₁, …, x₄ : G → {±1} with Σᵢ PAFᵢ = 4|G| at 0, −12 on K \ 0 and +4 off K.
  - Equivalently Σᵢ |x̂ᵢ(χ)|² = 4n for χ nontrivial on K, and 12 for the three χ nontrivial on Q and trivial on K,
    with n = |G| + 3.

**Settled, and usable freely:**
1. **Theorem B [P, M]:** every such core completes to a Hadamard matrix of order 4n with one fixed border certificate
   per Q (`logs/openai_4p3/border_certificates.json`).
   - Verified on real cores of orders 1916 and 3964 (`report/openai_4p3_review.md` §2).
   - Assembly and verification take seconds (`theory/4p3_review_border.py`, `theory/4p3_review_hverify.py`).
   - **A core is all that is needed.**
2. **Normal form [P, M]:** the Q-compressions zᵢ(q) = Σ_{k∈K} xᵢ(q + k) are, up to row operations:
   - one constant row (a, a, a, a);
   - three rows (tⱼ + 2, tⱼ, tⱼ, tⱼ);
   - with m = a² + Σⱼ tⱼ(tⱼ + 1).
3. **What fails [P, M + E]:**
   - **Search over small families.** Families with polynomially many seeds and shifts per row fail. Each of the
     about 2m remaining correlation equations costs about 5–6 bits, against O(log m) bits of freedom, so they fall
     below one expected solution by m ≈ 20 (`report/openai_4p3_review.md` §3.3–3.4, §4).
   - **Your two earlier families** (fixed-index cyclotomic seeds with small repairs, and elliptic seeds) are closed.
   - **Multiplier slices** with up to 160 orbits are closed at m = 149 and 179.
   - Obstructions: Theorem α and Lemma R (`report/b3_push_answer.md` §3).

**So the requirement is an identity.** A viable proposal must make the correlation sums constant **by algebra**, as
the quadratic character makes Paley's PAFs constant, leaving at most O(log m) conditions to search or check.

## 2. The 38 open orders it would cover (n ≡ 7 mod 8, so m = (n − 3)/4 is odd)

| 4n | n | m = \|K\| |
|---|---|---|
| 2396 | 599 | 149 |
| 2876 | 719 | 179 |
| 3004 | 751 | 11·17 |
| 3452 | 863 | 5·43 |
| 4124 | 1031 | 257 |
| 4156 | 1039 | 7·37 |
| 4348 | 1087 | 271 |
| 4412 | 1103 | 5²·11 |
| 4892 | 1223 | 5·61 |
| 5276 | 1319 | 7·47 |
| 5308 | 1327 | 331 |
| 5468 | 1367 | 11·31 |
| 5692 | 1423 | 5·71 |
| 5756 | 1439 | 359 |
| 5884 | 1471 | 367 |
| 5948 | 1487 | 7·53 |
| 6172 | 1543 | 5·7·11 |
| 6236 | 1559 | 389 |
| 6268 | 1567 | 17·23 |
| 6332 | 1583 | 5·79 |
| 6652 | 1663 | 5·83 |
| 6876 | 1719 | 3·11·13 |
| 7132 | 1783 | 5·89 |
| 7324 | 1831 | 457 |
| 7388 | 1847 | 461 |
| 7484 | 1871 | 467 |
| 7516 | 1879 | 7·67 |
| 8156 | 2039 | 509 |
| 8252 | 2063 | 5·103 |
| 8348 | 2087 | 521 |
| 8380 | 2095 | 523 |
| 8828 | 2207 | 19·29 |
| 8860 | 2215 | 7·79 |
| 9148 | 2287 | 571 |
| 9308 | 2327 | 7·83 |
| 9532 | 2383 | 5·7·17 |
| 9596 | 2399 | 599 |
| 9692 | 2423 | 5·11² |

K may be any abelian group of order m: cyclic, ℤ₅² × ℤ₁₁ at 1103, ℤ₁₁² × ℤ₅ at 2423, and so on.

## 3. Known cores of this family (data for pattern-finding)

- **1916** (n = 479, m = 119 = 7·17, Q = ℤ₂²): the published core of Alpöge, Voinov and Reynolds-Haertle,
  `data/core_v479_b3.txt`. It is invariant under an order-12 group after a translation; see `theory/bt_controls.py`
  `c1916` and `report/b3_push_answer.md` §1.
- **3964** (n = 991, m = 247 = 13·19): Kumar's certificate, `hadamard/3964/certificate.json` at
  `https://github.com/vatskishan03/hadamard`. Its S has order 27; see `theory/bt_controls.py` `c3964`.
- **Small m:** `theory/bt_brute_count.py` counts cores exactly for small structures, and
  `logs/bt/calib/calib_small.jsonl` holds exact counts for 561 small structures (n ≤ 31). The elliptic controls give
  cores at m = 3, 5 and 7 (`logs/openai_4p3/`).

Both 119 = 7·17 and 247 = 13·19 are products of two primes. The search-found cores may hide an algebraic pattern,
for example cyclotomy of the ring ℤ_{pq}, or a product of two prime structures. Studying them is a legitimate
starting point.

## 4. Directions worth considering (not exhaustive)

1. **Two-prime cyclotomy.** For m = p·r, use the unit group of ℤ_m and generalised cyclotomic classes, in the style
   of twin-prime difference sets, combined over Q.
   - Many of the 38 orders have m = p·r: 751, 863, 1039, 1223, 1319, 1367, 1423, 1487, 1567, 1583, 1663, 1783, 1879,
     2063, 2207, 2215 and 2327.
   - Test the idea on 1916 and 3964 first: are their cores of this form?
2. **Finite fields and Galois rings:**
   - K = F_q additive, where m is a prime power (prime m here);
   - or K a Galois-ring quotient (Kumar's 839 used ℤ₄³ × ℤ₁₃);
   - with Q realised inside a ring extension, so that the character sums become Gauss or Jacobi sums with known
     absolute values.
3. **Relative difference sets** (or semi-regular ones) relative to K in G, and the Fourier dictionary between them
   and these PAF conditions.
4. **Product or recursive identities:** cores at m₁ and m₂ giving a core at m₁m₂, or at another combination. A
   product theorem would turn the small cores (m = 3, 5, 7, 119, 247, …) into many orders.
5. **Obstructions:** a proof that a natural identity class cannot work (for example, every Gauss-sum construction
   over F_m with Q = ℤ₂²) is a useful result.

## 5. Proof obligations and screens before any claim

- **Small-m test:** every proposal must be built and verified by exact integer PAF at **all** small m in its class,
  up to about m ≤ 99, or the first ten instances.
  - Record holds / fails per m.
  - A construction that fails at small m is not a family.
- **Bordered matrix:** for at least one small m, assemble the matrix with Theorem B and verify H Hᵀ = 4n·I with code
  sharing nothing with the construction.
- **Coverage:** for each of the 38 orders of §2, say whether the construction's conditions hold (holds / fails /
  unknown), with the computation.

## 6. The answer file: `report/openai_identity_answer.md`

It has these sections:
- **Bottom line:** the best proposal, its coverage of the 38 orders, and confidence.
- **Analysis of the known cores** (1916, 3964, small m): any algebraic structure found, or its absence.
- **One section per proposal:**
  - the definition, with a reference builder `theory/openai_identity_build.py`;
  - the proof, or the remaining condition;
  - the small-m results;
  - the coverage;
  - a test for Claude with pass and kill criteria fixed in advance.
- **Ideas considered and rejected,** one line each.
- **What a Claude instance should run first.**
- **Unfinished items.**

Commit with a message naming the proposals, and push to `codex/order-3452-compression`.
