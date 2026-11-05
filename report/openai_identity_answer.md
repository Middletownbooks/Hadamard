# Algebraic identities for the odd-|K| k4 bordered family

Session: 2026-11-05. Authorized push branch: `codex/order-3452-compression`.
Claim tags: [P] proof, [M] exact computation with named script, [E] estimate.

## Bottom line

[P] **A natural tensor-Paley class is impossible at every odd m > 1, on both Q.**
Theorem T below permits arbitrary binary fillings on all coordinate hyperplanes,
arbitrary balanced punctured factors, and independent translations of every Q-cell.
It therefore rules out direct twin-prime-character products with arbitrary axis repairs.
Its obstruction is an exact four-unit excess in a necessary spectral average.
It excludes this construction class at **38/38 targets**, not the target cores themselves.
No new family or new Hadamard matrix is claimed.

## Analysis of the known cores

Unfinished: inspect the 1916 and 3964 cores for two-prime multiplicative structure; compare their partial compressions with small examples.

## Proposal 1: balanced punctured tensors with arbitrary boundary filling

### Definition and reference builder

Write K = K_1 x ... x K_s, with odd d_j = |K_j| > 1, m = product d_j,
and L = product (d_j - 1). In each of the 16 cells f = x_i(q, .), independently,
choose h_j : K_j -> {0,+1,-1} satisfying h_j(0)=0, |h_j(t)|=1 for t != 0,
and sum h_j=0. Require

    f(k) = epsilon product_j h_j(k_j), when every k_j != 0.

On the union of the coordinate hyperplanes choose **any** values in {+1,-1};
there is no bound on their number or any required symmetry. Independent choices
of factors, signs, and distinguished points in every cell are allowed. Translations
move the distinguished points to zero without changing the spectral averages used below.
The decomposition of K is shared; its factors can be arbitrary abelian groups.

`theory/openai_identity_build.py` builds examples of this *rejected* class. For
prime factors it uses the quadratic character; for other odd cyclic factors it
uses a balanced punctured sign vector. `--filling constant` attains the sharp
lower bound; `--filling varied` supplies another boundary choice. The builder
does not claim to produce valid cores.

### Theorem T [P]: no member of this class is a core

Let P_j = I - J/d_j be the orthogonal projection onto mean-zero functions on K_j,
and P = tensor_j P_j. These are precisely the characters nonprincipal on **every**
factor. For a cell f define its average energy on those L characters by

    E(f) = (1/L) sum_{chi_j != 1 for every j} |fhat(chi)|^2
         = (m/L) ||P f||^2.                                  (T1)

The core equations, summed over the four quotient characters and using Parseval
on Q, require

    sum_{i,q} E(x_i(q,.)) = 16m + 12.                         (T2)

This uses only characters nontrivial on K, so it holds for either Q.

Put H = epsilon tensor h_j and B=f-H. Then PH=H, ||H||^2=L,
and H and B have disjoint supports. Thus <H,PB>=<H,B>=0 and

    E(f) = m + (m/L)||P B||^2.                               (T3)

We show that the second term is at least 1, **for every binary boundary filling**.
For each factor let D_j send a vector v to (v(t)-v(0))_{t != 0}, and let
D = tensor D_j. Each coordinate of DB, indexed by an interior point, is an
alternating sum over its 2^s cube corners. The all-nonzero corner contributes 0;
the other 2^s-1 entries are each +1 or -1. Hence every coordinate of DB is an
odd integer, and ||DB||^2 >= L.

Now D_j P_j = D_j and D_j D_j^T = I + J of size d_j-1, whose largest eigenvalue
is d_j. Consequently D=DP and ||D||^2 = product d_j=m. Therefore

    L <= ||DB||^2 = ||D P B||^2 <= m ||P B||^2,

proving E(f) >= m+1. The 16 cells give at least 16m+16, contradicting (T2).
This proves Theorem T. No estimates, character-sum conjectures, or search assumptions occur.

The bound is sharp for an individual cell: B=1 on the boundary and 0 in the
interior gives E(B)=1. It still exceeds the available **total** energy by four.

### Scope and coverage [P]

- The ordinary quadratic character of any odd finite field is an eligible h_j.
- No Gauss-sum magnitude formula is needed: arbitrary balanced punctured factors work.
- All prime m, all two-prime m, and all more general product decompositions are covered.
- Per-cell translations, factor automorphisms, whole-row operations, and either Q
  do not evade the proof.
- The theorem does **not** cover a nonseparable or unbalanced interior sign table,
  higher-order cyclotomy in general, or repairs that change the nonzero interior.
- All 38 listed orders: **fails for this class**. General existence remains unknown.

The classical twin-prime-power difference set uses a quadratic-character product
on nonzero pairs with prescribed axes: see O Cathain and Stafford, *On twin prime
power Hadamard matrices*, Definition 14,
https://maths.nuigalway.ie/~padraig/Docs/OCathainStafford.pdf.
Theorem T concerns placing such sign templates in this k4 core; it does not deny
the classical difference-set construction at its own parameters.

### Small-m checks and Claude pass/kill criteria

`theory/openai_identity_verify.py` imports no builder. It computes every PAF by
integer shifts and computes (T1) independently using tensor integer centering
matrices d_j I-J, keeping the energy as a rational number.

```
python theory/openai_identity_build.py --suite --out /tmp/openai_identity_tensor_suite.json
python theory/openai_identity_verify.py /tmp/openai_identity_tensor_suite.json \
  --expect reject-tensor --out logs/openai_identity/tensor_small.json
```

The suite covers every odd m=3,...,99, both Q, the cyclic single-factor presentation
and the prime-factor product presentation (deduplicated), and two boundary fillings.
These are checks of representative rejected candidates, **not exhaustive enumeration
of the class**; the universal conclusion is the proof above. The expected result is
296 rejected candidates and no energy below 16m+16. A binary candidate in the stated
class with even one cell E(f)<m+1 kills the theorem; a valid core kills it outright.
No parameter search over this class is recommended.

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
