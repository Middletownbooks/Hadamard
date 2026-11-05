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

[P, M] The broader **three-stratum product class** (arbitrary signs on all nine
rectangles, not necessarily a character product) cannot give an infinite family:
its two factor orders satisfy min(p,r)<=35 and max(p,r)<=16 min(p,r)-9,
so m<=19285. Exact necessary screens, including compression-conditioned subgroup
energies, rule out this one- or two-factor class at **38/38** targets. This is not
a nonexistence claim for unrestricted cores.

[P] Two additional direct approaches fail for every m>1: a character-independent
scalar norm product of two cores, and taking the four negative sets to be relative
transversals of K (even allowing complements). Confidence is proof-level for the
stated classes, with independent exact checks described below. No claim relies on [E].

## Analysis of the known cores

[M] `theory/openai_identity_known.py`, output `logs/openai_identity/known_cores.json`,
checks both full cores by integer PAF, converts each cell to CRT coordinates, and
examines **every** independent cell translation (119 or 247 centers).

| property | 1916 | 3964 |
|---|---|---|
| K | C7 x C17 | C13 x C19 |
| verified coupled symmetry order | 12 | 27 |
| order of its projection to the unit group of K | 12 | 27 |
| index of that projection in the unit group | **8** | **8** |
| size of the kernel acting trivially on Q | 4 | 9 |
| rowwise Q-translations making the supplied generators invariant | 0,3,1,0 | 0,0,0,0 |
| cells with some balanced rank-one punctured interior after translation | 0 of 16 | **1 of 16** |
| minimum flips to an arbitrary quadratic 3x3 sign table, per cell | 24–42 | 3–98 |

[M] These are orders of the **supplied verified subgroups**, not assertions about
the full automorphism groups. The index-eight coincidence is a structural observation,
not a construction theorem. Independent quadratic-class invariance is absent in
every cell, even after translation. The 3964 cell (row 0,q=0) does have a pure
quadratic-product interior, but its axes need three changes to become a full
quadratic-class table. Other cells are far from that model.

[M] The original Q-compressions are:

```
1916: (3,1,1,1), (1,1,1,3), (9,11,9,9), (5,5,5,5)
3964: (-11,-13,-13,-13), (-5,-7,-7,-7), (1,-1,-1,-1), (7,7,7,7)
```

[M] In 3964 the product-interior cell has primitive-character energy E=259=m+12,
while the cell (row 3,q=0) has E=86=m-161. The 16 energies sum to 3964=16m+12.
Theorem T would force each to be at least 248. The successful core uses low-energy
cells outside that tensor class to compensate for high-energy ones; filling axes
alone cannot do this. The log contains all 16 rational energies and all cellwise
distances for both cores.

[M] `theory/openai_identity_controls.py` verifies the six previously reviewed
cores at m=3,5,7 on both Q by the new independent integer PAF checker. All pass.
For m=3 it uses the reviewed Theorem B assembly to produce two order-60 matrices;
`theory/4p3_review_hverify.py`, in a separate process sharing no builder or assembly
code, verifies every row-pair dot product by integer popcounts. Both pass, and
corrupted cores and matrices are rejected. These are **existing positive controls**,
not matrices from the rejected proposals. Results: `positive_small.json` and
`border_controls.json` in `logs/openai_identity/`.

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
[paper](https://maths.nuigalway.ie/~padraig/Docs/OCathainStafford.pdf).
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
of the class**; the universal conclusion is the proof above. The observed result [M] is
296 rejected candidates and no energy below 16m+16. A binary candidate in the stated
class with even one cell E(f)<m+1 kills the theorem; a valid core kills it outright.
No parameter search over this class is recommended.

[M] The boundary inequality was also checked exhaustively, independently of the
candidate builder, in `theory/openai_identity_obstructions.py`:

| factor orders | boundary points | all binary boundary fillings | minimum E(B) |
|---|---:|---:|---:|
| 3 | 1 | 2 | 1 |
| 3 x 3 | 5 | 32 | 1 |
| 3 x 5 | 7 | 128 | 1 |
| 5 x 5 | 9 | 512 | 1 |
| 3 x 3 x 3 | 19 | 524288 | 1 |

These exhaustions verify the sharp lower bound on the boundary term, not a search
for cores. The proof of T supplies the conclusion for every size and filling.

## Proposal 2: arbitrary three-stratum product sign tables

### Definition

Let K=A x B, with odd p=|A|, r=|B|, both at least 3. In every Q-cell partition
each factor into one distinguished point and two equal halves of its complement.
Assign arbitrary signs to the resulting nine rectangles. Partitions and distinguished
points may vary independently with the cell. Quadratic residues/nonresidues in
two finite fields are a special case; the proof permits **any** balanced partitions.
This is broader than Proposal 1: the four interior rectangle signs need not be
separable or balanced. Single-factor three-stratum templates are included separately.

Reference builder: `theory/openai_identity_build.py --kind three-stratum --factors p,r`.
`--table-masks` takes 16 comma-separated integers in 0..511, in row-major (i,q)
order. Bit 3*a+b is the sign on stratum (a,b), with 0=distinguished point,
1=positive half, 2=negative half; bit 1 means +1. The default masks are fixed
diagnostic examples, not claimed core parameters. Whole-cell translations can
be applied to the emitted arrays; the obstruction is translation invariant.

### Subgroup fibre identity F [P]

For **any** subgroup N<=K of odd order d, summing the PAF target over N gives

    sum_{i, c in G/N} (sum_{x in c} x_i(x))^2
       = 16m - 12(d-1).                                    (F1)

Suppose each N-fibre has one distinguished point and two equally large sign-constant
classes. Its sum has absolute value 1, d-2, or d. Let A count the fibres of type d,
and B those of type d-2 (at d=3 the latter may be indistinguishable from type 1).
There are 16m/d fibres. Subtracting the baseline 1 per fibre in (F1), then dividing
by d-1, gives the exact necessary condition

    (d+1) A + (d-3) B = 16(m/d)-12,
    A,B >= 0, A+B <= 16m/d.                                (F2)

For a single factor (m=d), this is (d+1)A+(d-3)B=4, impossible for odd d>=9.
Unlike a probability estimate, F2 excludes every choice of signs and fibre centers.

### Theorem C [P]: no infinite two-factor three-stratum family

For two factors, an unbalanced p-fibre over a nonzero r-class occurs at least
(r-1)/2 times, because that class is a rectangle of the table. Equation F2 implies

    (p-3)(r-1)/2 <= 16r-12                                 (C1)

if such a fibre exists. When p,r>=37, the left side is at least 17(r-1), which
exceeds 16r-12. Thus every p-fibre over a nonzero r-point is balanced on its
nonzero p-points. Interchanging the factors gives the same property in the other
direction. The four interior signs are necessarily e times (1,-1,-1,1), namely
the balanced punctured tensor of Theorem T. The axis fillings are arbitrary
and T already excludes them. Hence min(p,r)<=35.

Also F2 for the larger factor d=max(p,r), with s=min(p,r), has a positive right
side and therefore at least one unbalanced fibre. It implies d-3<=16s-12, or
d<=16s-9. Consequently m<=35*551=19285. **There is no infinite family in this
class, even if all partitions, signs and centers vary with p,r.** If both factors
are finite fields, their orders are odd prime powers; the largest allowable smaller
order is 31, improving this bound to 31*487=15097. These are safe upper bounds,
not assertions of existence at their endpoints.

### Exact normal-form gate [P, M]

A cell sum is a signed sum of the nine rectangle sizes

    (1,(p-1)/2,(p-1)/2) tensor (1,(r-1)/2,(r-1)/2).

The script enumerates all 512 sign tables and obtains the set Z of possible sums.
It then enumerates a>0 and triples of odd t_j, up to permutation, satisfying

    a, t_j, t_j+2 in Z,
    a^2 + sum_j t_j(t_j+1) = pr.

It uses the reviewed compression normal form, so an empty list rules out the full
class on both Q. A nonempty list is **only a necessary-condition pass**. A second
gate is F2 for each factor.

[M] `theory/openai_identity_obstructions.py` checks all arithmetic two-factor
presentations of every target, including presentations that might not occur for
a particular abelian K; this only enlarges the class being screened. For prime m
it uses the single-factor F2. The separate gates exclude 35 targets. Their three
remaining cases, before the stronger gate below, are:

| order | m | surviving factor orders | possible (a;t1,t2,t3), up to t permutations |
|---|---:|---|---|
| 3004 | 187 | 11 x 17 | (1;-9,-9,-7), (5;-9,-1,9), (5;-7,5,9) |
| 6268 | 391 | 17 x 23 | (17;-7,5,5) |
| 9532 | 595 | 17 x 35 | (17;-1,-1,17) |

The 17 x 35 presentation uses an arbitrary balanced partition of a group of order
35; it is not a pair of finite fields. C says nothing about a full 27-rectangle
model on three factors (5,7,17). Similar higher-factor table classes must not be
confused with the balanced tensors already excluded by T.

### Compression-conditioned subgroup gate: all 38 targets excluded [P, M]

The separate gates discard which cell sums and fibre energies can occur together.
Retain that association. For each 3x3 sign table F put

    w=(1,(p-1)/2,(p-1)/2), v=(1,(r-1)/2,(r-1)/2),
    z = w F v^T,
    c_p = [sum_b v_b (sum_a w_a F_ab)^2 - r]/(p-1),
    c_r = [sum_a w_a (sum_b v_b F_ab)^2 - p]/(r-1).

These are nonnegative integers. For any core, its 16 cells must have total costs
(16r-12,16p-12). Enumerate all 512 tables, group their (c_p,c_r) pairs by z, and
for each normal form join the 16 required lists on that target. The dynamic program
keeps every attainable pair and discards only coordinates exceeding the nonnegative
target. It is an **exact necessary gate**, not an estimate or a restricted seed search.

[M] All profiles for the last three cases are empty at the target. Small explicit
certificates, from the table lists after each cell's nonnegative cost caps, are:

| order | normal-form profile | incompatible cost |
|---|---|---|
| 3004 | (1;-9,-9,-7) | total c_r has minimum 104; any increase is at least 140, but 164 is required |
| 3004 | (5;-9,-1,9) | total c_r has minimum 74; any increase is at least 140, but 164 is required |
| 3004 | (5;-7,5,9) | total c_r is at least 184, but 164 is required |
| 6268 | (17;-7,5,5) | total c_p is 268+308j, j>=0, never the required 356 |
| 9532 | (17;-1,-1,17) | total c_p is 140+476j, j>=0, never the required 548 |

[M] An independent implementation, `theory/openai_identity_gatecheck.py`, imports
neither the builder nor the screen. It expands the tables on the full groups and
checks the subgroup PAF sums by literal integer shifts; it uses an 8+8 join with
divide-and-conquer sum sets rather than the main screen's 16-step dynamic program.
It audits 7680 tables: ten small presentations, these three final cases, and both
known-core sizes. It confirms all five latter sizes have no table core. Its output
is `logs/openai_identity/conditioned_gate_audit.json`.

### Small checks and Claude test

[M] The same script records the exact gates for every odd p,r>=3 with p<=r and
pr<=99 (30 presentations), and single-factor prime-power sizes through 99. These
entries say `fails` or `unknown`, never `holds` merely because a gate passed.
This is a nonexistence result, not a proposed positive construction. The full
PAF checks in Proposal 1 test a common subclass at every odd m<=99; the valid
controls test the checker independently. In addition, the reference table builder
was run at all 30 small two-factor presentations, on both Q; all 60 diagnostic
candidates fail exact full PAF. That observation is **not** the nonexistence proof
for the broad class or for the small sizes at which the gates pass.

```
python theory/openai_identity_obstructions.py
python theory/openai_identity_build.py --kind three-stratum --suite --out /tmp/openai_identity_tables.json
python theory/openai_identity_verify.py /tmp/openai_identity_tables.json \
  --expect reject-candidate --out logs/openai_identity/three_stratum_small_paf.json
```

**Pass:** reproduce F2 from subgroup sums, reproduce the 512-table sum sets and
38/38 exclusions, and verify the C1-to-T implication. **Kill:** any genuine core
in this class with p,r>=37, or max(p,r)>16min(p,r)-9; or a table core at a target
marked `fails`. Small cases labelled `unknown` by the necessary gates are not
excluded; a positive example there would not contradict the general finite bound.
No long search is requested.

## Proposal 3: character-independent scalar norm products

### Definition and theorem N [P]

Suppose cores on Q x K1 and Q x K2, of odd subgroup sizes a,b, are combined into
four sequences on Q x K1 x K2 by a rule with the exact Fourier norm identity

    N_out(eta,chi1,chi2) = c N_a(eta,chi1) N_b(eta,chi2),     (N1)

where c is independent of the character (it may depend on a,b), and the same
quotient character eta is used in both factors. This includes a proposed norm
composition **only if it actually has N1**; no blanket assertion is made about
all quaternion or Turyn constructions.

For eta!=1, chi1=chi2=1, each input norm is 12 and the output must be 12, forcing
c=1/12. At the principal character the input norms are 16a+12 and 16b+12, whereas
the required output is 16ab+12. The excess is

    (16a+12)(16b+12)/12 - (16ab+12)
       = 16ab/3 + 16(a+b) > 0.                             (N2)

Thus no such rule works, even before imposing binary output coefficients. This
does not exclude character-dependent projections, auxiliary sequences, or a
product that changes how quotient characters are identified. Such a rule needs
a different norm identity and an explicit binaryness proof.

### Computation, coverage and Claude test

[M] `theory/openai_identity_obstructions.py` evaluates N2 as an exact rational
number on the ten unordered pairs with repetition from a,b in {1,3,5,7}.
All ten fail. This checks the algebraic norm obstruction, not PAFs of nonexistent
output sequences. No binary reference builder can satisfy N1, by N2; the norm
reference implementation is `norm_product(a,b)` in that script.

[P] Coverage: **fails for the stated product rule at all 38 targets**, for every
factorization. At prime subgroup orders there is no nontrivial factorization,
and allowing a factor of size 1 still gives the contradiction.

**Pass:** the special quotient norm forces c=1/12 and the principal excess equals
N2. **Kill:** a valid output core satisfying N1. A different product identity is
outside this theorem and remains open.

## Proposal 4: genuine relative-difference-set blocks

### Definition, dictionary and theorem RDS [P]

Allow each x_i, up to negation, to be G-2D_i, where D_i has no nonzero differences
in K. Equivalently D_i meets every K-coset at most once. Ordinary relative
difference sets with forbidden subgroup K are a special case; uniformity of their
differences outside K is not assumed. Complements are covered by the allowed
negation of the whole sequence. Write v=4m, k_i=|D_i| and use group sums G,K.

Expanding the exact PAF identity gives the useful general difference-family dictionary

    sum_i D_i D_i^(-1) = (v+3)e - 4K + lambda G,
    lambda = sum_i k_i - v + 1.                            (R1)

This dictionary holds for every binary core, not only the transversal class.
If m>1, a nonzero element of K exists; its left coefficient is zero in the
transversal class, forcing lambda=4. Therefore sum_i k_i=v+3=4m+3.

But a transversal cell has compression m or m-2. The reviewed normal form has
exactly one constant row and three singleton rows. A constant row has 0 or 4
occupied K-cosets; a singleton row has 1 or 3 occupied K-cosets. Hence

    sum_i k_i <= 4 + 3*3 = 13 < 4m+3  (odd m>=3),          (R2)

a contradiction. Row negations used to reach normal form do not alter the number
of equal entries, so the singleton/constant classification applies to the original
transversal rows too. It follows that four genuine relative-transversal blocks
cannot be the core at any odd m>1.

### Computation, coverage and Claude test

[M] `transversal_q_gate` in `theory/openai_identity_obstructions.py` enumerates
all 16 occupancy vectors and all 3876 four-row multisets, for every odd m=3,...,99
and both Q. All 98 cases have zero solutions to the necessary quotient PAF equations.
This is exhaustive for the quotient of the stated class; it is not a full-core
search. Arbitrary signs within occupied cosets do not arise for indicator sets,
and their point locations cannot change these quotient vectors.

[M] `theory/openai_identity_build.py --kind transversal --suite` also builds an
explicit diagnostic quadruple at every odd m<=99 on each Q: one negative set of
size 4 and three of size 3, each a transversal. All 98 examples fail full integer
PAF under the independent checker (`transversal_small_paf.json`). This is a
reference builder for the rejected proposal, not a positive family.

[P] Coverage: **fails at all 38 targets**. This does not rule out relative
difference sets used as auxiliary pieces of a more elaborate construction, or
supplementary divisible difference families with nonzero K-differences in R1.

**Pass:** verify R1 and the size contradiction R2; reproduce 98 empty quotient
enumerations. **Kill:** a core whose negative sets, after optional row negations,
all meet each K-coset at most once. The enumeration is the reference implementation;
there is no successful binary builder for this class.

## Coverage of all 38 supplied targets

[P, M] `fails` means the **defined construction class** cannot supply a core, not that the Hadamard order is impossible. T, N and RDS are universal proofs; C combines the proof and the exact necessary gates in `theory/openai_identity_obstructions.py`. Both Q are covered. Arithmetic two-factor presentations over-enumerate the possibilities for noncyclic groups as well. No unrestricted target is solved.

| order | n | m | T: balanced tensor | C: one/two-factor table | N: scalar norm | RDS: transversal |
|---|---:|---:|---|---|---|---|
| 2396 | 599 | 149 | fails | fails | fails | fails |
| 2876 | 719 | 179 | fails | fails | fails | fails |
| 3004 | 751 | 187 | fails | fails | fails | fails |
| 3452 | 863 | 215 | fails | fails | fails | fails |
| 4124 | 1031 | 257 | fails | fails | fails | fails |
| 4156 | 1039 | 259 | fails | fails | fails | fails |
| 4348 | 1087 | 271 | fails | fails | fails | fails |
| 4412 | 1103 | 275 | fails | fails | fails | fails |
| 4892 | 1223 | 305 | fails | fails | fails | fails |
| 5276 | 1319 | 329 | fails | fails | fails | fails |
| 5308 | 1327 | 331 | fails | fails | fails | fails |
| 5468 | 1367 | 341 | fails | fails | fails | fails |
| 5692 | 1423 | 355 | fails | fails | fails | fails |
| 5756 | 1439 | 359 | fails | fails | fails | fails |
| 5884 | 1471 | 367 | fails | fails | fails | fails |
| 5948 | 1487 | 371 | fails | fails | fails | fails |
| 6172 | 1543 | 385 | fails | fails | fails | fails |
| 6236 | 1559 | 389 | fails | fails | fails | fails |
| 6268 | 1567 | 391 | fails | fails | fails | fails |
| 6332 | 1583 | 395 | fails | fails | fails | fails |
| 6652 | 1663 | 415 | fails | fails | fails | fails |
| 6876 | 1719 | 429 | fails | fails | fails | fails |
| 7132 | 1783 | 445 | fails | fails | fails | fails |
| 7324 | 1831 | 457 | fails | fails | fails | fails |
| 7388 | 1847 | 461 | fails | fails | fails | fails |
| 7484 | 1871 | 467 | fails | fails | fails | fails |
| 7516 | 1879 | 469 | fails | fails | fails | fails |
| 8156 | 2039 | 509 | fails | fails | fails | fails |
| 8252 | 2063 | 515 | fails | fails | fails | fails |
| 8348 | 2087 | 521 | fails | fails | fails | fails |
| 8380 | 2095 | 523 | fails | fails | fails | fails |
| 8828 | 2207 | 551 | fails | fails | fails | fails |
| 8860 | 2215 | 553 | fails | fails | fails | fails |
| 9148 | 2287 | 571 | fails | fails | fails | fails |
| 9308 | 2327 | 581 | fails | fails | fails | fails |
| 9532 | 2383 | 595 | fails | fails | fails | fails |
| 9596 | 2399 | 599 | fails | fails | fails | fails |
| 9692 | 2423 | 605 | fails | fails | fails | fails |

## Ideas considered and rejected

- [P] Tensor Paley products with arbitrary coordinate-hyperplane fillings: T gives a four-unit spectral excess.
- [P, M] One- or two-factor three-stratum tables: C gives a finite bound; the exact gates exclude all 38 targets.
- [P] Character-independent scalar norm multiplication on a shared Q: special and principal spectra require incompatible scalings.
- [P] Four genuine relative-transversal negative sets: the difference-family and compression size counts conflict.
- [M] Explaining the two known cores as translated independent quadratic sign tables: every cell requires changes; the 512-table gates exclude those two sizes as well.
- Previously reviewed small cyclotomic-repair and elliptic search families: remain closed as specified in the prompt; not rerun here.

## What a Claude instance should run first

Run from the repository root with Python 3 and NumPy (the existing independent
matrix verifier uses NumPy's `bitwise_count`). Each invocation is a short exact
check, not a search campaign. The output JSON files are committed.

```bash
# 1. Audit the crucial finite gate with independently expanded subgroup PAFs.
python theory/openai_identity_gatecheck.py

# 2. Reproduce the 38-target coverage, boundary exhaustions and elementary obstructions.
python theory/openai_identity_obstructions.py

# 3. Reject the tensor candidates at every odd m through 99, both Q.
python theory/openai_identity_build.py --suite --out /tmp/openai_identity_tensor_suite.json
python theory/openai_identity_verify.py /tmp/openai_identity_tensor_suite.json \
  --expect reject-tensor --out logs/openai_identity/tensor_small.json

# 4. Check all 60 small diagnostic three-stratum candidates by full integer PAF.
python theory/openai_identity_build.py --kind three-stratum --suite --out /tmp/openai_identity_tables.json
python theory/openai_identity_verify.py /tmp/openai_identity_tables.json \
  --expect reject-candidate --out logs/openai_identity/three_stratum_small_paf.json

# 5. Positive controls: six known cores, then two independent order-60 matrix checks.
python theory/openai_identity_controls.py

# Explicit relative-transversal diagnostic candidates (98 full PAF checks).
python theory/openai_identity_build.py --kind transversal --suite --out /tmp/openai_identity_transversals.json
python theory/openai_identity_verify.py /tmp/openai_identity_transversals.json \
  --expect reject-candidate --out logs/openai_identity/transversal_small_paf.json

# 6. Optional reproduction of the real-core structural analysis.
curl -L --fail https://raw.githubusercontent.com/vatskishan03/hadamard/main/hadamard/3964/certificate.json \
  -o /tmp/openai_identity_3964.json
python theory/openai_identity_known.py --kumar /tmp/openai_identity_3964.json
```

[M] Expected counts: 15 independent table-audit cases (7680 tables), 38 excluded
targets for T/C/N/RDS, 296 rejected tensor examples, 60 rejected diagnostic tables,
98 rejected explicit transversal examples,
6 passing old cores, and 2 passing independently verified order-60 matrices.
The matrix checks include corrupted negative controls. The source 3964 certificate
used here has Git blob SHA `414fec749756d5cf597cb22a161cc16ff3be9314`; a later file
should be checked against that provenance before comparing detailed fingerprints.

Before attempting any new identity, check its cells against T and C and its norm
against N. A promising escape must change nonzero interiors, use more than two
independent sign-table factors, use genuinely coupled character classes, or employ
a norm identity with auxiliary/projected pieces. These are scope boundaries, not
positive proposals. Do not spend time searching the classes proved empty here.

## Unfinished items

- No positive algebraic identity producing a new core has been found.
- General coupled cyclotomy (including the index-eight symmetry seen in the controls),
  higher-order characters, general three-or-more-factor sign tables, Galois-ring
  constructions, and product identities outside N1 remain unresolved.
- C's finite bound is not a classification of all smaller sizes: the committed
  small-size gates deliberately retain `unknown` where necessary conditions pass.
- The mathematical proofs and new scripts await the requested external Claude review;
  they have not been labelled as already reviewed by Claude.

## Provenance and setup

The bordered construction is Theorem 7 of Egan–Kotsireas–London–Ó Catháin, in Kumar's form. Theorem B and the compression normal form are the prior OpenAI results reviewed in `report/openai_4p3_review.md`.

[M] The two requested branches were fetched and merged cleanly. Direct Git pushing
lacked credentials, so the commits were published through the GitHub connector.
The initial published checkpoint `9ec96bea` is a two-parent merge of the prior
Codex head `27cfc29c` and Claude head `9ae002e8`, with the answer skeleton added;
both histories are preserved. All updates are descendants of that commit and
target only `codex/order-3452-compression`. No pull request was opened.
