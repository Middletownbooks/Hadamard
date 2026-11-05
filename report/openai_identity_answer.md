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
so m<=19285. Exact necessary screens rule out this class at **35/38** targets.
Only 3004, 6268 and 9532 remain undecided **within that class**. This is not a
claim that those three admit cores, nor that the other orders have no cores.

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
it uses the single-factor F2. Exactly 35 targets are excluded. The survivors are:

| order | m | surviving factor orders | possible (a;t1,t2,t3), up to t permutations |
|---|---:|---|---|
| 3004 | 187 | 11 x 17 | (1;-9,-9,-7), (5;-9,-1,9), (5;-7,5,9) |
| 6268 | 391 | 17 x 23 | (17;-7,5,5) |
| 9532 | 595 | 17 x 35 | (17;-1,-1,17) |

The 17 x 35 presentation uses an arbitrary balanced partition of a group of order
35; it is not a pair of finite fields. C says nothing about a full 27-rectangle
model on three factors (5,7,17). Similar higher-factor table classes must not be
confused with the balanced tensors already excluded by T.

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
35/38 exclusions, and verify the C1-to-T implication. **Kill:** any genuine core
in this class with p,r>=37, or max(p,r)>16min(p,r)-9; or a table core at a target
marked `fails`. A valid core at one of the three `unknown` targets would not
contradict the theorem. No long search is requested.

## Proposal 3: norm products and genuine relative-difference-set blocks

Unfinished: record the two short obstructions, their exact checks, and their scope.

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
