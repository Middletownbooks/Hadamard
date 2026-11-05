# Algebraic proposals for the k4 bordered family at n = 4p + 3

Status: completed short-session report, 2026-11-05. Two ranked proposals, exact small controls,
and reproducible scripts are supplied. No Hadamard matrix at an open order is claimed.

Tags: [P] proof supplied; [M] exact computation with the script and scope named; [E] estimate with method named. Untested proposals are explicitly conditional.

Credit: the bordered construction is Theorem 7 of Egan–Kotsireas–London–Ó Catháin (EKLO), in Kumar's form (https://github.com/vatskishan03/hadamard). The census, pre-existing obstructions, calibration and border tooling are by Claude Code in this repository.

Setup: GitHub's comparison found codex/order-3452-compression at 5b26dae60414599ccccf4c692220f616fc2aae24 to be an ancestor of claude/focused-planck-j1645f at a87392b4d91fcaafe0bc61d85435b4ff0a06c878 (nine commits ahead). The target branch was fast-forwarded to that reference commit with an expected-head check. The first checkpoint was pushed while the shell environment was unavailable. A shell subsequently became available; a separate checkout of the target branch was cloned, both branches fetched, and the merge reported already up to date. Only codex/order-3452-compression is written. No PR is opened.

## Bottom line

**The strongest result is an obstruction, not a new open-order construction.**

1. **Rank 1, high confidence [P, M]:** an exact compression normal form for both quotient groups rules out
   small repairs to balanced character seeds. Across the 16 targets, at least **29–61 total sign changes**,
   or a common per-cell allowance of at least **4–7**, are necessary. These minima concern compression,
   not full cores. More generally, **fixed cyclotomic index and a fixed number of exceptional positions per
   cell cannot give an unbounded prime family**, even with independently translated cells and arbitrary
   class unions. The proof and finite arithmetic gate are below.
2. **Rank 2, conditional and low confidence as an infinite family [M]:** cubic character seeds with three
   root changes give explicit cores for p=3,5 on both groups and p=7 on Z_2^2. Independent integer code
   verifies their bordered matrices at orders **60,92,124**. The complete seed gate is empty at p=11,17,19,
   so this is **not** a construction for every prime. All 16 targets pass compression only; full existence
   is **unknown at every one**. No target-prime correlation search was performed.
3. **Supporting result [P, M]:** the normal form and a finite certificate establish a border for **every
   valid core of this family over either Q**, after row sign/order/translation equivalences. This extends
   the existing tested border classes to a proof covering all odd compressions; it uses the existing
   Claude Code border reduction and does not supply the missing cores.

All computations were deterministic and used integer arithmetic. The longest completed invocation here
was the p=13 selected-seed phase test, about 44 seconds. Scripts and exact commands are in the final section.

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


## Proposal 1 (rank 1): a compression normal form and an exact obstruction to small repairs

### 1. Definition and scope

A cell starts with a translated and dilated cyclotomic sign sequence of index e, e | p-1.
Choose a primitive root g, arbitrary signs s_0,...,s_(e-1), a zero sign eta, a shift h and a nonzero dilation u.
At k, use eta if u(k-h)=0; otherwise use s_j when u(k-h) belongs to g^j<g^e>.
Then toggle an arbitrary set F of at most c positions. All 16 cells may use different signs, shifts, dilations and repair sets.
This is substantially broader than common multiplier invariance: different shifts usually destroy it.
The reference implementation is `theory/openai_4p3_build.py --ansatz cyclotomic`.
Its JSON parameter object is `{"cells": [[...four cells...], ...four rows...]}`. Each cell has
`e`, `signs` (e signs), `zero` (+/-1), `shift`, `scale` (nonzero mod p), and `flips`
(distinct original coordinate positions). The primitive root is chosen as the least primitive root.
The theorem and gate below use one common e and common upper budget c; the builder also permits
different e in different cells, which is a broader class and requires its corresponding sum sets.

The particularly broad balanced version allows *any* initial binary cell of sum +1 or -1, including arbitrary
permutations of completed Legendre sequences, followed by arbitrary sign changes.

### 2. Proof

**Compression normal form [P].** For either Q, after permuting whole rows and negating whole rows, every odd
4-by-4 compression satisfying the Q equations consists of one constant row (a,a,a,a) and three rows
which have three entries t_j and one entry t_j+2, with a and all t_j odd. The singleton positions may differ.
The exact scalar condition is

    p = a^2 + t_1(t_1+1) + t_2(t_2+1) + t_3(t_3+1).       (NF)

Conversely these conditions imply all four Q-level equations.

For Q=Z_2^2, set w_i,eta=(1/2) sum_q eta(q)z_iq. Each nonprincipal column has squared norm 3, hence three
entries +/-1 and one 0. Since every z_iq is odd, the three nonprincipal w entries in a given row have the
same parity: for example, w_i,eta-w_i,theta is a difference of two odd z entries. Thus their zero occurs
in the same row. That row is constant. In every other row the three entries are +/-1; inverse Walsh
transform gives one exceptional entry differing by 2 from the other three. Negate that whole row if needed.

For Q=Z_4 put A_i=(z_i0-z_i2)/2, B_i=(z_i1-z_i3)/2 and
C_i=(z_i0-z_i1+z_i2-z_i3)/2. The nonprincipal energies give
sum_i(A_i^2+B_i^2)=3 and sum_i C_i^2=3. Also C_i = A_i-B_i (mod 2).
Each of the three rows with C_i odd therefore uses at least one unit of the first energy budget.
Consequently each has exactly one of A_i,B_i equal to +/-1, the other zero; the remaining row has A=B=C=0.
Solving these relations gives the same constant/singleton form. Finally, the total zero-lag compression
energy is 4a^2 + sum_j(4t_j^2+4t_j+4)=4p+12, giving (NF).
For a singleton row the off-zero PAF on either Q is 4t^2+4t, proving the converse.

**Sharp compression repair gate [P].** A balanced cell needs at least
d(z)=(|z|-1)/2 sign changes to reach odd sum z. This number is attainable at the sum level
when the initial sum's sign can be chosen independently. For a fixed initial sum s in {-1,1},
the exact cost is |z-s|/2; allowing either s gives the stated lower bound.
Hence the least *possible at the compression level* total repair cost is

    min over (NF) of 4d(a) + sum_j [3d(t_j)+d(t_j+2)].

The least common per-cell cap is the minimum of max d(z_iq) over those forms. These are lower bounds
for actual cores, not a claim that the PAF equations can be satisfied at these budgets.

For the general cyclotomic rule put s=(p-1)/e. Its unmodified cell sums are
B={eta+s(e-2j): eta in {-1,1}, 0<=j<=e}.
From a base b, flipping u minus signs up and v plus signs down gives z=b+2(u-v), where
u+v<=c, u<=(p-b)/2 and v<=(p+b)/2. This gives an exact finite set S(p,e,c) of attainable cell sums.
The ansatz passes compression iff (NF) has a solution with a,t_j,t_j+2 all in S.
Necessity follows from the normal form; sufficiency at this level follows by selecting each cell independently.
The set S is closed under negation, so the row sign normalization loses nothing.

**Fixed-budget consequence [P].** If e is even, put B0=2c+1.
When 2(p-1)/e-B0 > sqrt(4p+12), any attainable sum within the target energy bound has |z|<=B0.
Then no core exists if p>4B0^2-3. If e is odd, every base has |b|>=(p-1)/e-1.
If L=(p-1)/e-1-2c>0 and 16L^2>4p+12 there is no core.
Thus fixed index e and fixed repair budget c cannot give an unbounded prime family.
This covers cell-dependent affine twists and different class unions.

### 3. Exact checks

[M] The embedded HOST_JS in `theory/openai_4p3_compression.py` was executed using V8, with integer arithmetic
(all values are far below 2^53), no random seed. logs/openai_4p3/compression_host.json stores the complete
minima and witnesses for all primes from 5 through 97 and the 16 targets. After the runtime became available,
the pure-Python implementation reproduced every target minimum (`compression_python.json`).
`theory/openai_4p3_check.py` checks the witnesses by direct quotient correlations on both Q and independently
enumerates the cyclotomic builder's attainable sums with zero or one repair at p=5,7,13, for applicable
e in {2,4,6}. Results: `logs/openai_4p3/checks.json`.

[M] For each target, every applicable e in {2,4,6,8,12,16} has the same minimum per-cell budget shown
in the table. Thus these common-index classes with c<=3 fail at all 16 targets, including cell-dependent
translations. This statement is a compression computation, not a rerun of the multiplier census.

### 4. Coverage

At all 16 targets: **fails** for the balanced-base ansatz below the respective budgets in this table.
At or above a budget: **unknown** for full cores. The two minima optimize different objectives and need
not be attained by the same compression.

| 4n | p | minimum total changes | minimum common per-cell cap |
|---|---:|---:|---:|
| 2396 | 149 | 29 | 4 |
| 2876 | 179 | 29 | 4 |
| 4124 | 257 | 37 | 5 |
| 4348 | 271 | 38 | 5 |
| 5308 | 331 | 41 | 5 |
| 5756 | 359 | 43 | 5 |
| 5884 | 367 | 39 | 6 |
| 6236 | 389 | 41 | 6 |
| 7324 | 457 | 47 | 6 |
| 7388 | 461 | 47 | 6 |
| 7484 | 467 | 43 | 6 |
| 8156 | 509 | 44 | 7 |
| 8348 | 521 | 49 | 6 |
| 8380 | 523 | 51 | 6 |
| 9148 | 571 | 53 | 6 |
| 9596 | 599 | 61 | 7 |

### 5. Test for Claude

Run python3 theory/openai_4p3_compression.py --out logs/openai_4p3/compression_python.json.
Compare the two minimum columns with compression_host.json, and independently verify all Q-level witnesses.
A discrepancy kills the implementation or proof and must be resolved before using this screen.
A zero from the gate kills only the specified (p,e,c) class.
A positive gate is permission to test correlations, not evidence that a core exists.


## Proposal 2 (rank 2): elliptic-character seeds with three root changes

### 1. Explicit construction

Use rows i=0,1,2,3. Let chi(0)=0 and chi be the quadratic character of F_p. In row i choose
lambda_i in F_p minus {0,1}, epsilon_i in {-1,1}, a nonzero u_i, and three root signs
eta_i0,eta_i1,eta_i2. Define f_lambda(v)=v(v-1)(v-lambda) and

    b_i(k) = epsilon_i chi(f_lambda_i(u_i k)) when u_i k is not a root;
             eta_ij                         when u_i k is root j in (0,1,lambda_i).

For each i=1,2,3 choose a root position r_i with b_i(r_i)=-1, and set

    b'_i = b_i + 2 delta_(r_i).

Choose four shifts h_iq for each row. The 16 binary cells are exactly

    a_0q(k) = b_0(k-h_0q),
    a_i0(k) = b'_i(k-h_i0)        (i=1,2,3),
    a_iq(k) = b_i(k-h_iq)         (i=1,2,3; q=1,2,3).

This rule works as a *candidate generator* on both Q. All h_i0 may be set to zero, since independent
whole-row translations in K preserve PAFs. Dilation u_i is common across the four cells of a row;
allowing different dilations within a row is outside this proposal and outside its negative results.

Implementation: `theory/openai_4p3_build.py --ansatz elliptic`. Parameters are a list of four objects.
`lam`, `sign`, `scale` are lambda, epsilon, u. `mask` is an integer 0..7, with bit j=1 giving eta_ij=-1.
Rows 1..3 also have `root`, the original k-coordinate r_i; `shifts` is [h_i0,...,h_i3]. The builder validates
the chosen root and outputs four sequences in q*p+k order. A complete p=5 parameter file is
`logs/openai_4p3/params_p5_v4.json`. Built-in `--example` controls cover every success listed below.

### 2. Reduced conditions and the step still open

**Compression [P].** Let T_lambda=sum_k chi(f_lambda(k)); this integer is the negative elliptic trace
in the usual point-count convention, since the projective curve y^2=f_lambda(x) has
1+sum_x(1+chi(f_lambda(x))) = p+1+T_lambda points. No estimate for this trace is needed here.
Then c_i=sum b_i=epsilon_i T_lambda_i+eta_i0+eta_i1+eta_i2. The compression test is exactly

    c_0^2 + c_1(c_1+1) + c_2(c_2+1) + c_3(c_3+1) = p.       (C)

This bypasses the balanced-seed obstruction by permitting nonzero cubic character sums.

**Seed gate [P].** Put P_i(t)=PAF(b_i)(t) and P'_i(t)=PAF(b'_i)(t). The equations on K are independent
of every shift h_iq. They hold iff, for each t=1,...,(p-1)/2,

    P_0(t) + sum_(i=1..3) (3P_i(t)+P'_i(t))/4 = -3.         (K)

All numerators are divisible by 4: every binary PAF at length p is p mod 4 (the number of disagreements
around a permutation cycle is even). Equivalently, directly expanding the one-position change gives

    P'_i(t)=P_i(t)+2[b_i(r_i+t)+b_i(r_i-t)]+4 delta_0(t).

These identities follow by multiplication; hence (K) is necessary and sufficient for every K-lag,
including lag zero automatically. It also implies (C) by summing over t. The separate (C) test is cheaper.

For a character-sum implementation, write F_i(k)=epsilon_i chi(f_lambda_i(u_i k)) and let R_i be its
three root positions. Then exactly

    P_i(t) = sum_k chi(f_lambda_i(u_i k) f_lambda_i(u_i(k+t)))
             + sum_(r in R_i) eta_ir [F_i(r+t)+F_i(r-t)]
             + sum_(r,s in R_i) eta_ir eta_is 1_(s-r=t).

Thus the unresolved identities are explicit degree-six character sums with known root terms, not an
unspecified numerical objective. No formula making their sum constant for an unbounded prime class
has been found.

**Phase gate [P].** Define u_iq=b_i except u_i0=b'_i for i>0. For every d!=0 and t in F_p, the remaining
necessary and sufficient condition is

    sum_(i,q) PCF(u_iq,u_i,q+d; t+h_iq-h_i,q+d) = 4.        (S)

Here PCF(v,w;t)=sum_k v(k)w(k+t), with q+d in the chosen quotient. Substitution k'=k-h_iq proves
the formula. Conditions (K) and (S) therefore certify a core exactly. They are reduced tests,
**not a proof that suitable parameters exist**.

`theory/openai_4p3_search.py` enumerates all lambda, twist, root signs, nonzero dilation and negative
root choices for (K), then does a full-tuple 2+2 join of the integer vectors. A zero is rigorous for
this entire seed class. For a positive it retains one seed representative and tries shifts using
another exact 2+2 join. Failure there concerns only the retained quadruple; other seeds with the
same PAF vector can have different cross-correlations.

### 3. Exact checks [M]

`logs/openai_4p3/small.json`, generated by `theory/openai_4p3_search.py` (5.4 s):

| p | distinct ordinary / switched seed vectors | complete seed gate | phase result |
|---:|---:|---|---|
| 3 | 2 / 3 | holds for explicit parameters | cores on both Q |
| 5 | 4 / 10 | holds for explicit parameters | cores on both Q |
| 7 | 9 / 35 | holds for explicit parameters | core on Z_2^2; chosen seeds fail on Z_4 |
| 11 | 45 / 175 | **fails, exhaustive** | impossible in this ansatz on either Q |
| 13 | 87 / 264 | holds for explicit parameters | chosen seeds fail on both Q, follow-up below |
| 17 | 156 / 496 | **fails, exhaustive** | impossible in this ansatz on either Q |
| 19 | 204 / 612 | **fails, exhaustive** | impossible in this ansatz on either Q |

At p=13 the first run hit its 2,000,000-pair cap. A separate bounded retry with a 5,000,000-pair cap
completed in 43.5 s (`logs/openai_4p3/phase13.json`). It rules out shifts for the selected seed quadruple,
not the whole p=13 ansatz. The phase lists were [343,2197,2197,2197] for Z_2^2 and
[322,2197,2197,2197] for Z_4. No search was left running.

All five positive cores were checked at every group lag by `theory/openai_4p3_verify.py`, which imports
no construction code. `theory/openai_4p3_check.py` then assembled the five matrices using the existing
border tooling and passed each to the independent verifier's Python-integer row-dot-product check of
**all N^2 entries of H H^T=N I**. Orders were 60 (both Q), 92 (both Q), and 124 (Z_2^2).
One-sign corruptions of each core and matrix were rejected. Logs: `checks.json` and
`matrix_p5_v4_check.json`; a complete small matrix is `matrix_p5_v4.json`.
These are known-order controls, not new existence results for those orders.

### 4. Coverage of all 16 target orders [M]

`theory/openai_4p3_search.py --coverage` computes all cubic traces and all attainable root-completion
sums, and solves (C) exactly. It supplies parameters for each witness in `logs/openai_4p3/coverage.json`.
The following tuple is (c_0,c_1,c_2,c_3); `openai_4p3_check.py` builds all 32 quotient choices and verifies
their quotient correlations. These are not full-core checks.

| 4n | p | compression witness | (C) | full (K),(S) |
|---:|---:|---|---|---|
| 2396 | 149 | (-11,-5,-3,1) | holds | unknown |
| 2876 | 179 | (-13,-3,1,1) | holds | unknown |
| 4124 | 257 | (-15,-5,-3,-3) | holds | unknown |
| 4348 | 271 | (-15,-7,1,1) | holds | unknown |
| 5308 | 331 | (-17,-7,-1,-1) | holds | unknown |
| 5756 | 359 | (-17,-5,-5,5) | holds | unknown |
| 5884 | 367 | (-19,-3,-1,-1) | holds | unknown |
| 6236 | 389 | (-19,-5,-3,1) | holds | unknown |
| 7324 | 457 | (-21,1,1,3) | holds | unknown |
| 7388 | 461 | (-21,-5,-1,-1) | holds | unknown |
| 7484 | 467 | (-21,-5,-3,-1) | holds | unknown |
| 8156 | 509 | (-21,-7,-5,-3) | holds | unknown |
| 8348 | 521 | (-21,-9,-3,1) | holds | unknown |
| 8380 | 523 | (-21,-7,-5,-5) | holds | unknown |
| 9148 | 571 | (-23,-7,-1,-1) | holds | unknown |
| 9596 | 599 | (-23,-5,-5,5) | holds | unknown |

### 5. Test for Claude and fixed pass/kill rules

First reproduce the small seed census and all five matrices with the commands below.
Expected exact zeros: p=11,17,19. Any positive there refutes this implementation's exclusion and requires
an independently checked counterexample. Any reconstructed example failing integer PAF or H H^T kills
the corresponding claimed control.

Next useful test: complete the seed gate at p=23,29,31, with phase enumeration disabled and a 540-second
limit. A complete empty gate kills this cubic-root-switch ansatz at that p, for both Q. A returned
quadruple passes only (K), and should motivate symbolic analysis of the displayed degree-six sums.
`capped` is **unknown**, never a zero. Failure of one selected phase join never kills the full ansatz.
The zeros at 11,17,19 already kill the claim of an all-prime family. A useful positive claim at a target
must still pass the separate integer core verifier, followed by independent bordered matrix verification.

## Supporting theorem: a border for every valid odd compression [P, M]

By the normal form, permute the four core sequences, negate whole sequences if needed, and translate
each on Q so that the three exceptional compression entries sit at q=0. These operations preserve
the sum of PAFs. The compressions become

    z_0=c_0 1,  z_j=c_j 1+2 delta_0 (j=1,2,3).

Use the existing Kumar/Claude Code border reduction with z*=(0,2 delta_0,2 delta_0,2 delta_0).
`theory/openai_4p3_border.py --certificates` supplies explicit X (12 by 12), Rh (12 by 16), Ch (16 by 12)
for both Q in `logs/openai_4p3/border_certificates.json`. The existing deterministic frame search finds
each at frame index 75 (the 76th frame). It checks, by integer arithmetic,

    X X^T=12 I,   Rh Rh^T=16 I,
    Ch Ch^T=I_4 tensor Gamma_Q,
    X Ch^T + Rh Mc(z*)^T=0,

with all entries +/-1 and each four-entry block of every Rh row balanced.
Here Gamma_Q has diagonal 12 and off-diagonal -4; Mc is the existing compressed GS array.

**Proof of transfer [P].** Mc(z)-Mc(z*) has constant blocks. It therefore annihilates Rh^T, so the linking
identity and Ch are unchanged by any c_i. Lift Rh and Ch constantly along K. The top row Gram is
12 I+16p I=(16p+12)I. The lower Gram is (16p+12)I by the core identity plus the displayed Ch identity;
the cross Gram vanishes by the linking identity. Thus the assembled binary matrix is Hadamard whenever
the input is a valid core. All possible odd Q-compressions are covered by the proved normal form.
This is a finite-certificate proof of border availability, not a proof of core existence.

## Ideas considered and rejected

- Pure translated or signed zero-completed quadratic characters in every coset: compression-energy obstruction above [P].
- Independently shifted low-index cyclotomic cells with only zero to three repairs each: the exact common-index gate excludes every applicable tested index at all 16 targets [M].
- A fixed-index, fixed-repair family for arbitrarily large primes: the explicit inequalities in Proposal 1 rule it out [P].
- The cubic-root-switch rule for every prime: the complete seed gate at 11,17,19 rules it out [M]; only conditional subfamilies remain.
- Treating a selected-seed phase failure as nonexistence: rejected because equal PAF signatures do not determine cross-correlations [P, follows from (S)].
- Direct group-homomorphic transport of the Z_522 core for H(2092) to the p=523 axis: every homomorphism Z_522 to Z_523 is trivial, since its image order divides both 522 and 523 [P]. This does not exclude a different product construction using that matrix.
- A new Galois-ring or conference-matrix product identity: none was obtained; no exclusion of those directions is claimed.

## What a Claude instance should run first

Python 3 and numpy suffice. Run from the repository root. All commands are deterministic; seed is null.
Each search invocation has a hard logical budget of at most 600 seconds and reports caps as unknown.

```sh
# 1. Arithmetic obstruction and target-only compression coverage.
python3 theory/openai_4p3_compression.py --out logs/openai_4p3/compression_python.json
python3 theory/openai_4p3_search.py --coverage --seconds 540 --out logs/openai_4p3/coverage.json

# 2. Full small seed census; first selected quadruple only for the phase step.
python3 theory/openai_4p3_search.py --seconds 540 --out logs/openai_4p3/small.json

# 3. Independent PAF and H H^T checks of every positive control, plus negative controls.
python3 theory/openai_4p3_check.py

# 4. Standalone construction and matrix verification, p=5, Q=Z_2^2.
python3 theory/openai_4p3_build.py --p 5 --q v4 --ansatz elliptic --params logs/openai_4p3/params_p5_v4.json --out logs/openai_4p3/core_p5_v4.json
python3 theory/openai_4p3_verify.py logs/openai_4p3/core_p5_v4.json
python3 theory/openai_4p3_border.py logs/openai_4p3/core_p5_v4.json --out logs/openai_4p3/matrix_p5_v4.json
python3 theory/openai_4p3_verify.py logs/openai_4p3/matrix_p5_v4.json --out logs/openai_4p3/matrix_p5_v4_check.json

# 5. Universal finite border certificates, already checked in this session.
python3 theory/openai_4p3_border.py --certificates --out logs/openai_4p3/border_certificates.json

# Optional reproduction of the completed p=13 selected-seed rejection (about 44 s here).
python3 theory/openai_4p3_search.py --primes 13 --max-pairs 5000000 --phase-max-p 13 --seconds 540 --out logs/openai_4p3/phase13.json

# Next unrun experiment: seed gate only; no large phase search.
python3 theory/openai_4p3_search.py --primes 23,29,31 --phase-max-p 0 --seconds 540 --out logs/openai_4p3/next_seeds.json
```

Expected outcomes of steps 1–5: the two minimum columns agree with `compression_host.json`; all 16 cubic
compressions pass; the small census matches the table; all five matrix controls pass with orders
60,60,92,92,124; both certificates have frame index 75. Any mismatch must be resolved before using
these claims. The next unrun experiment has the pass/kill criteria in Proposal 2; no outcome is asserted.

## Unfinished items

- No unbounded positive algebraic family, target-order core, or target-order Hadamard matrix was found.
- The cubic ansatz's degree-six character-sum identities have not been solved in closed form.
- Alternative seed quadruples and alternative PAF-equivalent representatives at p=7 (Z_4) and p=13 remain untested.
- Full seed and phase gates at the 16 target primes remain unrun; their compression witnesses are not candidates certified at the full group level.
- No obstruction to unrestricted cores, or to elliptic/quartic variants outside the two precisely defined classes, is claimed.
- The p=23,29,31 follow-up above is deliberately left to the testing instance; no job remains running.
