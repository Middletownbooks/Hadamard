# Algebraic proposals for the k4 bordered family at n = 4p + 3

Status: initial checkpoint, 2026-11-05. Work in progress; no new core or Hadamard matrix is claimed.

Tags: [P] proof supplied; [M] exact computation with the script and scope named; [E] estimate with method named. Untested proposals are explicitly conditional.

Credit: the bordered construction is Theorem 7 of Egan–Kotsireas–London–Ó Catháin (EKLO), in Kumar's form (https://github.com/vatskishan03/hadamard). The census, pre-existing obstructions, calibration and border tooling are by Claude Code in this repository.

Setup: GitHub's comparison found codex/order-3452-compression at 5b26dae60414599ccccf4c692220f616fc2aae24 to be an ancestor of claude/focused-planck-j1645f at a87392b4d91fcaafe0bc61d85435b4ff0a06c878 (nine commits ahead). The target branch was fast-forwarded to that reference commit with an expected-head check. The shell execution environment was unavailable, so this used the GitHub API, not a local git checkout. No PR is opened and no other branch is written.

## Bottom line

[P] Proposal 1 rigorously excludes small repairs to balanced character seeds and gives an exact compression gate for independently translated cyclotomic cells. [M] At the 16 targets the balanced version needs at least 29–61 total sign changes, or a common per-cell cap of at least 4–7. Full-core existence remains open. A second, conditional elliptic-character proposal is in progress.

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
The planned builder will implement this rule literally.

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
d(z)=(|z|-1)/2 sign changes to reach odd sum z. This number is attainable at the sum level.
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

[M] The embedded HOST_JS in theory/openai_4p3_compression.py was executed using V8, with integer arithmetic
(all values are far below 2^53), no random seed. logs/openai_4p3/compression_host.json stores the complete
minima and witnesses for all primes from 5 through 97 and the 16 targets. The pure-Python implementation is supplied
for independent reproduction; it has not yet run because this session has no shell runtime.

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

In progress. Small cases and a shift-independent seed gate are being checked. No open-order success is claimed.

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
