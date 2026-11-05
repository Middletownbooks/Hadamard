# Beating the deficit floor: block-exchanging symmetry, with exact gates

Date: 2026-11-05. Branch: `codex/order-3452-compression`.

**Status:** two completed research proposals: one gated constructive ansatz and one proved algorithm-class barrier.
No core or Hadamard matrix is claimed. No long search was run.

**Credit.** The bordered construction is Kumar's form of Theorem 7 of Egan–Kotsireas–London–Ó Catháin [EKLO].
The 863 setup and first certificate are by OpenAI Codex. The 523 solver and 863 lift analysis are by Claude Code.
This report uses the latter's `h863_analog.Host`, `hD_pair`, and `gauss_bits_mean` diagnostics.

**Tags:** [P] proof supplied; [M] finite exact computation, script identified; [E] estimate with assumptions.
Sampled Gaussian ranks and covariances are diagnostics, not exact ranks or lattices of the whole feasible set.
All costs below count pair visits or equivalent signature probes; conversion to time needs new benchmarking.

## Bottom line

1. **Best constructive proposal: exchange one pair of blocks under the S14 involution, instead of fixing all four.**
   Use `(A,B,C,tau C)` with only A and B individually tau-invariant. It is a different subset of S7 from the
   empty S14 subset. It passes an exact simultaneous Q-compression/lag-zero gate at 863: **40 surviving Q tuples**
   [M, `theory/openai_deficit_symmetry.py`]. No full H-level solution has yet been obtained.
2. **Conditional 863 exponent: 45.1–45.6 [E] across two synthetic seeds, not a demonstrated algorithm at that cost.** On one synthetic
   lag-zero-compatible y, the diagnostic gives D_def=23.71, final S-merge bits F=49.16, and D-generator entropies
   (8.72,8.93,17.18). If real stage-1 solutions retain these numbers, if one parent bucket supplies two sides
   of 2^29 usable distinct S entries, and if their generation is sufficiently cheap, then
   `2^23.71 * (2^(49.16+1-29) + 2^17.65 + 2^17.18)` has exponent **45.08**.
   Parent construction, signature costs and memory overhead are **additional**, with a remaining factor of
   about 1.9 before 2^46. Passing the gates below is mandatory. Confidence in feasibility is low; confidence in
   the algebraic definition and finite screens is high.
3. **At 519 this exact ansatz is empty already at the simultaneous Q/lag-zero gate [M]. Its true-target cost
   exponent is infinity.** A synthetic mean estimate of deficit 10.06 is therefore misleading. Use 519 for
   planted correctness tests, not a claim of a known core in this slice. Known existence of H(2076) says nothing
   about this restricted host. This is a useful early rejection, not evidence for a positive cost claim.
4. **A second result rules out a particular shortcut [P, proposal 2]:** no more than eight independent
   rectangular switches with a common nonzero mixed Fourier coefficient can leave the full correlation system
   affine in their switch bits. This is not a lower bound for all S/D splits or all nonlinear algorithms.

## Proposal 1 — block-exchanging equivariance, rather than four invariant blocks

### 1. Definition and exact algebra

[P] Let `tau(q1,q2,u,t)=(q2,q1,-u,t)`, and let `tau C(g)=C(tau g)`. Restrict the four sequences to

```
x0=A, x1=B, x2=C, x3=tau C;   tau A=A, tau B=B.
```

All three independent sequences remain S7-invariant. No condition `tau C=C` is imposed. Thus this is not S14;
S14's emptiness does not imply emptiness here. Conversely neither subset contains the other.

[P] Write R_x for PAF_x. A change of summation variable gives `R_(tau C)(g)=R_C(tau g)`. Hence the exact target is

```
R_A(g)+R_B(g)+R_C(g)+R_C(tau g)=F(g).
```

The left side is tau-invariant and so is F. Each omitted lag is therefore a copy of a retained equation. There
are no cross-block products: the third independent contribution is `R_C+tau R_C`, not `2 R_C` at a general lag.
Never substitute the latter. This identity is the reason the ansatz removes equations legitimately.

[P] The z=(1,1) translation commutes with tau. Therefore the usual sum/difference split is valid inside this
ansatz. For each independent generator X, set `S_X(g)=X(g)+X(g+z)` and `D_X(g)=X(g)-X(g+z)` on a transversal.
The fourth block's S and D are determined by the third. At the S-level use the exact coarse lag map
`R(g)+R(g+z)`. For the parent level sum this further over u in Z_m. The implementation constructs both maps
explicitly. The D-check is a **three-list** equality, contributions A, B, and C+tau C. Enumerate A+B versus C;
do not enumerate a fictitious independent fourth block.

[P] At the y-level, A and B have 11 independent H-cell values apiece; C has 20 and determines the fourth block.
There are 42 independent cells, with energy weights four 1s and thirty-eight 2s, rather than S14's 44 cells.
The class alphabet is `y=e+7s`, e=±1, s in {-6,-4,-2,0,2,4,6}. Signs within a fixed y-cell have exactly
`binomial(6,(s+6)/2)` possibilities. Fixing the original Codex y is **not** part of this proposal: it re-splits
at a newly solved y-level. The original compressed certificate and all its files are preserved.

### 2. Measured gates, estimates and total cost

[M] `theory/openai_deficit_symmetry.py` performs these finite, deterministic tests (no seed):

* Q character coefficients have squared sums `(4v,12,12,12)`. Every nontrivial coefficient is even and hence
  belongs to {-2,0,2}. Enumerate every allowed nonnegative trivial coefficient up to floor(sqrt(4v)), invert the
  4x4 Walsh transform, and require odd integral Q-cell sums. Independent block sign reversals justify the
  nonnegative trivial coefficients. The two fixed blocks are interchangeable. A,B require equality of Walsh
  coefficients 01 and 10; C's partner has these two coefficients exchanged.
* Enumerate all such triples A,B,C, by an exact join of their four squared-coefficient signatures. This yields
  48 tuples at 863 and 24 at 519, in this normalization (not a count modulo every symmetry).
* Check the complete alphabet against the weighted H-level lag-zero energy by finite subset-sum. Both orders
  pass this weaker test, with witnesses in the logs.
* For each Q tuple, impose its row sums **and** its energy. A fixed block's q=00,11 rows have weights
  `(1,2,...,2)`; its q=01 row has m free entries, repeated at q=10 in reversed order. A free C has four ordinary
  m-entry rows and contributes twice its energy. Dynamic programming records every reachable (sum,energy),
  then joins the three independent blocks' energy sets. Excess energy is 2856 at 863 and 1512 at 519.
  Result: **40/48 survive at 863, 0/24 at 519**. A second implementation in `openai_deficit_controls.py` exhaustively enumerates alphabet row tuples instead of using the reachable-state DP and agrees on every surviving tuple. These are necessary conditions, not full stage-1 decisions.

[E] `theory/openai_deficit_estimate.py` samples uniformly in a synthetic fixed-y fibre satisfying the ansatz
and lag-zero energy, with all nonsmall letters of absolute value 13. It does **not** satisfy the other H-level
lags or Q conditions. Read its Gaussian at the sample mean, not at the true target. Exact row enumeration gives
conditional Shannon entropy H(D|S), not the log-size of every S-solution's D-fibre.

| diagnostic, seed 863, N=4000 | value | status |
|---|---:|---|
| log2 fibre given synthetic y | 177.786 | [M] exact product of row-option counts, log rounded |
| C(all), sampled rank | 144.93, 32 | [E] lattice-Gaussian at mean |
| C(S), sampled rank | 86.39, 18 | [E] |
| C(parent), sampled rank | 37.23, 6 | [E] |
| F=C(S)-C(parent) | 49.16 | [E] chain-rule diagnostic |
| hD(A), hD(B), hD(C) | 8.719, 8.933, 17.175 | [M] finite row-pair entropy sums, rounded |
| D_def | 23.713 | [E] 144.93-86.39-34.827 |
| expected cores per synthetic y | 2^32.856 | [E] inapplicable to existence until actual y is found |

[E] A second synthetic placement (seed 864, N=4000) gives fibre exponent 177.371, C_all=145.04, C_S=86.24, C_parent=37.26, h_D=34.373, D_def=24.427, F=48.98, and log2 expected cores=32.331. The same n=29 idealized substitution costs **2^45.60**, leaving only a factor 1.32 before the budget. It fails the desired one-bit uncertainty margin even before parent overhead. These two seeds illustrate sensitivity; they do not bound variation over real y.

[E] For **real** S-solutions measure their list sizes L_A,L_B,L_C, their D-completion probability p_D, the parent
cost P per S-solution, the effective residual match bits F, and usable side exponents n_L,n_R. With balanced
n_L=n_R=n, the proposed cost formula is

```
T_core = (1/p_D) * [ P + 2^(F+1-n) + L_A*L_B + L_C ].
```

The two D lists use O(L_A L_B) memory and exact integer signature comparison; hash collisions require full
rechecking. The final S-join formula counts both sides, assumes uniform effective signatures in the chosen
parent bucket, and requires generation of enough distinct candidates. Repeated representations count only
once. For unbalanced lists use `(2^n_L+2^n_R)/2^(n_L+n_R-F)`, not the balanced formula.

[E] Substitution of the synthetic diagnostics and **assumed n=29, P=0** gives exponent **45.08** at 863. This
is a conditional feasibility target, not an achieved upper bound. With n=28 it is about **46.00**, leaving
essentially no overhead. At n=30 it is about **44.25** but memory becomes harder. An 8-byte hash plus 8-byte
reconstructible candidate identifier at n=29 uses 8 GiB for the stored side, plus sorting/workspace; stream
its counterpart. If identifiers cannot reconstruct all choices, charge their actual size. Partitioning or
regeneration is not free: include every repeated pair visit. A few nanoseconds per lookup is not assumed.

[E] The proposed gain spends roughly 161.3 bits of the original 339.1-bit y-fibre while reducing the conditional
constraint diagnostic by roughly 130.7 bits. It still leaves roughly 32.9 expected cores rather than 63.5.
These compare different synthetic/actual fibres, so only the method, not their difference to a decimal, is
meaningful. No production cost follows from this comparison alone.

### 3. Which floor assumption changes

[P] The deficit estimates in the existing report concern four independently chosen S7 blocks, or four
individually fixed S14 blocks. Here an automorphism exchanges blocks and preserves their **sum of PAFs**.
Individual-block invariance is not needed for aggregate invariance. This changes both the equations and the
independent lists. It does not claim an entropy discount while retaining the original 2^174 fibre.

[E] The potential gain is a D-deficit near 24, a roughly 49-bit final S match, and a D check with two sides near
2^18, instead of the original 210-bit mixed match. Actual S-solution conditioning could destroy this gain:
the 523 energy-pinning measurements are precisely why the synthetic estimate is insufficient.

### 4. Cheap tests, pass and kill criteria

**First gate, 863, at most two CPU-hours; individual solver calls at most 300 seconds.**

[P, specification] Iterate `joint_gate(5,q_gate(5))` in its returned order. For each tuple build an integer
stage-1 model with 42 variables from the exact alphabet. Expand into four length-20 arrays by the ansatz,
fix the four Q row sums of each block, and impose on **all 20 H lags**

```
sum_i sum_(q,u) y_i(q,u)*y_i(q xor a,u+b) =
    2936 if (a,b)=(0,0); -516 if a=0,b!=0; 172 otherwise.
```

All products are exact integer multiplication constraints (CP-SAT AddMultiplicationEquality, or exhaustive
lists as in `h863_s14_stage1_exact.py`). Use min-energy bounds from the finite DP. Allow one worker, seed
863+i for tuple i; log status, elapsed time, tuple, and every certificate. Independently verify an emitted y
using `h863_lib.check_y`, its alphabet and the block relation. This completely specifies the first experiment;
no implementation choice should relax the constraints. A 300-second UNKNOWN is not INFEASIBLE.

* **Pass to costing:** at least one independently verified real y. Prefer 4 distinct y if inexpensive.
* **Kill the ansatz rigorously:** all 40 tuples decided INFEASIBLE, with exhaustive or solver certificates/logs.
* **Operational stop:** no y within two CPU-hours; record UNKNOWN and stop solver construction for this ansatz.
  Do not report an emptiness theorem from timeouts.

**Second gate, real-y costing, at most two CPU-hours.**

[E, specification] Replace synthetic E,S in `openai_deficit_estimate.py` by E=y mod 14 (representatives ±1),
S=(y-E)/7. Retain its generator/orbit sampler. Use seeds 863,864,865 and 8000 samples each; read all levels at
the true targets with `q523_gauss.gauss_bits`. Check small-dimensional tails against exact hit counts. Use the
z=(1,1), then Z5 parent maps already implemented. Enumerate or sample the parent triples by exact signature
joins: pair A-parent+B-parent against the contribution C-parent+tau C-parent. For each fixed matching middle
signature, collect descendant S candidates on both sides. Keep the origin identifiers, deduplicate, and
measure effective side sizes and generation visits. This is a three-generator join, not the existing four-list
kernel with C inserted twice.

* Require real-y estimated log2 expected cores >=10, D_def <=25, and stable tail fits (within 1 bit where
  exact hits are numerous). These are screening thresholds, not sufficient conditions.
* Ultimately require `log2(T_core)<=46` **including** measured parent/list construction, D checking and memory
  overhead, with at least a one-bit uncertainty allowance for the extrapolated probabilities. Thus the
  45.08 synthetic illustration barely has the desired margin and is not itself a pass.
* Measure available descendants before assuming n=29. If any necessary list capacity or cost bound fails,
  kill this solver design even if the ansatz remains mathematically viable.
* When sampling S-solutions, use their actual D-list sizes. Do not substitute unconditional H(D|S) for them.

**Planted controls (run here).** [M] `openai_deficit_controls.py` constructs random ansatz sequences with nine
independent rectangular switch bits, sets targets from a uniformly selected point of the 512-point face,
and enumerates the entire face without injecting a privileged candidate. Both y and w stay fixed. Integer
PAFs at every S7 lag class agree with FFT-rounded values for all 512 points. The plant is recovered uniquely
at m=5 (mask 477, seed 863) and m=3 (mask 115, seed 519). This checks exact algebra and recovery, **not**
full-fibre complexity or a true-target core. A future production join must repeat planted tests with plants
not forced into privileged middle buckets; report inclusion probability separately from recovery conditional
on inclusion.

**519 true-target test (run here).** [M] The complete Q/energy gate is empty. Kill this *specific* ansatz for
519 immediately. Do not spend hours on its synthetic promise. An unrestricted 519 S7 solver is still a valid
control for general machinery, but cannot validate this restricted ansatz's true-target rate.

### 5. Transfer

[P] The algebra applies to any target-preserving involution tau normalizing the multiplier group and commuting
with the chosen order-two translation z. Start with the 863-shaped groups in `h863_lift_answer.md` §4 using
q-swap and u inversion. Recompute the exact alphabet, quotient PSDs and joint energy gate for each p,f,m.
No claim that those gates pass is made. The existing 38–90-bit deficits cannot simply be halved numerically.

[E] Priority transfer hosts: 1583 and 2063 (also six classes), then 583 as a second small analogue. For the rich
hosts in `b3_push_answer.md` §8.1 (643,739,751,839,955,971), first enumerate target-preserving normalizer
involutions. If there is no suitable commuting z, this particular staged solver does not transfer unchanged.
This finite screen is preferable to repeating an unmodified annealer on those hosts.

## Proposal 2 — a barrier to independent-switch linearisation

### 1. Definition: the algorithm class being ruled out

[P] Consider a Boolean face of the fixed-(y,w) table fibre parameterized by independent signs d_j in {±1}:

```
x(d) = x_bar + sum_(j=1..k) d_j v_j.
```

Each direction is a nonzero rectangular table switch: in one block and one q-cell it has opposite signs on
rows u!=v and columns a!=b of the m x 6 table. Class entries expand to the seven corresponding t-coordinates;
singletons are unchanged. Disjoint supports and alternating entries make all 2^k orientations valid signs
with the same row and column sums. The proposed shortcut would fix a background and solve the complete
correlation constraints by linear algebra in all k switch signs, with no quadratic products or extra
quadratic compatibility check. We rule out k>8. More generally the proof applies to any affine cube whose
directions all have nonzero Fourier image at one common character.

### 2. Proof of the eight-switch barrier

[P] For any one mixed character chi, collect the four block transforms in C^4, viewed as R^8. Put
`a=hat(x_bar)(chi)` and `b_j=hat(v_j)(chi)`. The target includes

```
||a + sum_j d_j b_j||^2 = 3452.
```

After using d_j^2=1, its coefficient of d_j d_l for j!=l is `2 Re <b_j,b_l>`. The Boolean monomials
`1,d_j,d_j*d_l` are linearly independent as functions on the entire cube (multiply by a monomial and average
uniformly to recover its coefficient). Thus the norm is affine on the cube only if the nonzero real vectors
b_j are pairwise orthogonal. There are at most eight such vectors in R^8. If all moves lie in a single block,
there are at most two in that block's R^2 plane. Consequently a wholly affine four-block rectangle scheme has
at most two switches per block, at most eight total. This is a necessary condition, not a sufficiency theorem.

[P] Every rectangle in this problem has a nonzero transform at **every** mixed character. Up to a nonzero
constant its transform is

```
psi(q) * (zeta_m^(r*u)-zeta_m^(r*v)) * (eta_a(s)-eta_b(s)).
```

For m=5 or 3 and r!=0, the first difference is nonzero. For s!=0 mod43, if the period difference were zero,
the nonzero integer polynomial
`P(X)=sum_(t in C_a) X^t - sum_(t in C_b) X^t`
would vanish at the primitive root zeta_43^s. It has degree <=42 and constant term 0. The minimal polynomial
is `Phi_43=1+X+...+X^42`; hence P would be a rational constant multiple of Phi_43. Its zero constant term
would force P=0, contradicting disjoint nonempty classes. So all rectangles share a nonzero mixed character,
and the dimension argument applies irrespective of which rows, columns, blocks or q-cells they use.

[P] This remains a barrier if equations are replaced by invertible linear combinations: if every new
constraint is affine, every old one is affine. It does **not** rule out eliminating some variables using
nonlinear relations, changing variables nonlinearly, branching on quadratic products, using several
backgrounds, exploiting Fourier-disjoint non-rectangular directions, or an algorithm that handles a
non-affine residual. The lag-zero norm alone would give only disjoint-support orthogonality; the mixed
character norm is essential to this stronger bound.

### 3. Cost and the assumption it exposes

[P] A valid k<=8 affine switch face contains at most 256 configurations. Exhaustively covering an arbitrary
fixed fibre of size 2^h by such faces requires at least `ceil(2^h/256)` faces, even with overlap allowed. This
is a counting bound on this covering strategy, not a general computation lower bound. At the supplied
863 median h=173.9 it is at least 2^165.9 faces per (y,w); at the 519 diagnostic h=70.7, at least 2^62.7.
The measured 863 fibres vary, so use their actual h in an implementation.

[E] Under the repository's homogeneous random-target approximation, a face with at most 2^8 points has
expected hits at most `2^(8-C)` if its placement has no extra target bias. The most optimistic face-decision
cost is then `T_core ~ 2^(C-8) * t_face`: about **2^202 t_face at 863**, using C=210, and **2^84.4 t_face at
519**, using C=92.4. These are estimates for untargeted face sampling, not a proved lower bound on an
adaptive algorithm's success. Reusing backgrounds or selecting target-correlated faces requires a new
probability/cost analysis. Merely saying "the D problem is linear" does not supply one.

[P] The attractive but invalid assumption is that many independent margin-preserving switches can be moved
into D while retaining an affine D system. Quadratic cross terms prevent that. The result is useful because
it rules out this entire family of unrestricted Boolean-cube linearizations before solver engineering.
It does not prove that no split of S7 has deficit below 30.

### 4. Cheap tests and pass/kill criteria

[M] The two planted controls already supplied each enumerate a nine-switch face using exact integer PAFs.
Each finds nine nonzero mixed second-difference vectors among 36 switch pairs. Their values are computed by
`V(0)-V(e_j)-V(e_l)+V(e_j+e_l)`; nonzero is an exact witness that the corresponding Boolean product cannot
be omitted. The tests recover the planted target, so they do not confuse failure of linearization with
absence of solutions. Their full PAF vectors also independently check the FFT sampling convention.

[P, specification] To test any proposed implementation of this shortcut, choose 20 seeds, generate independent
rectangles in the actual 863 or 519 fibre using `h863_lib.TableSampler` (or the analogous 3x6 sampler), retain
only disjoint supports, and evaluate these integer second differences at all lag classes. For a proposed
k>=9 wholly affine cube, **kill the proposal on the first nonzero second difference**. A zero result for all
pairs at k>=9 means either a direction is zero, the directions are not independent, the test omitted lags,
or the implementation is wrong, by the proof above. For k<=8, zero second differences pass only the algebra
check; require the cost inequality above or a measured target-bias improvement before proceeding.

[M] At the real 519 target the proposal-1 ansatz is empty, but the affine-switch theorem also applies to the
unrestricted 519 S7 tables. Its planted controls remain valid. No benchmark is allowed to infer a true-target
hit rate from forced inclusion of a planted candidate.

### 5. Transfer

[P] For four real sign sequences the mixed-character space is always R^8. The same result applies to any
prime m and prime p with distinct nonempty multiplier classes, including 1583 (p=79), 2063 (p=103), and the
other Z2^2 x Z_m x Z_p hosts of the transfer table. The polynomial proof uses p prime, not p=43 specifically.
For composite m, restrict to rectangles whose row difference is detected by a common chosen character.
For rich hosts in `b3_push_answer.md` §8.1 with different group shapes, the common-frequency nonvanishing
hypothesis must be checked; it must not be assumed. For b blocks the dimension bound is 2b at a complex
character, or b at a real character with nonzero directions.

### Related representation accounting (limited negative result)

[P] Writing x=a+b gives the exact identity
`PAF_x=PAF_a+PAF_b+CC(a,b)+CC(b,a)`. Thus an additive join on the two PAF signatures alone discards a required
term. At a character that term is `2 Re(hat(a)*conj(hat(b)))`. Rectangular overlaps do not make it vanish;
the theorem explains when a whole family could make such terms vanish simultaneously.

[P] If a physical solution has R representations and a filtering rule retains each representation with
probability q, the probability of retaining the solution is at most min(1,R*q), by the union bound, without
assuming independence. Counting R as a free R-fold gain after choosing filters of probability 1/R is wrong.
A legitimate representation algorithm must state its filters, cross-term signatures, representation
multiplicities, and actual list-generation work, and count each reconstructed physical core once.

This does **not** rule out BJMM/May–Ozerov-style improvements with a correctly accounted representation
scheme. No such complete low-cost scheme is established here; it is omitted from the ranked proposals
rather than assigned an unsupported favorable exponent.


## Ideas considered and rejected

* Four individually tau-invariant blocks: the repository's S14 stage-1 emptiness already decides them.
* The original y->w lift with a faster generic inner join: does not remove its 210-bit residual match.
* Counting x3=tau C as an independent block in a Wagner tree: invalid independence and duplicate completions.
* Calling the 519 synthetic estimate evidence of success: the exact Q/energy gate disproves that inference.
* Two block pairs exchanged by the same involution: at its fixed nontrivial Q character the target would be 2a^2+2b^2=12 with a,b even, impossible modulo 8 [P].
* A 3-cycle of blocks induced by an order-three automorphism at 863: the trivial-character equation gives 3a^2+b^2=863 after dividing row sums by two, impossible modulo 3 [P].
* A 4-cycle of blocks: equal row-sum magnitudes force their square to equal v, impossible for v=863 [P].
* Lag-zero pass alone: necessary but weak; even our 519 ansatz passes it and fails the joint gate.

## What a Claude max-effort instance should run first

Run the exact 863 stage-1 model of proposal 1, section 4. Do not implement a large lift solver first.
The expected useful outcome is either a real y certificate or a complete emptiness decision. UNKNOWN is a
bounded unsuccessful experiment. The exact 40-case specification above is the handoff.

Reproduce the already-run cheap screens from the repository root:

```bash
python3 theory/openai_deficit_symmetry.py --m 5
python3 theory/openai_deficit_symmetry.py --m 3
python3 theory/openai_deficit_estimate.py --m 5 --N 4000 --seed 863
python3 theory/openai_deficit_estimate.py --m 5 --N 4000 --seed 864
python3 theory/openai_deficit_estimate.py --m 3 --N 4000 --seed 519
python3 theory/openai_deficit_controls.py --m 5 --seed 863
python3 theory/openai_deficit_controls.py --m 3 --seed 519
```

Logs are under `logs/openai_deficit/`. Commands, seeds and output paths are indexed in `logs/openai_deficit/manifest.json`.

## Unfinished items

* No full real-y gate for the new 863 ansatz has been run; this is the deliberately short-session handoff.
* No actual S-parent capacity, conditioned D-completion probability, or complete solver timing has been measured.
* No other transfer host has been screened with this ansatz. No Hadamard matrix was found or asserted.
