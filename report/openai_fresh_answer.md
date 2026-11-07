# A faster single-correlation oracle, and the limits of the resulting search

Final research report, 2026-11-07. Answers `report/prompt_openai_fresh_idea.md`.

**Result:** one algorithmic candidate is developed below, with a proof, executable
code and negative calibration evidence. The useful result is an exact
integer-field optimizer for one lag, reducing its cost from quadratic to
near-linear in the sequence length. The associated whole-design search has **not**
earned a production run at an open order. No new open-order Hadamard matrix,
Paley-type identity, or unconditional sub-2^46 solution-cost claim is made.
There are not three credible proposals to rank; padding the list would conceal
this outcome. Rank 1 below is a research candidate, not an endorsed solver.

[P] means a proof given here or explicitly attributed. [M] means an exact
computation, with its script and log identified; wall-clock measurements are
observations, not complexity proofs. [E] marks an estimate, reconstruction or
untested hypothesis. Float-driven searches are heuristic even when their outputs
receive exact verification.

The designated branch was fetched and the Claude branch merged cleanly before
the first checkpoint. Initial and experimental checkpoints were published on
`codex/order-3452-compression`; no other branch or pull request is a deliverable.
The source audit is against Claude commit
`a202d5184c5fe6bec0f9f671c2630a56a2247281` and the merged repository.

## 1. Rank 1: exact lag constraints through cyclic domain walls

### 1.1 What is searched

Let G be a finite abelian group of size m and let four functions
x_b:G -> {−1,+1} have prescribed total periodic autocorrelation T(d):

    R(d) = sum_b sum_g x_b(g) x_b(g+d) = T(d).

For each nonzero lag representative modulo the common multiplier group and
inversion, give every full group coordinate its own local binary copy. The
translation g -> g+d partitions G into cycles. Across all four blocks, let w_d
be the number of sign changes around these cycles, including each closing edge.

**[P] R(d)=4m−2w_d.** Thus the local lag constraint is exactly a prescribed
number w_d=(4m−T(d))/2 of domain walls. Each cycle has an even number of walls.
For an order-two translation both directed edges count; this convention is
essential and is used by the code.

The search separates two sets:

- A: every lag's local copy is binary and has its exact wall count. Separate
  binary cardinality factors impose specified block sums or compressed sums.
- B: all copies of a coordinate agree, coordinates in the same multiplier orbit
  agree, and optionally their full Fourier vectors satisfy the joint spectral
  norm constraints. The last option is the final tested formulation.

The independent local copies make projection onto A tractable. Concurrence in B
restores the shared signs. Their intersection gives the desired core. There is
no enumeration of compressed completions and no assumption that single-sign
changes decrease a global energy.

The supplied hosts are fully specified:

| H order | Group G, m | Common multiplier | Target T | Additional restriction |
|---|---|---|---|---|
| 4412 | Z_1103, m=1103 | unique subgroup of order 19 in F_1103^* | T(0)=4412, zero elsewhere | signed block sums (1,37,39,39) |
| 3452 | Z_2^2 x Z_5 x Z_43, m=860 | order-7 subgroup acting only on Z_43 | T(0)=3440; −12 on K\{0}; +4 off K, K={0}^2 x Z_5 x Z_43 | Z_43-compression y from `data/h3452_stage1.json` |

[M] Both orders are in `data/open_orders_lt10000.json`. The 4412 row profile is
the only absolute profile surviving the multiplier row-sum gate; see
`openai_fresh_audit.py`. This gate is necessary, not an existence certificate.
For 3452 the supplied constructor checks the exact compressed PAF identity
before using y. Its fibre has 560 orbit signs; its estimated richness is inherited
from `report/h863_lift_answer.md`, not newly proved here.

### 1.2 The exact optimizer and its proof

Consider any disjoint union of cycles with M vertices. For integer fields h_j,
write

    F(k) = max { sum_j h_j x_j : x_j in {−1,+1}, w(x)=2k }.

**Lemma [P]: F is discretely concave on its feasible integer interval.**
Take maximizing assignments x,y with wall counts 2(k−1),2(k+1). On each cycle,
the set where x differs from y is a union of consecutive runs, or the whole
cycle. Swapping the x and y values on one such run leaves their pointwise sum
unchanged. Only its two boundary edges change; the wall count of x changes by
−2,0 or +2 and that of y by the opposite amount. A whole-cycle swap changes
neither wall count. Since w(y)−w(x)=4, the sum of the run contributions is +4,
so some run contributes +2. Swapping it gives two assignments with wall count
2k and the same total linear objective. Hence

    2 F(k) >= F(k−1) + F(k+1).

The proof also works across multiple cycles, with no assumption on the fields.
All k from zero through the sum of floor(cycle_length/2) are feasible.

**Exact algorithm [P].** Put S=sum|h_j|. The slopes F(k)−F(k−1) are integers
between −2S and 2S. For an integer penalty lambda, solve

    max_x [ 2 sum_j h_j x_j − lambda w(x) ].

This is a two-state dynamic program on each cycle: fix its first sign, propagate
the best objective for each last sign, and close the final edge. Among objective
ties, compute both the minimum-wall and maximum-wall optimizers. Binary-search
lambda in [−2S−1,2S+1]. Concavity and integral slopes guarantee a penalty at
which the desired wall count lies between those two optimizers' counts.

If the counts straddle the desired count, swap +2-wall disagreement runs from
the smaller-count optimizer toward the larger. Each swap preserves the sum of
the two penalized objectives. Since both were global maximizers, neither new
assignment can have a smaller objective without the other exceeding the maximum;
both remain maximizers. There are enough +2 runs because the total count
difference is positive. Stop at the prescribed even count. This yields an
exact maximizer of the original constrained linear objective.

A cycle oracle takes O(M) time and O(M) memory. The complete integer algorithm
uses **O(M log(1+S)) time and O(M) memory**, versus O(Mw) time and traceback
storage for the counted dynamic program when w is of order M. This is a proved
reduction in the polynomial cost of a local projection, **not** a proved change
in the exponential cost of finding a complete Hadamard core.

`theory/openai_fresh_parametric.py` implements the integer oracle in embedded
C++17, using checked int64 arithmetic. It bounds the fields conservatively so
all penalty scores stay far from overflow and the unreachable-state sentinel.
`theory/openai_fresh_cycle.py` contains a separately implemented counted dynamic
program used as an independent optimization reference.

**Real-field qualification [P].** The reflection search supplies real fields v.
The fast implementation uses h=round(qv), with
q=2^20/max(1,max|v_j|/1024). The integer oracle is exact for h, not necessarily
for v. If each rounding error is at most 1/(2q), the linear-objective regret is
at most M/q, and the squared-distance regret at most 2M/q. FFTs and floating
arithmetic add their usual numerical errors. No claim of exact real Euclidean
projection is made for this implementation.

### 1.3 The outer iteration, completely specified

Initialize every replica with an independent standard normal value, using
NumPy `default_rng(seed)`. Apply

    a = P_A(z)
    b = P_B(2a−z)
    z = z + beta (b−a),             beta = 1.

P_A applies the cycle oracle to each lag, and assigns +1 to the largest k
coordinates in a cardinality factor requiring k positive signs.

For the final spectral P_B, average all occurrences of each orbit variable,
expand to four arrays on G, take their FFT, and at every character scale the
four-component complex vector to norm sqrt(FFT(T)). At a zero vector choose
the positive real first component. Fix the DC vector to the requested row sums
in the b=0 case. In the 3452 case fix the whole Z_43-trivial frequency slice to
FFT(y), with the FFT on Z_2^2 x Z_5. Invert, average within orbits to remove
roundoff, and copy back to all replicas.

[P] Before numerical error, this P_B is a nearest-point projection for the
replica metric: every full coordinate has the same occurrence count, checked by
the constructor; common multiplier invariance and conjugate symmetry survive
the radial Fourier projection. Prescribed DC/compression values already have
the required norms. The target squared Fourier norms at 3452 are 3452 except
at the three nonprincipal characters of G/K, where they are 12; the principal
one is 3452.

The earlier linear P_B is also retained: weighted concurrence followed by the
exact weighted affine row/compression projection. It did not solve the control
problem either. The weak endpoint-only relaxation in `openai_fresh_concur.py`
is preserved as an unsuccessful intermediate experiment, not a second proposal.

Every 10 iterations, round the orbit averages of a and inspect the PAF residual.
A zero float residual is merely a trigger for a direct integer PAF test. The
script records valid b=0 cores outside the prescribed signed row profile
separately; these are still valid GS cores and must not be discarded as failed
Hadamard constructions. Full matrix verification is a separate program.
A fixed point of the iteration lies in A intersect B; there is **no** convergence
or positive success-probability theorem.

### 1.4 What changes, and what does not

[E] The motivation is to satisfy a whole lag exactly, including its cycle
closure, before reconciling lags. This exposes a stronger local constraint than
independent edge products, while the parametric algorithm makes it affordable.
It does not pay the repository's exact staged lift deficit on every attempt and
can use a larger fixed fibre. These are reasons to test it, not guarantees that
its basins of attraction contain useful solutions.

This specialization and integration are new to the audited repository. General
parametric Ising methods, constrained segment optimization, divide-and-concur,
and reflection algorithms are established ideas; global priority is **not**
claimed. In particular the joint Fourier-sphere projection is already
Proposition 3.9 of Artacho, Campoy, Kotsireas and Tam, cited below. Merely replacing
`scripts/s4_rrr.py`'s moving per-block targets with that fixed sphere does not
meet the user's novelty requirement. The cycle oracle is the substantive
algorithmic ingredient retained after that literature check.

## 2. Exact controls, measurements, and unsuccessful search tests

[M] `logs/openai_fresh/parametric_controls.json` records 300 exhaustive
integer-field projections, 84 complete concavity tables, 60 comparisons against
the independent counted DP on larger cycles, and rejection of an unsafe integer
range. All passed. `cycle_controls.json` records 45 additional exhaustive
small comparisons for the floating-field counted DP, with a 1e−9 tolerance.
These validate implementation; the general theorem is the proof above.

[M] `logs/openai_fresh/audit.json` compares the two optimizers on identical
integer fields at **every lag factor** in each open-order benchmark. All optimum
values agree. Timings exclude compilation/model setup and are short samples on
this machine, not portable speed guarantees.

| Host | Orbit signs | Replicas | Lag factors | Counted-DP updates | Parametric updates | Counted / fast projection seconds | Complete fast iteration seconds |
|---|---:|---:|---:|---:|---:|---:|---:|
| 4412 | 236 | 132,360 | 29 | 1,128,500,896 | 20,876,288 | 1.706 / 0.203 | 0.239 |
| 3452 | 560 | 247,680 | 71 | 1,614,648,064 | 37,089,280 | 2.134 / 0.366 | 0.382 |

The update count drops about 54x and 44x; elapsed projection time improves about
8.4x and 5.8x. At tiny lengths the parametric method can instead be slower than
the quadratic DP, because its repeated oracle calls outweigh the asymptotic gain.

Selected search results follow. All runs in this session were under ten minutes;
the longest listed search took about 98 seconds. Seeds are 0 unless the log says
otherwise. Trials within a run use successive draws from one seeded generator.

| Method and log prefix | Control / iterations per trial | Result |
|---|---|---|
| Fixed joint sphere, A-first (`sphere_Afirst_43`) | p=43, d=3, rows (11,5,5,1); 12 x 20,000 | 8 cores |
| Fixed joint sphere, A-first (`sphere_Afirst_61_profile`) | p=61, d=3, rows (11,11,1,1); 12 x 100,000 | 1 core |
| Fixed joint sphere, A-first (`sphere_Afirst_131`) | p=131, d=5, rows (21,9,1,1); 12 x 20,000 | 0 cores; this specific profile was not established as a positive control |
| Fixed joint sphere, B-first (`sphere_k4_863`) | supplied 3452 fibre; 4 x 10,000 | 0 cores |
| Endpoint replicas (`concur_43`, `concur_61`) | p=43 fixed rows / p=61 free rows; 4 x 20,000 each | 0 hits in the requested constraints; earlier p=43 out-of-profile cores were not retained |
| Exact counted cycles, spectral concurrence (`cycle_spectral_43`, `cycle_spectral_61`) | same fixed p=43/p=61 rows; 4 x 10,000 each | 0 prescribed-row hits; p=43 attained PAF zero outside that profile |
| Parametric cycles, spectral concurrence (`parametric_43_recorded`) | p=43 fixed rows; 4 x 10,000 | 0 hits with the prescribed signed rows; **2 valid GS cores outside those rows** |
| Parametric cycles, spectral concurrence (`parametric_61`) | p=61 fixed rows; 4 x 10,000 | 0 cores |

The two retained p=43 row vectors are (11,1,7,1) and (11,5,5,−1). The second can
be normalized to the requested rows by a whole-block sign flip. The first is
also a perfectly valid order-172 GS core. The `solved` counters use the stricter
signed-row test; they must not be read as the number of all Hadamard cores found.
`parametric_43_recorded.json` repeats the same seeds as `parametric_43.json` only
to preserve these outputs; it is **not** four additional independent trials.

[M] `theory/openai_fresh_verify.py` imports none of the search modules. It
recomputes correlations directly, constructs the GS matrix independently, checks
**HH^T=nI in int64 arithmetic**, and confirms that flipping one matrix entry is
rejected. `independent_matrices.json` records eleven passing matrix checks at
known orders 172 and 244: eight sphere p=43 outputs, one sphere p=61 output,
and both retained parametric p=43 outputs. This is eleven verified outputs,
not a claim of eleven inequivalent designs.

**Decision:** the new local optimizer is validated; the full iteration is not
validated for scaling. Four failures at p=61 do not prove failure probability
one or nonexistence. They do give no basis for assuming an adequate success
rate at p=1103 or in the 3452 fibre. Stop production scaling at this checkpoint.

## 3. Cost per solution: conditional arithmetic, not a promised budget

Let I be a restart cap and p_I its probability of finding a valid core. With
independent restarts and per-iteration work c, the expected work is at most
c I / p_I, conditional on p_I>0. Early stopping can lower it. Neither p_I nor
its dependence on m is known. In particular a positive expected core count in
a fibre is not an estimate of p_I.

For L lag factors and fields bounded by B, the projection work is
O(L m log(mB)); the previous counted DP is O(L m^2). Storage for the replicas,
working vectors and parametric traceback is O(Lm). The measured hosts use only
a few megabytes per vector; [E] under 1 GB per worker is a generous allowance,
so four workers fit 15 GB without storing candidate pools.

[E] Charge 16 elementary operations per DP state update, plus Fourier,
concurrence, cardinality and bookkeeping overhead. Rounding upward gives a
working estimate c=2^29 at 4412 and c=2^30 at 3452 for the measured fields.
Use **2^30 and 2^31** respectively as a factor-two planning allowance.
These are explicit operation-accounting estimates, not measured CPU instruction
counts or rigorous worst-case constants. Integer penalty-search bounds also
remain logarithmic when the fields grow; quantization caps each field near 2^30.

| Assumptions using planning c | 4412 expected operations | 3452 expected operations |
|---|---:|---:|
| I=10,000, p_I=1/4 | about 2^45.29 | about 2^46.29 |
| I=100,000, p_I=1/10 | about 2^49.93 | about 2^50.93 |
| Minimum p_I to meet 2^46 with I=10,000 | 0.153 | 0.306 |

Thus even the optimistic 25%-within-10,000 assumption misses the planning budget
at 3452. At the measured 0.239/0.382 seconds per iteration it would imply about
2.7/4.2 CPU-hours per success before the factor-two allowance; four cores change
elapsed time, not the operation budget. These examples are sensitivity tests,
**not fitted predictions**. The current data justify no finite numerical
prediction of work per open-order solution. They certainly do not justify
advertising either target as reachable within 2^46.

## 4. A precise first Claude experiment and blind pass/kill gates

The recommended immediate action is to reproduce the proof controls and stop
unless a separate, predeclared calibration is wanted. No long open-order run is
justified by this report. Use the repository root, Python 3 with NumPy, and g++
with C++17 support. No new Python package or external solver is required.

```bash
python theory/openai_fresh_parametric.py --controls --out logs/openai_fresh/recheck_parametric.json
python theory/openai_fresh_cycle.py --controls --out logs/openai_fresh/recheck_cycle.json
python theory/openai_fresh_audit.py --benchmark --out logs/openai_fresh/recheck_audit.json
python theory/openai_fresh_verify.py logs/openai_fresh/sphere_Afirst_43.json logs/openai_fresh/sphere_Afirst_61_profile.json logs/openai_fresh/parametric_43_recorded.json --out logs/openai_fresh/recheck_matrices.json
```

[M/E] The retained controls took seconds; allow one minute for compilation and
machine variation. **Kill immediately on any optimum mismatch, invalid local
wall count, compressed-target mismatch, or failed integer Gram check.**

If a larger statistical test of the same outer iteration is desired, freeze the
algorithm, beta=1, quantization, check interval 10 and the following schedule
before running it. Do not tune on the holdout seeds or warm-start from known
signs. Accept any independently verified b=0 core, including one stored under
`outside_requested_sums`; sign/permutation normalization is allowed. Compression
must still be respected in any future 3452 experiment.

| Gate, in order | Blind setup | Trials / cap | Pass criterion |
|---|---|---|---|
| A | p=61,d=3, rows (11,11,1,1) | seeds 1001..1008, four trials each; 10,000 iterations/trial | at least 16/32 independently verified cores |
| B | p=131,d=5, rows (19,9,9,1) | seeds 2001..2004, four trials each; 30,000 iterations/trial | at least 8/16 cores |
| C | p=223,d=3, rows (11,11,11,23), and p=311,d=5, rows (21,19,21,1) | each: seeds 3001..3002, four trials each; 100,000 iterations/trial | at least 2/8 cores at **each** published-order control |

The p=131 profile is supported by `logs/w1_v131_d5.jsonl`; it differs from the
unverified (21,9,1,1) exploratory profile. Gate C uses only the published group,
multiplier and row data, never the published signs. Those originals are for
verification/provenance, not solver initialization.

Run one process per core, at most four. Give each process 5,400 seconds and the
whole calibration 12 CPU-hours. A trial interrupted by either cap counts as a
failure for this decision, rather than being silently removed. Stop at the
first failed gate; no substitution of profiles, seeds or beta after seeing
results. [E] The present evidence makes passing A uncertain and a full pass
unlikely; the thresholds are deliberately demanding requirements for scaling,
not success rates inferred from four exploratory trials.

For example, the first holdout process is exactly:

```bash
python theory/openai_fresh_cycle.py --kind b0 --p 61 --d 3 --profile 11,11,1,1 --parametric --spectral-concur --beta 1 --seed 1001 --trials 4 --iters 10000 --seconds 5400 --check 10 --out logs/openai_fresh/holdout_61_1001.json
python theory/openai_fresh_verify.py logs/openai_fresh/holdout_61_1001.json --out logs/openai_fresh/holdout_61_1001_verified.json
```

[E] Expect roughly two minutes for that four-trial p=61 process on this machine;
verify timing locally. Repeat with the exact seeds in the table and count trial
successes, not matrices or duplicate normalized outputs within a trial.

Only if A, B and C all pass is a **bounded exploratory pilot**, not a
budget-certified production search, warranted at 4412: four one-trial processes,
seeds 4001..4004, each at most 10,000 iterations and 5,400 seconds, with
`--p 1103 --d 19 --profile 1,37,39,39`; all other flags match the example.
Stop after these four trials. A failed pilot gives no production authorization
from this proposal; a success requires the independent full Gram check before
any matrix claim. Do not extrapolate cyclic b=0 success rates to 3452: a blind
known k4 control in the same group/fibre formulation is still missing. The 3452
code and timings are an oracle benchmark, not a completed calibration ladder.

## 5. What can actually be reconstructed about the August methods

### 5.1 The five unbordered examples

[M] `openai_fresh_audit.py` reads the originals directly, checks all core PAFs
with integer arithmetic, and computes flatness without relying on float FFTs.
For a block of length p, row sum r and autocorrelations C(t),

    Phi = [p sum_(t!=0) C(t)^2 − (r^2−p)^2] / [(p−1)p^2].

This is [P] Parseval applied to the definition in EKLO, with the principal
character omitted. Numerators and denominators are retained in `audit.json`.

| p | Four measured Phi values, rounded |
|---|---|
| 223 | 0.08983, 0.08983, 0.08983, 0.09490 |
| 283 | 0.08446, 0.08927, 0.08206, 0.08170 |
| 311 | 0.10472, 0.09358, 0.10472, 0.09627 |
| 487 | 0.13556, 0.13629, 0.12340, 0.11881 |
| 491 | 0.12529, 0.11375, 0.12620, 0.09566 |

Together with the repository's near-local-minimum fingerprints, this supports
[E] a procedure which finds individually flat blocks, then couples or joins
them. It does **not** recover the historical algorithm, its join procedure, or
its running time. The gap report's strong wording that it was the same
flat-pool k-tree family is an interpretation, not an exact computational result.
At p=283 that report itself estimates lists near 2^40 for a plain four-list
k-tree, while estimating the entire relevant flat set near 2^33 with large
uncertainty. This tension prevents a confident reconstruction. Bigger flat
pools would also repeat an existing repository route.

### 5.2 Orders 668 and 716

Set m=2h, x=(u,v), c=(u+v)/2 and d=(u−v)/2 on the first h coordinates.
[P] c,d are disjoint ternary sequences and

    sum_b PAF(c_b) = (2h−1) at 0, −2 at other lags;
    sum_b NAF(d_b) = (2h+1) at 0, 0 at other lags,

where NAF includes a minus sign on wraparound. These are exactly the periodic
and negacyclic components of the length-2h target, not two unrelated heuristics.
Once sigma=u*v fixes the supports, the two sign problems are independent.

[M] The originals at h=83 and 89 satisfy both identities at every lag. Each
block's sigma is an anchored two-step arc with one or two exceptions; there are
five exceptions total in each quadruple. The exact boundaries and exceptional
indices are in `audit.json`; they reproduce the earlier forensic finding.

[E] A support-template selection followed by independent complementary-sequence
searches is therefore a credible structural reconstruction. It does not reveal
whether the historical sign solver used end-to-middle search, a join or another
method. The calibrated counts in `report/alpoge_gap_analysis.md` turn negative
around h corresponding to m near 400, well below m=598,718,738,750. Extrapolating
those fixed templates to the open orders has no supporting cost model.

### 5.3 The proposed p=599 multiplier shortcuts are impossible

[P/M] For a prime p and a multiplier of odd order d, an invariant block has
row sum e+d(2a−(p−1)/d), e in {−1,+1}. A b=0 core requires the four squared row
sums to total 4p. At p=599, the allowed absolute sums below sqrt(2396) are
{1,25,27} for d=13 and {1,45,47} for d=23. Exhaustive four-square enumeration
finds no profile (`audit.json`). This is a fresh verification of an existing
obstruction, not a fresh empty-slice search.

[P] Allowing these subgroups to permute and sign-change the four blocks does
not rescue them: the signed permutation group has order 2^4*4!=384, with no
13- or 23-torsion. Any homomorphism from a subgroup of order 13 or 23 is trivial.
The same row obstruction therefore applies to such an equivariant extension
containing either subgroup. This excludes only this particular class of
four-block signed-permutation actions, not every construction at 2396.

## 6. A separate no-go lemma: low-rank contraction

**[P] Lemma.** Delete s rows and s columns from a Hadamard matrix H of order
N+s, leaving an N-by-N matrix A. If A+E is a Hadamard matrix of order N, then
rank(E) >= N−s.

**Proof.** Let C be the s-by-N portion in the deleted rows and retained columns.
Orthogonality gives A^T A=(N+s)I−C^T C. If v lies in ker(E) intersect ker(C),
then Av=(A+E)v. The first identity says its squared norm is (N+s)||v||^2;
orthogonality of A+E says it is N||v||^2. Hence v=0. C is therefore injective
on ker(E), so dim ker(E)<=s and rank(E)>=N−s. No symmetry or entrywise
assumption on E was needed beyond A+E being a Hadamard matrix. QED.

For example, contracting order 2400 to 2396 by deleting four rows and columns
requires a correction of rank at least 2392. A bounded-rank border repair
cannot do it. This does not forbid dense full-rank switching, and no claim of
historical priority for the lemma is made.

## 7. Ideas considered and rejected

- Fixed joint Fourier spheres alone: already published; short control success
  falls sharply with size. Replacing moving per-block targets is not new enough.
- Endpoint-only divide-and-concur: weak local constraints; no useful retained
  scaling evidence. Replaced by whole translation cycles, not ranked separately.
- Larger flat pools plus Wagner/k-tree joins: already developed here; the memory
  and pool-size obstruction is not removed by renaming the join.
- Half-shift arcs with a few exceptions at open lengths: calibrated support
  counts are unfavorable; re-running fixed-template searches is not justified.
- Reopening p=599 multiplier 13/23 or a signed four-block version: exact row-sum
  obstruction, including the equivariant extension proved above.
- Smaller staged k4 or block-exchange slices: already decided empty; no new
  ansatz was found that changes those conclusions.
- Tensor-Paley, three-stratum, elliptic or cyclotomic identities: already closed
  in `report/openai_4p3_review.md` and related reports.
- Ordinary quaternary-conference search at 646/940/1042/1108/1192: Kumar's route
  and its ingredient targets are established; no new competitive solver emerged
  from this audit, so repeating his construction is not ranked as a new idea.
- Low-rank repair of a nearby larger Hadamard matrix: ruled out by Section 6.
- Learned moves or generic SDP/lattice recovery: no specified training advantage,
  rounding probability or cost improvement was established; slogans are not proposals.

## 8. Files, reproduction scope and credit

The scripts are independent of the earlier S4 search, except for explicitly
loading the existing compression certificate and original data:

| Script under `theory/` | Role |
|---|---|
| `openai_fresh_sphere.py` | fixed Fourier-sphere baseline, group/orbit definitions and integer core checks |
| `openai_fresh_concur.py` | unsuccessful endpoint-replica baseline |
| `openai_fresh_cycle.py` | counted cycle DP, replica model and outer search |
| `openai_fresh_parametric.py` | faster integer cycle optimizer and exhaustive controls |
| `openai_fresh_audit.py` | original-data audit, row-sum gates and open-host projection benchmarks |
| `openai_fresh_verify.py` | independent GS construction and exact full matrix Gram check |

All experimental logs are under `logs/openai_fresh/`. Early sphere logs without
a `direction` field used B-first. Search logs saying “direct integer PAF only”
retain their original certification scope; the subsequent separate
`independent_matrices.json` adds full Gram verification for the eleven specified
small outputs. No search was run long enough to justify an open-order
nonexistence claim.

Primary sources and attribution:

1. Egan, Kotsireas, London and Ó Catháin (EKLO),
   [survey v2](https://arxiv.org/abs/2610.27175v2): arrays, Theorem 7 and the
   flatness analysis. The current open-order membership is taken from this
   repository's JSON, not inferred from a changing survey count.
2. Alpöge, Voinov and Reynolds-Haertle: the twelve announced matrices, with the
   original files preserved in `forensics/originals/`. This report audits their
   outputs; it does not attribute our hypothetical procedure to them as fact.
3. Kishan Kumar, [hadamard](https://github.com/vatskishan03/hadamard), Apache-2.0:
   constructions, ingredient routes, `FRONTIER.md`, `search/NOTES.md`, and the
   primary descriptions under `hadamard/2572`, `4204`, `6796`. Their staged,
   support-biased and quaternary-conference successes motivate looking beyond
   exhaustive slice searches. No Kumar implementation was copied here.
4. F. J. Aragón Artacho, R. Campoy, I. S. Kotsireas and M. K. Tam,
   [A feasibility approach for constructing combinatorial designs of circulant
   type](https://arxiv.org/abs/1711.02502), especially Proposition 3.9:
   the joint Fourier projection is prior work.
5. V. Elser, [The complexity of bit retrieval](https://arxiv.org/abs/1601.03428):
   binary/Fourier feasibility context and reflection methods.
6. S. Gravel and V. Elser,
   [Divide and concur: A general approach to constraint satisfaction](https://arxiv.org/abs/0801.0222):
   the replica/concurrence framework is prior work.
7. This repository: existing constructions, staged-search successes and exact
   closures, especially `report/alpoge_methods_forensics.md`,
   `report/alpoge_gap_analysis.md`, `report/open_screen_answer.md`,
   `report/h4892_answer.md`, `report/h863_lift_answer.md`,
   `report/blockswap_answer.md`, `report/openai_4p3_review.md`, and
   `report/tmat_routes.md`.

The unresolved issue is a whole-core search with a measured, adequate success
rate. The local optimization theorem and its implementation are reusable;
the present outer iteration has not supplied that missing result.
