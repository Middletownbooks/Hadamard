# Williamson matrices of order 65: pre-registration

Written 2026-09-30, before any Williamson code exists in this repo. Commit time stamps it.

## Question

Do Williamson matrices of order n = 65 exist? That is, four symmetric ±1 sequences A, B, C, D on Z_65 with
PAF_A(s) + PAF_B(s) + PAF_C(s) + PAF_D(s) = 0 for every s ≠ 0.

Status going in: exhaustive searches cover every odd n ≤ 59 (none at 35, 47, 53, 59). Turyn's construction gives
61 and 63. Kharaghani–Mohammadian–Tayfeh-Rezaie–Zaitsev (arXiv 2605.08661, v2 Sept 2026) list 65 among orders
where Williamson matrices "have not yet been found". 65 is the smallest undecided odd order.

## Method (fixed now; changes after this point are logged as amendments)

n = p·q with gcd(p, q) = 1; for 65, p = 5 and q = 13. Index i ↔ (i mod p, i mod q), so a block is a p×q ±1 array
M with M[−a,−b] = M[a,b].

**Normalisation.** Negating a block or permuting blocks preserves the property, so every Williamson quadruple is
equivalent to one with row sums s_A ≥ s_B ≥ s_C ≥ s_D > 0. N(n) is the number of such normalised quadruples
(labelled; no further reduction). Σ s_X² = 4n.

1. **Stage A, compressions.** For m ∈ {p, q} enumerate all symmetric compressed sequences (entries of parity
   n/m, |entry| ≤ n/m, given sum) whose DFT satisfies |X̂(k)|² ≤ 4n at every nonzero k. Join them into
   quadruples with Σ PAF = 0 off zero. Exact: no sampling, no caps.
2. **Stage B, consistency.** Combine a Z_p quadruple and a Z_q quadruple with the same row sums; keep a
   combination only if every block's margin pair (row sums, column sums) has at least one symmetric ±1 lift.
   Reduce combinations modulo the multiplier group (Z_n)* acting on all four blocks together (exact: a
   multiplier maps solutions to solutions).
3. **Stage C, lifts.** For each block and margin pair, enumerate every symmetric ±1 array with those margins,
   keeping those with |X̂(k)|² ≤ 4n at all k ≠ 0.
4. **Stage D, join.** Meet in the middle on (A,B) against (C,D). The key is the PAF vector together with the
   product vector a_j b_j (j ≠ 0). This uses the Williamson product condition: a_j b_j c_j d_j = ε for all
   j ≠ 0, with ε fixed by the row sums (the parity identity checked on 1,000 random symmetric sequences on
   2026-09-30).

Every hit is checked by an independent verifier that computes PAFs directly from the ±1 sequences.

## Predictions (single threshold each; anything that cannot be evaluated counts as FAIL)

- **W0, build validation.** For n ∈ {15, 21}, the pipeline's N(n) equals N(n) from an independent brute force.
  The brute force uses no compression, no multiplier reduction and no product condition. Both counts come
  from the same normalisation. The pipeline must also give the same N(n) with the multiplier reduction on and
  off. Every brute-force solution must satisfy the product condition, and its compressions must appear in the
  pipeline's stage-A list.
- **W1, negative control.** n = 35 (= 5·7; nonexistence known, Đoković 1993): the pipeline reports 0 solutions
  with every stage exact.
- **W2, positive controls.** n ∈ {45, 51, 55, 57, 63}, all composite with a coprime split. Williamson matrices
  exist at each by Turyn (q = 89, 101, 109, 113, 125). The pipeline finds and verifies at least one quadruple
  at each of the five. A miss at any one is FAIL.
- **W3, cost gate for 65.** Before the stage-D join at 65, compute its exact work, Σ over surviving
  combinations of |L_A||L_B| + |L_C||L_D|. Convert it to core-hours with the rate measured in W2. Run 65 only if
  the projection is ≤ 500 core-hours (this container plus Kaggle). Otherwise stop and report the number.
- **W4, reading the result at 65.**
  - **Found:** a quadruple that passes the independent verifier. It is reported regardless of W0–W2.
  - **Not found:** a claim of nonexistence is made only if all of the following hold. Otherwise the result is
    "no solution found in an incomplete search".
    - W0, W1 and W2 all PASS.
    - No stage hit a cap or truncation.
    - The stage-A quadruple lists at 65 match a second, independent stage-A implementation.
    - The completeness ledger shows a stage-D record for every canonical combination.

## What does not count

- Our Gaussian cost model. It predicts solutions at 35, 47, 53 and 59, where there are none, so it is not
  used for any decision here.
- Any sampled stage.

## Amendments

**A1 (2026-09-30, before any W2 result was read): compiled stage-D join.** Same algorithm and same keys as
`join4`: a meet in the middle on the PAF sum plus the product-condition mask. It is now in numba
(`scripts/wil_join_nb.py`), with a 64-bit linear PAF hash so each pair's key costs O(1). Every hit is still
re-verified exactly.
- It gives the same solution sets as the Python join at n = 15, 21 and 35.
- The W2 controls are rerun with it, so the W3 rate comes from the join that would run at 65.
- A partial n = 45 run with the Python join (case 1: 2 solutions, 304 s) is kept as
  `logs/wil/n45_pyjoin_partial.log`. The compiled join gives the same 2 solutions in 15.5 s.

**Note, no change to method:** the mod-4 compression condition implied by the product condition (the four blocks'
coset sums ≡ fixed residue mod 4) holds for every complementary compression quadruple at n = 35, 45 and 65. It is
implied by complementarity, so it has no filtering power and is not used.

**A2 (2026-09-30, before any n = 55, 57 or 63 result was read): pair-PSD filter in the compiled join.** A pair
(A,B) with PSD_A(f) + PSD_B(f) > 4n at some f ≠ 0 cannot be part of a solution, because the four PSDs are
nonnegative and sum to 4n. Such pairs are skipped on both sides.
- **Filter power at 65:** it passes about 3.4% of pairs (12 random combinations from case (9,9,7,7)).
- **Implementation:**
  - each side's partner list is sorted by PSD at every frequency;
  - for each element, the frequency admitting the fewest partners is chosen, and only that prefix is scanned;
  - the hash table is sized to the surviving pairs.
- **Validation:** solution sets are identical to the Python join at n = 15, 21 and 35. The PSD-filter-only
  variant also gave identical sets at n = 45.
- **Rate at n = 65 sizes:** about 7.5×10^7 pair operations per second steady-state with a machine load of about 3,
  against 1.6×10^7 for A1 under the same load. The work unit is unchanged: |L_A||L_B| + |L_C||L_D|, counted before
  the filter.
- **W2 rerun:** all five controls are rerun on this kernel. The A1 results (n = 45: 2 representatives and 12
  quadruples; n = 51: 4 representatives and 32 quadruples) are kept as `logs/wil_n{45,51}_A1.json`.

**Sharding (implementation only):** `--case c` and `--shard i/K`. Shard i takes the Z_p representatives with index
≡ i mod K. The shards partition the ledger exactly (checked at n = 21 and 35).

**Implementation note (no change to method):** the lift cache is now bounded (`WIL_CACHE_GB`, default 2) and is
cleared when it passes the limit. Lifts are recomputed on demand, so results are unaffected. Checked at n = 21 and
35 with forced clears: identical solution sets. The unbounded cache reached about 7 GB at n = 55 and 57 and caused
an OOM kill, so n = 55, 57 and 63 were restarted from scratch with the bound. The killed partial logs are kept.

**A3 (2026-09-30, decided by the user on cost grounds): W2 at n = 63 is replaced by a targeted check.** An
exhaustive n = 63 run is about 8×10^13 pair operations for its first row-sum case alone (8.1M combinations),
and on the order of 2,000 core-hours overall, more than n = 65 itself. At amendment time the local n = 63 run
had covered about 2% of its first case with no solution (seen before this amendment; the log is kept). That
partial null carries no weight either way.
- **W2 (amended):** at least one verified quadruple at each of n = 45, 51, 55 and 57, found blind. (Status when
  amended: all four had passed.)
- **W2-T63 (new, targeted, not blind):** build a Williamson quadruple of order 63 independently (Turyn's family,
  q = 125) and check it with the independent verifier. Then show that the pipeline:
  - (a) keeps each of its Z_7 and Z_9 compression quadruples in stage A;
  - (b) keeps its canonical combination (after multiplier reduction) with every block liftable;
  - (c) returns a verified solution when stage D is run on that one combination, and that solution is in the
    multiplier orbit of the constructed one.

  PASS iff (a), (b) and (c) all hold. A construction that fails the verifier makes W2-T63 unevaluable, which
  counts as FAIL.
- The W4 condition "W0, W1 and W2 all PASS" now reads "W0, W1, W2 (amended) and W2-T63 all PASS".

## Addendum S (2026-10-06, after the W(65) null; registered before any SAT run): independent SAT cross-check

**Encoding** (`scripts/wil_sat.py`, new code that shares nothing with `wil_search.py`):
- **Variables:** one Boolean per half-position of each symmetric block (x = −1). That is 33 per block at n = 65.
- **Row-sum case:** fixed per instance (the same 5 cases, re-enumerated independently). The middle entry's parity
  is forced, and the half-weight is an exact cardinality constraint.
- **Williamson condition:** for every lag s, the four blocks together have exactly 2n disagreements x_i ≠ x_{i+s}
  over i ∈ Z_n. This is encoded with XOR auxiliaries and an exact cardinality constraint (totalizer), and is
  equivalent to Σ_X PAF_X(s) = 0.
- **Not used:** compressions to Z_5, lifts, PSD pruning, the product condition, hashing and the Z_5 multiplier
  reduction.
- **Splitting (only if needed for speed):** instances may be split by the Z_13 compression quadruple, taken from
  the independent stage-A implementation (`wil_stageA_alt.py`). Canonicalisation under (Z_13)*/±1 is reimplemented
  inside `wil_sat.py`.

**Predictions:**
- **S0 (validation):** by full enumeration with blocking clauses, SAT gives the same number of normalised quadruples
  as the pipeline:
  - N(15) = 24 and N(21) = 96, the brute-force values;
  - 0 at n = 35;
  - 12 at n = 45 (the pipeline's orbit-expanded count). This one is evaluated if it finishes within 6 core-hours,
    otherwise it is reported as unevaluated.
- **S1 (cross-check):** at n = 65 every instance returns UNSAT, with no instance left unsolved. A SAT result is
  checked by the independent verifier and would contradict W(65). An instance that does not finish leaves S1
  unevaluated, and that is reported.

### Amendment A-S1 (2026-10-06, after S0 at n ≤ 21 and before any n = 35, 45 or 65 result of the final encoding)

**Why:** without extra relations the plain encoding took about 10.6 s per compression combo at n = 35 and over 20 min
per Z_13 instance at n = 65, so S1 would be infeasible. On the user's instruction, the checker now encodes every
known necessary Williamson relation, not only the PAF condition. Each is a classical theorem, and each is coded
from its statement inside `wil_sat.py` without reusing code from `wil_search.py`.

**Added relations:**
- **Product condition** (Williamson):
  - a_j b_j c_j d_j = ε for j ≠ 0, with ε = (−1)^(W−n) and W the total number of −1 entries;
  - the opposite sign at j = 0.
  - These are encoded as 4-literal parity clauses. They were checked on all brute-force solutions at n = 15, 21
    and 25 before use.
- **Redundant pair-count form:** the number of positions where two blocks are both −1, encoded as an exact
  cardinality constraint (W − n) over AND auxiliaries.
- **PSD bound:**
  - PSD_X(k) ≤ 4n for every block X;
  - PSD_X(k) + PSD_Y(k) ≤ 4n for every pair of blocks;
  - both over all k.
  - It is enforced by an IPASIR-UP external propagator (CaDiCaL 1.9.5). The propagator only adds blocking
    clauses for complete blocks that violate the bound, so it prunes and never asserts.
- **Splitting:** by the pair (Z_5 compression, Z_13 compression), each canonicalised under its own multiplier
  group (the CRT action is independent on the two factors). The compressions are imposed as assumptions on unary
  totalizer counters.
  - This now uses the Z_5 reduction that Addendum S listed as unused. It is reimplemented, not shared.
  - Counts are therefore of reps under this reduction, and are compared with the pipeline's per-case rep counts.

**Validation already done with the final encoding:**
- At n = 21, per-case counts of 22, 6 and 4 equal the pipeline's.
- The S0 bars at 35 and 45 are unchanged.
- S1 is unchanged: every n = 65 combo is UNSAT, none unsolved, or else S1 is reported unevaluated with measured cost.
