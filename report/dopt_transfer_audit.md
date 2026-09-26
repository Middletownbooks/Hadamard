# D-optimal designs as a parallel track: transfer audit

Status legend: **[M]** measured here (script + output in `dopt/`), **[E]** estimate/model, **[R]** recalled or taken from
search-result snippets only (primary sources could not be fetched: arxiv.org, springer.com, adsabs and the maxdet
tables are blocked by this environment's egress proxy).
Predictions were registered in `dopt/PREREGISTRATION.md` before running; outputs are in `dopt/verify_output.json` and
`dopt/verify_p8.json`, produced by `dopt/verify.py`.

---

## 0. Verdicts at a glance

| # | question | verdict |
|---|---|---|
| 1 | §3.5 gate: do the open n ≡ 2 (mod 4) orders sit where the two-circulant route is live? | **YES → audit continues** [R] |
| 2 | §1 verified? | **Yes, with one registered FAIL (P1 counts off by exactly 2×, the (A,B)↔(B,A) swap)** [M] |
| 3 | §3.2: does the orbit reduction transfer? | **Yes for the current frontier (v ≲ 300); the k=2 coefficient crosses zero at v ≈ 326** [E, model checked against your MC] |
| 4 | module-by-module | table in §4 below |
| 5 | §5 positive control | **The join half already exists and passed (exhaustive, v = 7, 9, 11). Porting the local-search engine needs four named edits: an afternoon** [M for the join; E for the engine] |

Corrections to the prompt, in order of weight: §3.3 (two lists is *better* per unit of work, not worse);
§4 "spectral filter does not transfer" (it does, exactly: P8); §3.4 bites immediately (v = 99 is composite);
the "independent arrival at 90.9" is the same Gaussian model, not an independent one.

---

## 1. §3.5: the gate

- [R] Đoković–Kotsireas (2014/15, *D-optimal matrices of orders 118, 138, 150, 154 and 174*) built all five orders from
  two-circulant supplementary difference sets. Their abstract, as quoted in search snippets, says that D-optimal designs
  of order 2v < 200 (v odd) exist for every feasible order except v = 99, so **n = 198 is the smallest undecided order**.
  "Feasible" means 2v − 1 is a sum of two squares, which is exactly the two-circulant row-sum condition.
- n = 198: 2v − 1 = 197 = 14² + 1². The order is feasible and open, and the two-circulant (SDS) route is how every
  neighbouring order was closed. **Gate: PASS.**
- I could not confirm whether n = 198, or any order above 200, was settled after 2015. The table itself was
  unreachable. Treat "v = 99 is open" as [R], dated 2015.
- One non-two-circulant precedent: order 150 was first done by Holzmann–Kharaghani by another construction [R], so
  the two-circulant route is the main one but not the only one.

## 2. §1 checked

All by exhaustion, with every object rebuilt and checked in exact integer arithmetic:
R = [[A, B], [−Bᵀ, Aᵀ]], R·Rᵀ = diag(L, L), L = (2v−2)I + 2J, and |det R| = (2n−2)(n−2)^((n−2)/2) via fraction-free
Bareiss.

| check | result |
|---|---|
| P1 counts 392 / 1944 / 0 at v = 7 / 9 / 11 | **FAIL as registered.** Ordered pairs: **784 / 3888 / 0**, exactly 2× for v = 7 and 9. Since 4v−2 ≡ 2 (mod 4), \|a\| ≠ \|b\| always, so each unordered design is counted twice. The prompt's figures are consistent with unordered pairs. The v = 11 null matches. |
| P2 every pair is a true D-optimal design (Gram and Ehlich–Wojtas determinant) | PASS (784 + 3888 objects) |
| P3 negative tests: one flipped entry fails both checks; a wrong EW constant fails; the true pair passes | PASS |
| P4 PAF(s) ≡ v (mod 4) and Σ_{s≥1} PAF(s) = a² − v, all ±1 sequences, odd v ≤ 13 | PASS (0 counterexamples) |
| P5 smallest non-zero Σ_{s=1..(v−1)/2}(PAF_A + PAF_B − 2)² on the admissible manifold | PASS: 32 at v = 7 and v = 9 |
| P6 4v − 2 a sum of two squares among {487, 523, 571, 659, 739, 643} only for 643 | PASS. On the full 13-order worklist: **537, 573, 643** |
| P8 PSD_A(χ) + PSD_B(χ) = 2v − 2 at every non-trivial χ (v = 9, all pairs) | PASS (max deviation 1e−14; a non-pair deviates by 13.6) |

The Ehlich–Wojtas form and the necessary condition "2n − 2 is a sum of two squares" (n = 2v, so 4v − 2, equivalently
2v − 1) agree with the exhaustion at v = 11 (42 = 2·3·7, no designs) [M]. The attribution itself is [R].

## 3. §3.2: does the orbit reduction survive k = 2?

**Model.** For a sum of k PAF values with per-block variance ≈ v on a lattice of step 4:
P(sum = target) ≈ 4/√(2πkv), so **b(k, v) = ½·log₂(2πkv) − 2** bits per equation [E].

| (k, v) | your Monte Carlo b | model b | difference |
|---|---|---|---|
| (2, 67) / (2, 79) / (2, 103) / (2, 151) | 2.86 / 3.00 / 3.19 / 3.47 | 2.86 / 2.98 / 3.17 / 3.45 | ≤ 0.025 |
| (4, 67) / (4, 79) / (4, 103) / (4, 151) | 3.50 / 3.56 / 3.76 / 4.02 | 3.36 / 3.48 / 3.67 / 3.95 | ≤ 0.14 |
| (4, 523) | (needed 4.87 to reproduce 90.9) | **4.84** | 0.03 |

(P7: PASS.)

**Where the coefficient k − b/2 crosses zero:** b = 2k gives 2πkv = 2^(4k+4), so **v\* ≈ 326 for k = 2** and
≈ 41,700 for k = 4 [E].

**Where I think the prompt's reading needs correcting:**

1. **The live frontier sits below the crossing.** The first open order is v = 99, and the next ones are near
   100–150. There the k = 2 coefficient is +0.6 to +0.25 per orbit, which is small but positive. So the orbit
   reduction is *available* on today's D-optimal frontier. It becomes "expected empty" only past v ≈ 300, beyond the
   current frontier. The verdict is "engine yes, reduction yes but thin", not "reduction no".
2. **The 90.9 is not an independent arrival.** My estimator (`nbits.py`) and your Monte Carlo both measure the density
   of a lattice-Gaussian sum of k PAFs. Agreement confirms the arithmetic, not the model's physics.
3. **The model can be tens of bits off on structured instances, in the pessimistic direction.** Measured on this
   project:
   - v = 347: the independence estimate gives 2^−89 lift solutions for the real compression; the exact join finds 6.
   - The two-sided z-mod-3 stage on 347 predicts 2^−42.8 tuples; 7 are found.

   Near the crossing (|log₂E| of a few bits) the Gaussian model cannot decide existence either way. It should rank
   configurations, not kill them.
4. **§3.4 bites at the very first open order.** v = 99 = 9·11 is composite, so ⟨t⟩ does not act freely and
   r = (v−1)/d + 1 does not hold. The machinery that handles composite index groups is what this project built for
   the b ≥ 3 route: compression onto quotients (Z₉, Z₁₁, Z₃₃), exact stage A on the quotient, residue filters, and
   the two-sided lift for a product of two coprime factors. v = 99 = 9·11 is structurally the same situation as
   h = 85 = 5·17. I recall Đoković–Kotsireas using compression for composite v [R].

## 4. Module-by-module verdict

The transfer edits come in three kinds:
- **T1:** replace the target with 2 per lag class, i.e. the residual becomes Σ_k PAF_k − 2.
- **T2:** k = 4 → 2, which also drops the per-slot profile logic.
- **T3:** the row-sum set becomes a² + b² = 4v − 2.

| module | verdict |
|---|---|
| `redlib.py` (orbits, lag and frequency classes, prime v) | **as-is** for prime v; does not handle composite v (§3.4) |
| `hadlib.py` (PAF/PSD, multiplier stabilisers) | **as-is** |
| `engine.py` (W2 local search: incremental ΔPAF under a flip, S2/S3 modes, tabu) | **named edits:** T1 in the residual R_c; `random_y` hard-codes 4 blocks (T2); S2 mode divides the residual by 4 → by k; initial row-sum profile (T3). The ΔPAF algebra is per-block and unchanged. The joint modes use single-orbit flips; row-sum-preserving swaps exist in `flat_tabu` (single block) and carry over unchanged. |
| `nbits.py` (Gaussian n) | **named edit:** T2 + T3 (sum over signed (a, b) with a² + b² = 4v−2); T1 changes nothing in the Gaussian |
| `kt_gen.py` (flat-pool tabu generator) | **as-is** (single-block), with the row-sum targets from T3 |
| `kt_join.py`, `kt_run.py`, `kt_generic.py`, `kaggle/kt*_*.py` (4-list Wagner join) | **do not transfer, and are not needed.** With two lists a solution is an exact lookup (PAF_B = 2 − PAF_A): a hash join finds *all* M² pairs in O(M log M). `dopt/verify.py` already does this. |
| `kt_validate.py`, `rho_star.py`, `s1_*` (flat-pool density, join balance) | **named edit** T2 (a pair balance instead of a quadruple); the Wang–Landau density code is single-block and as-is |
| `w1_null.py`, `w1_analyse.py`, `dirichlet_null.py` (fingerprint nulls) | **named edit:** the per-class share statistic becomes Dirichlet(a, a) on 2 blocks, and the PSD sum is 2v−2 (P8) |
| `s4_rrr.py` (RRR / divide-and-concur with spectral sets) | **named edit:** the spectral projection uses PSD_A + PSD_B = 2v − 2 (P8), with 2 blocks |
| `plot_scaling.py`, growth-law fit, `grid.py` | **as-is** (protocol) |
| `b1*.py`, `kt_b1.py` (b = 1 bordered GS cores) | **do not transfer** (Hadamard-specific border) |
| `bge3_border`, `bge3_skeleton*`, `bge3_template347`, `bge2_*`, `b347/359/479_aut`, `census`, `extract_*` | **do not transfer** (reverse engineering of published matrices) |
| `bge3_cost.py`, `bge3_screen*.py` (cost model and exact quotient screen over Q × Z_h) | **named edits:** T1–T3, and drop Q (c = 1). The exact quotient-compression screen is directly useful for composite v (v = 99 onto Z₉, Z₁₁, Z₃₃). |
| `bge3_stageA.py`, `bge3_stageB_join.py`, `bge3_liftenum.py`, `bge3_residue.py`, `bge3_joint_residue.py`, `bge3_twosided*.py` (compression and lift, exact join, forced-residue filters, two-sided lift) | **named edits:** T1–T3; Q → trivial; the 4-list joins become 2-list hash joins (simpler). This is the most valuable transfer: it is the composite-v machinery (§3.4). Validated end to end on 347, including recovering the real core. |
| `kaggle/kt523_B2_*.py` (notebook harness: instance guard, `/tmp` caches, solutions printed to the cell) | **harness as-is**; payload replaced by the 2-list join |

## 5. §5 costed: positive control

- **Join half: done [M].** `dopt/verify.py` is an exhaustive two-list hash join over full space. It finds all 784
  (v = 7) and 3888 (v = 9) designs and proves v = 11 empty. Every object is rebuilt and checked (Gram + exact
  determinant), with negative tests that fail as they should. Full-space exhaustion is 2^v per list, so the same code
  reaches v ≈ 25–31 in minutes to an hour. That gives a graded ladder of known answers (v = 13, 15, 19, 21, 23, 25, …
  wherever 2v − 1 is a sum of two squares).
- **Engine half: an afternoon [E].** The four named edits to `engine.py`, then seven seeds at matched caps on the
  ladder above, with the prime-v orbit ansatz at d = 1 and 3. Success criterion, to be registered before running:
  the median seed finds a design, verified as an object, at every rung up to the first cap failure.
- **The growth-law re-fit (§5.2)** comes out of that same ladder run. I agree it is the most valuable output for the
  Hadamard track. Caveat: compare exponents in the same units (seconds against −log₂ density in constraint bits),
  since the k = 2 and k = 4 density scales differ.

## 6. Where the prompt is wrong or incomplete (summary)

1. **P1 counts are unordered.** The registered prediction fails; ordered counts are 2×.
2. **§3.3 is backwards on cost.** The two-list meet-in-the-middle covers M² pairs at *linear* cost, which is optimal;
   the four-list k-tree covers M⁴ combinations at super-linear cost. Hit counts scale as M² against M⁴, so matching
   the Hadamard hit rate needs pools of ≈ M², which is a memory bound, not a new move set.
3. **§4, spectral filtering does transfer exactly.** J only touches the trivial character, so PSD_A + PSD_B = 2v − 2
   at every other character (P8) and the per-block bound is PSD ≤ 2v − 2.
4. **§3.2:** the crossing is at v ≈ 326 for k = 2, above the current frontier. The 90.9 agreement is the same model,
   not independent. And the model is known to be pessimistic by tens of bits on structured instances.
5. **§3.4:** not hypothetical. The first open order, v = 99, is composite, and the composite-group compression and
   two-sided lift machinery transfers.
6. **Hadamard worklist ∩ D-optimal feasibility:** 537, 573 and 643 satisfy the two-squares condition, but the
   D-optimal orders 2v at those v are far past the D-optimal frontier (v ≈ 99). Instances are not shared, as §3.1
   says.
