# H(2588) (v = 647) in the 647 ×353 doubly-staged slice: exact cross-level solver, calibration, cost and run

Date: 2026-10-14. Branch `claude/focused-planck-j1645f`. Answers the 647 ×353 prompt (validate and calibrate the cross
level, cost it, run it, find a same-shape control).
Scripts: `theory/p647_*.py`, `scripts/verify_h2588.py`. Logs: `logs/p647/`.

**Tags.** [M] exact computation (script named); [P] proved here (proof in the text); [E] Gaussian estimate (estimator
and calibration named); [L] literature recollection, unchecked.

---

## 0. Bottom line and run status

*Status at 09:45 UTC 2026-10-14; this section is updated as the run proceeds.*

1. **(A) No H(2588) yet.** The production run is under way (§4); no core so far, against 0.98 expected.

2. **(B) The cross level is calibrated, and its cost per core is inside the budget** [M speeds × E counts, the counts
   exact]:
   - **Validated exact solver** (§1):
     - joint solution → per-layer lists → distinct σ (an additive base-3 code, P1) → exact S-join (a new
       radix-partitioned kernel) → D-check of every S-solution (the even-lag identity, P2) → full PAF;
     - plants: 12/12 at full size on real 647 joint solutions, 10/10 on 495;
     - against brute forces that never use the pipeline's signatures: 3/3 reduced plants and 5/5 reduced
       real-S-solution instances.
   - **Exact counts on 4 complete joint solutions** (§2):
     - 421,570–536,655 S-solutions each, so the distinct-σ Gaussian is +1.5–1.7 bits off;
     - pair-level yield (D-candidates) 2^28.0–2^29.1, so the pair-weighted Gaussian is 0.46–0.63 bit too high.
   - **The D-level holds to rank 9 of 11** (§2.3):
     - over all 2^30.65 candidates of the four runs: 44 hits against 44.0 expected at rank 7, 5 against 8.0 at
       rank 8, 1 against 1.15 at rank 9 (30.6 bits);
     - every 647 S-solution's candidates span the full D-lattice;
     - C(D|S) = 34.8 bits is used; the last ~4 bits are an extrapolation.
   - **Cost per core** (§3, §4):
     - 2^41.0 emitted pairs, about 15 core-hours per expected core, with energy-class selection (measured);
     - about 6–9 core-hours with class-merged units;
     - 2^42.6 (about 40 core-hours) without these levers.
   - **Positive control** (§5):
     - 495 ×83 is the first same-shape slice with a non-empty joint stage (y_B orbits with 2,400–6,912 joint
       solutions);
     - its S-level and pair-level counts reproduce the 647 calibration pattern;
     - its full D-level gives **0 cores at 7.4 expected**, and the cause is identified: rank-deficient
       per-S-solution D-lattices (median 24 candidates), the discrete-level failure mode;
     - this does not occur at 647 (median 840, full rank);
     - every smaller product slice tested is exactly empty.
     
     So no valid positive control exists at small orders for the 647 D-level, and the (B) claim rests on the exact
     calibration of §2 (to rank 9).

3. **(C) No obstruction.** S-solutions are abundant. The D-target lies in the D-lattice of every calibrated
   joint solution, and in every per-S-solution lattice tested at 647. The 647 deep tail matches the Gaussian
   wherever it can be counted.

**Run status** (`python3 theory/p647_tally.py`): see §4.3; cores found 0.

---

## 1. Tools and their validation [P, M]

### 1.1 Exact identities used by the solver [P]

Per block k and Z₂-layer e the labels of the 35 M-orbits of the layer form a 5 × 7 ±1 matrix with prescribed weighted
row sums (y_B) and column sums (y_A); `nd_cross.Cross.layer_labellings` lists them exactly (lists L0, L1).
With a₀, a₁ the two layers at the 35 orbit representatives of Z₃₂₃, σ = (a₀ + a₁)/2 and δ = (a₀ − a₁)/2.

- **P1 (additive σ code).** For a ∈ {±1}, (a₀ + a₁)/2 + 1 = (a₀ + 1)/2 + (a₁ + 1)/2. Hence
  code(σ) = Σ_b (σ_b + 1)·3^b = c₀(a₀) + c₁(a₁) with c_e(a) = Σ_b ((a_b + 1)/2)·3^b, and code is a bijection from
  {−1, 0, 1}³⁵ onto [0, 3³⁵) ⊂ [0, 2⁶³). So the distinct σ of a block are the distinct values of c₀[i] + c₁[j]: one
  int64 `unique` over n₀n₁ ≤ 2.1·10⁶ values, no row sorting.
- **P2 (D-level by even lags).** a₀ = σ + δ and a₁ = σ − δ give PAF(a₀) + PAF(a₁) = 2(PAF(σ) + PAF(δ)). So given
  an S-solution (Σ_k PAF(σ_k) = T_S), the D-equations Σ_k PAF(δ_k) = T_D are equivalent to
  Σ_k (PAF(a₀,k) + PAF(a₁,k)) = T_E = 2(T_S + T_D) = (2584, −4, …, −4), the even-lag PAF of x on Z₆₄₆. The D-check
  needs only the per-layer PAFs (n₀ + n₁ vectors per block), no per-pair D-signature.
- **P3 (block energy is fixed mod 24).** Write σ as a 5 × 7 matrix (rows {0, A₁..A₄}, columns {0, B₁..B₆}); its
  weighted margins are s_B(i) = σ_{i0} + 3Σ_{j≥1} σ_{ij} and s_A(j) = σ_{0j} + 4Σ_{i≥1} σ_{ij}, fixed by the joint
  solution. Since σ_{ij} ∈ {−1, 0, 1} and these are distinct mod 3 and mod 4, σ_{i0} ≡ s_B(i) (mod 3) and
  σ_{0j} ≡ s_A(j) (mod 4) fix the whole first row and column. Since σ ≡ σ² (mod 2), the number n₁₂ of nonzero cross
  entries satisfies n₁₂ ≡ Σ_{i≥1} (s_B(i) − σ_{i0})/3 (mod 2). With orbit sizes 1, 3, 4, 12,
  ‖σ‖² = σ₀₀² + 3n₃ + 4n₄ + 12n₁₂ is therefore fixed mod 24 per block. [M: every list of every calibration joint
  solution lies in one class mod 24; the S-solutions sit in left-energy classes 24 apart.] The join's energy split
  already respects this, so it costs nothing extra and saves nothing extra.

### 1.2 The exact S-join (`theory/p647_join.py`) [M]

- Four lists of distinct σ with 20-coordinate signatures (PAF at the 20 ⟨M, −1⟩-lag classes of Z₃₂₃; coordinate 0
  is the energy ‖σ‖²). Pairing (a, b | c, d).
- Every row gets a 64-bit additive hash and a residue mod a prime P = 4093 (both linear, so additive). For every left
  energy EL and every group of residues, the left and right pairs of that class are emitted (inner lists sorted by
  (energy, residue) with a CSR index), scattered into 2^rbits partitions by the top key bits, and each partition is
  joined in cache. The exact pair count of every (EL, residue) class is precomputed by a 2-D FFT, so passes have a
  known size. Packed 8-byte entries (24-bit fingerprint, 2 × 20-bit indices); a fingerprint match is confirmed on the
  64-bit key, a key match on the exact integer signature.
- **Self-test** (`python3 theory/p647_join.py --selftest`): 40/40 (16-byte kernel) and 40/40 (packed kernel) against a
  brute-force 2+2 dictionary join on random instances, and `nd_ejoin.ejoin4` agrees 40/40.
- **Speed** [M]: 22–26 ns per emitted pair on one energy class of a real joint solution (EL = 319 of j = 9902,
  1.05·10¹⁰ pairs, 151,950 S-solutions, identical hit sets for five kernel settings); about 30 ns under a 4-worker load.
  The per-x overhead and the scatter dominate; the 523 kernel ran at a similar 17–26 ns.
- The old `nd_ejoin.ejoin4` re-enumerates every pair in each of its up to 64 partition passes, which is why it could
  not run the cross level at this size.

### 1.3 The cross pipeline (`theory/p647_cross.py`) [M]

Lists → distinct σ (P1) → S-join → D-check of every S-solution (P2; a compiled 2+2 on 64-bit hashes over the product of
its σ-preimages, every match re-checked exactly) → x_k on Z₆₄₆ → full cyclic PAF check.

**Validation [M].**
- **Plants on real joint solutions** (`--plant`): one random (L0, L1) pair per block of a random joint solution of
  y_B #0; targets = the plant's own S-sums and even-lag sums; the S-join is run on the plant's own left-energy class
  (2.6·10⁸–1.6·10¹⁰ emitted pairs each). **12/12 recovered**, at the S-level and as the unique D-hit
  (`logs/p647/plants_el10.log`: 10/10, joint solutions 2800, 3095, 5761, 12389, 5800, 692, …; and 2/2 earlier).
- **Reduced plants against a brute force** (`--brute`): every kept pair → its sequence on Z₆₄₆ → its full 646-lag PAF
  by FFT; all quadruples with the plant's PAF sum by a 2+2 dictionary on the exact vector: 3/3 agree.
- **Real S-solutions against a brute force** (`--brute-real 4208 5`, `logs/p647/brute_real_j4208.log`): reduced lists
  around 25 real S-solutions plus 40 random σ per block; the pipeline join (real target) returns exactly the members
  of the complete S-solution set inside the lists, and so does an S-level brute force whose PAFs are computed from
  scratch (direct shifts on Z₃₂₃); for 3 of them per trial the D-check agrees with a full 646-lag brute force over all
  1.4k–3.2k candidates: **5/5**.
- Kernels: `p647_join.py --selftest` 40/40 (16-byte), 40/40 (packed), 40/40 (optimised), 30/30 (multi-left/right),
  `nd_ejoin.ejoin4` 40/40, all against a brute force; five kernel settings return identical hit sets on a real class.

---

## 2. Calibration of the cross level [M against E] (`theory/p647_calib.py`, `logs/p647/calib/`)

Four joint solutions of y_B #0 were solved exactly and completely (every left-energy class of the S-join, every
S-solution D-checked): two from the high-freedom class (H ≈ 79.4 bits) and two from the 78.3–78.4 class.

### 2.1 Three Gaussian predictions per joint solution [E] (`gauss_preds`, N = 20,000 sampled quadruples)

- **C_S^u**: S-level bits for σ drawn uniformly from the distinct σ of each block. Predicts the exact S-join count
  2^(H_S − C_S^u), H_S = Σ_k log₂ N_k.
- **C_S^p**: S-level bits for (L0, L1) pairs drawn uniformly (σ weighted by multiplicity). Predicts the pair-level
  S-solutions 2^(H − C_S^p), i.e. the D-candidates Σ_{S-sol} Π_k mult(σ_k), which is what the D-check sees.
- **C**: S and D together, pairs uniform (rank 23, Mahalanobis 0.48–0.53). Predicts 2^(H − C) cores, and
  C(D|S) = C − C_S^p.

### 2.2 Exact counts [M]

| j | H class | N_k (distinct σ) | H_S | S-solutions exact | 2^(H_S − C_S^u) | error | D-candidates exact | 2^(H − C_S^p) | error | D-hits |
|---|---|---|---|---|---|---|---|---|---|---|
| 4208 | 78.33 | 212k, 229k, 389k, 129k | 71.05 | 421,570 = 2^18.69 | 2^17.16 | **+1.52** | 2^28.00 | 2^28.53 | **−0.53** | 0 |
| 9902 | 78.40 | 212k, 227k, 274k, 205k | 71.20 | 465,177 = 2^18.83 | 2^17.35 | **+1.48** | 2^28.12 | 2^28.75 | **−0.63** | 0 |
| 11514 | 79.42 | 327k, 227k, 182k, 310k | 71.83 | 536,655 = 2^19.03 | 2^17.38 | **+1.66** | 2^29.07 | 2^29.55 | **−0.48** | 0 |
| 16846 | 79.39 | 212k, 227k, 389k, 206k | 71.71 | 525,492 = 2^19.00 | 2^17.32 | **+1.69** | 2^29.07 | 2^29.53 | **−0.46** | 0 |

- **The distinct-σ Gaussian undercounts S-solutions by 1.5–1.7 bits.** Its sampled quadruples have mean energy
  Σ_k ‖σ_k‖² ≈ 760 against the target 645 (Mahalanobis 1.8–2.0, in the tail), and at low energy the other PAF
  coordinates are less spread than the joint Gaussian assumes. Same sign and size as the 523 joins (+1.55–1.75).
- **The pair-weighted Gaussian, the one that sets the core rate, is 0.46–0.63 bit too high**, at Mahalanobis 0.5:
  under pair weighting the block energies average 160–168, so the target energy is typical.
- **S-solutions carry 2^9.2–2^10.1 D-candidates each on average** (2^12.4 at the lowest-energy classes, 2^9.1 at
  the centre), against a mean multiplicity product of 2^7.7 over the lists: the S-level energy equation selects
  low-energy, high-multiplicity σ.

### 2.3 The D-level [M against E] (`dtail`)

Every D-candidate of every S-solution (2^28.0–2^29.1 per joint solution) was enumerated. For subsets of the 12
cross lags, the exact fraction of candidates whose D-residual vanishes on the subset is compared with the lattice
Gaussian fitted to 200,000 sampled residuals of the same candidates:

| cross lags (Z₃₂₃) | j = 4208: exact / Gaussian bits | j = 9902 | j = 11514 | j = 16846 |
|---|---|---|---|---|
| 1 | 3.89 / 3.86 | 3.95 / 3.94 | 3.93 / 3.92 | 3.91 / 3.89 |
| 1, 2 | 7.74 / 7.70 | 7.83 / 7.81 | 7.79 / 7.76 | 7.75 / 7.70 |
| 1, 2, 3 | 11.55 / 11.50 | 11.61 / 11.60 | 11.59 / 11.53 | 11.49 / 11.44 |
| 1, 2, 3, 4 | 15.16 / 15.11 | 15.11 / 15.08 | 15.15 / 15.07 | 15.10 / 15.02 |
| 6, 7, 8, 9 | 14.87 / 14.74 | 14.86 / 14.78 | 14.87 / 14.75 | 14.85 / 14.75 |
| 1, 2, 3, 4, 6, 7 | 21.91 / 21.78 (68 hits) | 21.80 / 21.75 (80) | 21.81 / 21.72 (153) | 21.85 / 21.66 (149) |

- Agreement within 0.01–0.19 bit up to 22 bits, the Gaussian always slightly *under*estimating the bits (the
  candidates hit slightly less often than it says; the excess grows to ~0.1–0.2 bit at 22 bits).
- **The full conditional D-level: 34.49–34.68 bits** fitted directly to the candidates (rank 11, Mahalanobis
  0.36–0.46), against C − C_S^p = 34.55–34.79 from the unconditional model. The two estimates agree within 0.15 bit.
  The drift above (exact bits exceed the Gaussian by up to 0.19 bit at 22 bits, growing with depth) points to a
  slightly larger value. **The costs below use C(D|S) = 34.8 bits** [E, calibrated to 22 bits; fit 34.5–34.7].
- Expected cores of the four runs: Σ_j 2^(Y_j − 34.8) = 2^−6.8 + 2^−6.7 + 2^−5.7 + 2^−5.7 ≈ 0.056; found 0.
- **The deep tail, to rank 9 of 11** (`theory/p647_deeptail.py`, `logs/p647/deeptail_calib.log`): all 2^30.65
  D-candidates of the four complete calibration runs, exact hits on nested subsets of the cross lags against the
  Gaussian fitted to 300,000 of the same residuals:

  | lags (rank) | exact hits | Gaussian expected | exact / Gaussian bits |
  |---|---|---|---|
  | 1, 2, 3, 4 (4) | 47,038 | 48,986 | 15.13 / 15.07 |
  | + 6 (5) | 3,693 | 4,012 | 18.80 / 18.68 |
  | + 7 (6) | 450 | 484 | 21.84 / 21.73 |
  | + 8 (7) | **44** | **44.0** | 25.19 / 25.19 |
  | + 9 (8) | **5** | **8.0** | 28.33 / 27.65 |
  | + 12 (9) | **1** | **1.15** | 30.65 / 30.45 |
  | + 16 (10), + 24 (11 = all) | 0, 0 | 0.21, 0.06 | — / 32.90, 34.63 |

  **The Gaussian holds to 30.6 bits (rank 9 of 11).** Only the last two dimensions, about 4 bits, are extrapolated.
- **Every S-solution's candidates span the full D-lattice** (`theory/p647_perS_rank.py`): for 599 of 600 sampled
  S-solutions of j = 11514 the D-residual differences of its candidates have rank 11 (one has rank 10); median 840
  candidates per S-solution. So no per-S-solution coset can exclude the target in a direction the pooled Gaussian
  cannot see. This is the condition that fails for the 495 control (§5.2).
- **The central classes are not special** (`theory/p647_deeptail_el.py`, `logs/p647/deeptail_calib_el.log`): split
  by energy class, the two central classes (2^29.79 candidates) give 28 hits against 24.9 expected at rank 7 and
  4 against 4.6 at rank 8; the other classes give 16 against 18.7 and 1 against 3.4. The Gaussian bits differ by
  less than 0.15 between the two groups.
- **Live calibration in production** (each class's candidates are counted on the nested lag subsets after its
  D-check, `p647_prod.deep_counts`; `p647_tally.py` sums them against the calibration Gaussian). First 2^31.34
  candidates: rank 4: 79,256 against 79,268; rank 5: 6,544 against 6,492; rank 6: 800 against 784. Deeper hits
  come in **clusters**: in one class all 8 rank-7 hits are also rank-8 hits, from one S-solution whose candidates
  agree on the first 8 lags. So deep counts have far more than Poisson variance (rank 7: 40 against 71 after two
  classes), and cores, if they come, come in clusters too.

### 2.4 Energy structure and class selection [M] (`theory/p647_restrict.py`, `theory/p647_buckets.py`)

The S-join visits only left-energy classes EL ≡ c (mod 24) (P3), 18 of them per joint solution. Their exact
yield (D-candidates) against their exact work (emitted pairs), j = 4208 / j = 9902:

| EL | share of yield | share of work | yield/work gain |
|---|---|---|---|
| 271 | 0.053 / 0.070 | 0.118 / 0.130 | 0.45 / 0.54 |
| 295 | 0.162 / 0.191 | 0.121 / 0.119 | 1.34 / 1.60 |
| **319** | **0.281 / 0.294** | 0.105 / 0.100 | **2.67 / 2.93** |
| **343** | **0.276 / 0.256** | 0.101 / 0.100 | **2.73 / 2.55** |
| 367 | 0.156 / 0.128 | 0.114 / 0.112 | 1.37 / 1.15 |
| 391 | 0.050 / 0.037 | 0.117 / 0.106 | 0.43 / 0.35 |
| the other 12 | 0.022 / 0.024 | 0.324 / 0.333 | < 0.12 |

- The two classes nearest the pair-weighted energy centre ((ē_a + ē_b) + (645 − ē_c − ē_d))/2 (332 and 326 here)
  hold **55% of the yield for 20% of the work: a 2.7× gain**. The rule picks {319, 343} for both runs, so it can be
  applied to any joint solution in advance (`p647_prod.central_els`).
- Restricting the σ-lists instead (multiplicity floors, energy caps, a score mult ≥ 2^(β(e − ē))): best 1.79×.
  Splitting classes further into (e_a, e_c) energy buckets: no gain beyond the class choice (frontier 2.6× at 50%).
  The yield and the work sit in the same few buckets.

---

## 3. Cost per core [M speeds × E counts, counts calibrated exactly] (answer to (B))

### 3.1 Per joint solution

| mode | emitted pairs | D-candidates | expected cores (C(D\|S) = 34.8) | time (1 core) |
|---|---|---|---|---|
| full join, H 79.4 class (calibration j = 11514, 16846) | 1.3·10¹¹ = 2^36.9 | 2^29.07 | 2^−5.73 | 2,930–2,990 s join, 22.5–23.3 ns/pair (4 workers) |
| full join, H 78.3/78.4 class (j = 4208, 9902) | 1.0·10¹¹ = 2^36.5 | 2^28.0–28.1 | 2^−6.8 | 2,490–2,550 s, 24.4–25.1 ns/pair |
| **two central classes, H 79.4 class (production, 6 joint solutions)** | **2.2·10¹⁰ = 2^34.4** | **2^28.15–28.23** | **2^−6.63** | **460–600 s join + 70–100 s lists** |

### 3.2 Per expected core

- Two central classes, H 79.4 class: about 99 joint solutions per expected core, i.e. **2^41.0 emitted pairs and
  about 15 core-hours per expected core** (single-worker speed; about 20% more under a full 4-worker load).
- Full joins: 2^42.6 pairs and about 40 core-hours per expected core (H 79.4 class); 2^43.3 for the 78.x class.
- Both are inside the (B) budget of 2^45 operations / 100 core-hours. **This is the calibrated cost**: the speed is
  measured, the yield (D-candidates) is exact for every processed joint solution, and the only Gaussian factor
  left is C(D|S) = 34.8 bits, checked against exact counts up to 22 bits (§2.3).
- The D-check costs 2–4 s per joint solution, which is negligible.

### 3.3 The levers, and what they are worth [M]

1. **Energy-class selection** (2.7×, §2.4): the two classes nearest the pair-weighted energy centre.
2. **Symmetry** (6×): the stabiliser of y_B #0 in the slice's group has 18 elements on Z₃₈ (6 after dividing by M),
   acting freely on the joint solutions: 20,928 → **3,488 orbits** (`p647_common.stabiliser`, `canonical_joint`).
   One member per orbit suffices.
3. **Class merging** (2.78×, §4.1): the joint solutions of y_B #0 form complete bicliques under the stage-1A
   signature, so one join of a union of lefts against a union of rights serves every (left, right) pair.
4. List restriction by σ-multiplicity or energy (≤ 1.8×) and energy buckets (no gain over 1): not used.
5. Kernel tuning (branch-free emission, bitset prefilter): no gain (22.5 against 22.4 ns/pair). The kernel is
   bound by the scatter's memory traffic.

---

## 4. Production [M] (`theory/p647_prod.py`, `theory/p647_tally.py`, `logs/p647/prod/`)

### 4.1 Class-merged units [P, M]

- **Bicliques [P].** A joint solution of y_B #0 is a choice (o₀, o₁, o₂, o₃) of one y_A option per block with
  Σ_k sig₁A(o_k) = T₁A (the stage-1A S- and D-signatures). With left ℓ = (o₀, o₁) and right r = (o₂, o₃) this reads
  sig(ℓ) + sig(r) = T₁A. So all lefts with one signature s and all rights with T₁A − s form a **complete biclique**:
  every (ℓ, r) in it is a joint solution [M: all 1,230 components of the 20,928 joint solutions are complete; 1,152 of
  them 4 × 4].
- **Merged join.** One 2+2 join of (the union of the left pair sets) against (the union of the right pair sets) finds the
  S-solutions of every (ℓ, r) of a biclique, and nothing else, at the cost of Σ|left| + Σ|right| emitted pairs
  (`p647_join.join4m`: packed entries with 22-bit global indices; self-test 30/30 against a brute force over every
  (ℓ, r); on a real 2 × 4 unit (σ with multiplicity ≥ 4, the two central classes) it returns exactly the union of the 8
  separate joins, 342,800 S-solutions, in 616 s against 1,329 s; `theory/p647_valm.py`,
  `logs/p647/validation_join4m_unit0.json`).
- **Units.** A 4 × 4 biclique holds only 8 orbits (a stabiliser element maps it to itself), so a 2 × 4 sub-biclique
  already covers them. 410 sub-bicliques (384 of shape 2 × 4, 26 of 4 × 4) cover all 3,488 orbits with 2,512 emitted
  sides instead of 6,976: **2.78× less work** (`logs/p647/prod/units.json`, high-freedom units first).
- **D-check** per (ℓ, r) with that joint solution's own four blocks; every D-hit is assembled on Z₆₄₆ and checked by
  the full cyclic PAF; a core is written at once to `logs/p647/CORE_*.json`.

### 4.2 Expected yield [E, calibrated]

With the two central classes, a joint solution yields 2^28.2 D-candidates (H 79.4 class) or about 2^27.2 (78.x
class), i.e. 2^−6.6 or 2^−7.6 expected cores. Over the 3,488 orbits of y_B #0: about 712 × 0.0101 + 2,776 × 0.0050 ≈
**21 expected cores**, for about 2,512 sides × 1.1·10¹⁰ pairs ≈ 2.8·10¹³ emitted pairs (~190 core-hours). That is
about **9 core-hours per expected core** in unit mode. y_B #0 is one of 495 y_B orbits.

### 4.2a Waiting time

- Cores arrive at the calibrated mean rate λ (about 0.46 per hour on three unit workers, 0.6 on four), but full-depth
  hits are likely clustered, like the rank 7–8 hits of §2.3. If clusters have typical size c, the chance of no core
  after expectation λ is about e^(−λ/c) rather than e^(−λ).
- With c = 1–4, a first core is expected within about 2–9 hours of unit-mode running. The H(2092) run found its
  core at 3.06 expected cores (probability 0.05).

### 4.3 Run log

- 06:45–07:37 UTC: per-joint-solution mode, three workers (`w*.jsonl`): 6 joint solutions finished (j = 1999, 1109,
  420, 675, 2483, 78), each 2.2·10¹⁰ pairs, 2^28.15–2^28.23 D-candidates, 0 D-hits; 3 more partial.
- From 07:37 UTC: unit mode (`u*.jsonl`), workers 0–2 (worker 3's core is lent to the control search, §5).

---

## 5. Same-shape positive controls [M] (`theory/p647_control.py`, `theory/p647_ctrl_screen.py`, `theory/p647_ctrl_full.py`)

**Shape.** G = Z₂ × Z_p × Z_q (m = 2pq, v = 2pq + 1, v need not be prime), M = ⟨μ⟩ = M_p × M_q with proper nontrivial
cyclotomic subgroups on both primes and gcd(|M_p|, |M_q|) = 1 (so M is a product and fixes Z₂), exactly as 647 ×353
(|M₁₇| = 4, |M₁₉| = 3). The earlier "controls" 239 ×115 and 443 ×225 were *not* of this shape: there
gcd(|M_p|, |M_q|) = 2 and ⟨μ⟩ is diagonal.

### 5.1 Exact screens

| v | p, q | \|M_p\|, \|M_q\| | μ | result [M] |
|---|---|---|---|---|
| 71 | 5, 7 | 2, 3 | 9 | stage 1 alive both ways; **joint stage empty** (2 y_B orbits along Z₇; 10,316 along Z₅; 0 joint solutions both ways) |
| 111 | 5, 11 | 2, 5 | 9 | stage 1 alive; **joint stage empty** (1 orbit along Z₁₁; 1,503 along Z₅; 0 both ways) |
| 155 | 7, 11 | 2, 5 / 3, 5 | 27 / 9 | **empty** (stage 1 along Z₁₁) |
| 155 | 7, 11 | 3, 2 | 65 | stage 1 alive both ways; stage 1B along Z₁₁ very large (not finished) |
| 183 | 7, 13 | 3, 4 | 109 | **empty** (stage 1 along both primes; 73 canonical α-quadruples, 0 solutions) |
| 191, 231, 239, 267 | various | | 9, 9, 81, 9, 55 | **empty** (one compression has no stage-1 solution) |
| 131, 183, 191, 239, 267 | | | 9, 55, 25, 49, 67, 83, 37 | undecided: one compression alive, the other too large for the α-enumeration |
| **495** | **13, 19** | **4, 3** | **83** | **stage 1 alive; joint stage non-empty** (§5.2) |
| 443 | 13, 17 | 3, 4 | 55 | stage 1 alive along both primes; joint stage not yet run |

### 5.2 The 495 ×83 control (v = 495 = 5·9·11, order 1980)

- The closest shape to 647 ×353: |M₁₃| = 4 (three classes on Z₁₃) and |M₁₉| = 3 (the same six classes on Z₁₉);
  per Z₂-layer a 4 × 7 label matrix (647: 5 × 7); y_B on Z₂₆, y_A on Z₃₈.
- Stage 1B along Z₁₉: 102 y_B orbits. Joint test: orbits #1, #2 and #9 have 2,400, 2,400 and 6,912 joint solutions
  (#0 and #3–#8 none); the stabiliser of y_B #1 acts with 18 elements on Z₃₈, giving 400 orbits.
- Cross level, exact, on every orbit representative (0.3–0.8 s each): the pair-level yield is **0.68 bit below**
  its Gaussian (647: 0.46–0.63) and the distinct-σ S-count 0.45 bit above it (647: 1.5–1.7). C(D|S) ≈ 22.8 bits.
- **Result for y_B #1, #2, #9 (all 1,936 orbits of their joint solutions): 2^25.58 D-candidates, 7.4 expected
  cores by the pooled Gaussian (C(D|S) = 22.7), 0 found** (Poisson probability 0.06%). The pipeline is not the
  cause: 10/10 plants on real 495 joint solutions are recovered (`theory/p647_ctrl_check.py`,
  `logs/p647/control/check_495_83.log`), and an independent candidate-by-candidate count also finds no full-depth hit.
- **Where the Gaussian fails** (400 joint solutions, 2^23.48 candidates, nested cross lags): the exact hits follow
  it to rank 6 (rank 3: 7,271 against 7,846; rank 5: 95 against 126; rank 6: 15 against 19.1). At rank 7 there are
  13 hits in 2^24.3 candidates. At rank 8 (the full condition) there are 0 against 1.6 here and 7.4 over the whole run.
  On the rank-7 hits the residual is a single integer direction t (lag 17 = t, lag 25 = −t), and t took the values
  −2, 4, 6, never 0.
- **Why: discreteness per S-solution** (`theory/p647_perS_rank.py`). A 495 S-solution has a median of 24
  D-candidates. Their residual differences span the full rank-8 D-lattice for only 472 of 1,200 sampled
  S-solutions; the rest have rank 1–7. Exact coset test: for 3,466 of 3,960 S-solutions the D-target lies outside the
  affine sublattice their candidates span. The pooled Gaussian averages over S-solutions and cannot see this.
  A per-S-solution Gaussian removes 0.8 bit; the rest sits in the low-rank cases where any Gaussian is crude.
  **This is the discrete-level failure mode of the design report's §9.4, at the D-level of a small slice.**
- **647 is on the other side of this line**: median 840 candidates per S-solution and full-rank (11) lattices for
  599/600, with the deep tail matching to rank 9 (§2.3). So the 495 control is not a proxy for the 647 D-level, and
  a same-shape slice small enough to run end to end is, for this reason, unlikely to be one. The 495 run continues over
  its remaining y_B orbits (`full_495_83_all.log`), looking for any core in that slice.

---

## 6. Reproduce

```
pip install numba ortools                                        # if missing
python3 theory/p647_join.py --selftest                           # 40/40 x3 + 30/30 against brute force; ejoin4 40/40
python3 theory/p647_cross.py --plant 10 --seed 11                # 10/10 plants on real joint solutions (own energy class)
python3 theory/p647_cross.py --brute 3                           # reduced plants against a 646-lag brute force
python3 theory/p647_cross.py --brute-real 4208 5                 # real S-solutions against brute forces (needs the calibration run)
python3 theory/p647_calib.py 4208 9902 11514 16846               # exact calibration runs (~45 min each, 4 in parallel)
python3 theory/p647_deeptail.py 4208 9902 11514 16846            # deep D-tail to rank 9 (seconds, from the saved S-solutions)
python3 theory/p647_restrict.py 4208; python3 theory/p647_buckets.py 4208 9902   # restriction and bucket analyses
python3 theory/p647_valm.py                                      # join4m equals the separate joins on a real unit
python3 theory/p647_prod.py --units --worker W --workers 4 --els 2   # production, W = 0..3 (restart-safe; touch logs/p647/prod/STOP)
python3 theory/p647_tally.py                                     # status: expected cores, live deep tail, cores found
python3 theory/p647_finish.py logs/p647/CORE_xxx.json            # a core -> data/, verify_h2588.py, pub/H2588/
python3 scripts/verify_h2588.py data/core_v523_b1_q523_p523f_w0.txt /tmp/H2092.npz   # verifier test: SHA-256 9594b975...
python3 theory/p647_control.py --list 700; python3 theory/p647_ctrl_screen.py 700  # product-shape slices, stage-1 screen
python3 theory/p647_control.py 70 9 5 7; python3 theory/p647_control.py 110 9 11 5  # 71 and 111: joint stage empty
python3 theory/p647_ctrl_full.py 494 83 13 19 --yb 1,2,9         # 495 x83 control (all orbits of three y_B)
python3 theory/p647_ctrl_check.py 494 83 13 19 --yb 9 --plants 10 --tail 400    # its plants and D-tail
python3 theory/p647_perS_rank.py                                 # per-S-solution D-lattice ranks, 647 vs 495
```

---

## 7. Hand-off

- **Running** (started by this session, restart-safe, resume by the same command):
  - unit workers 0–2 (`theory/p647_prod.py --units --worker W --workers 4 --els 2`): 410 units (`units.json`,
    high-freedom units first); worker 3's share waits until a fourth worker is started;
  - the 495 ×83 all-orbit control (`theory/p647_ctrl_full.py 494 83 13 19`, `logs/p647/control/full_495_83_all.log`);
  - a paused (SIGSTOP) stage-1 screen of the product slices (`p647_ctrl_screen.py`, PID in the session log).
- **On a core:** `python3 theory/p647_finish.py logs/p647/CORE_*.json` writes `data/core_v647_b1_p647_<tag>.txt`,
  checks it with `nd_doubly.bordered_gs_cyclic`, then with the independent `scripts/verify_h2588.py`
  (`data/H2588_packed.npz`, SHA-256 printed), and builds `pub/H2588/` (sequences, border, normalized matrix, verify.py).
- **If y_B #0 is exhausted without a core** (410 units, ~21 expected cores with two classes, ~190 core-hours):
  - process the remaining energy classes of the same units (`--els 4`, then all: about +60% yield for +120% work);
    the finished-unit record `u{w}.jsonl` must first be made per K, or the finished units are skipped;
  - then resume the joint enumeration over y_B orbits #1–#494 (`nd_joint_exhaust.py`), after adding a cap on the
    S-join (the #1 run died out of memory in the previous session).

