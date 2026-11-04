# H(3452), v = 863: the S7 host's G/ℤ₅ level is solved in a second, its lift is out of reach (no-go), and the S14 slice is empty

Date: 2026-11-03 to 2026-11-04. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_863_lift.md`.
Scripts: `theory/h863_*.py`. Logs: `logs/h863/`. Solutions: `data/h863_g5_solutions.json`.

**Credit.**
- The construction set-up at this order, the ℤ₄₃-compressed certificate y (`data/h3452_stage1.json`), the border
  (`logs/h3452/border.json`), the projection certificate and the scripts `theory/h3452_*.py` are by **OpenAI Codex,
  directed by the repository owner** (branch `codex/order-3452-compression`, commit 7fb0434). They were brought in
  unchanged in commit 6646b51.
- The bordered construction (`theory/bt_engine.py`, `theory/bt_border.py`) is Kumar's form of [EKLO] Theorem 7.
- The scripts `theory/h863_*.py` written in this session are new. They import neither `bt_engine` nor the Codex
  scripts; the exception is `h863_fibre_estimate.py`, from the previous session.

**Tags.**
- **[P]** proved here (proof in the text);
- **[M]** exact computation (script named);
- **[E]** estimate: the lattice Gaussian `theory/q523_gauss.gauss_bits` read at the true target unless stated
  otherwise, with the sampling distribution named.

UNKNOWN from a solver decides nothing. No core and no H(3452) was found.

---

## 0. Bottom line

1. **Q1: the G/ℤ₅ level inside the Codex fibre is not empty [M].**
   - **Solutions.** 39 distinct solutions w are kept (`data/h863_g5_solutions.json`). They come from a two-stage join
     over the characters of ℤ₂² (`theory/h863_g5_solve.py`, §2); the first solutions appeared after 1.3 s on one core.
   - **Checks.** Every solution was checked on all 172 lags of ℤ₂² × ℤ₄₃ inside the search. It was then checked
     again by `theory/h863_verify_g5.py`, which shares no code with the search
     (`logs/h863/verify_g5_solutions.json`, 39/39). That verifier:
     - builds explicit 5 × 6 sign tables;
     - builds explicit S-invariant sign sequences of length 860;
     - checks their ℤ₄₃-compression against y on all 80 cells;
     - checks their ℤ₅-compression on all 172 lags.
   - **The CP-SAT result.** The UNKNOWN after 1200 s (`logs/h863/g5_gate_codex_cert.json`) was a limit of the
     solver, not of the level.
   - **Abundance.** The level holds 2^99.0 distinct w per y [E], about 2^98.0 after the measured calibration (§2.3).
     A planted test recovered the planted w in 6 of 6 trials.
2. **Q2: the lift chain [M sizes, E counts].**

   | stage | log₂ size | equations: rank, bits | expected count |
   |---|---|---|---|
   | y-fibre (H = G/ℤ₄₃ fixed) | 339.11 [M] | 59, 275.5 [E] | 2^63.5 cores per y |
   | G/ℤ₅ level (w) | entropy 164.28 [M] (192.99 if every w is counted equally) | 12, 65.2 [E] | 2^99.0 w per y (2^98.0 calibrated) |
   | lift (y, w) → G | **172.7–175.0** [M, 39 real w; median 173.9] | **47, 209.8–210.3** [E, 8 real w]; target in the lattice [M] | **2^−35.0 … 2^−36.9 cores per (y, w)** |

   - **What the lift is.** The remaining fibre is, per block and Z₂²-cell, a 5 × 6 table with both margins fixed. Its
     47 equations are the characters nontrivial on both ℤ₅ and ℤ₄₃ (48, less one Parseval relation).
   - **The deficit.** In the 523 form this stage has D_def = C(mixed | y, w) − log₂ fibre(y, w) = 35.0–36.9, mean
     35.9. That is about 10 bits worse than 523's 26.1. Worse, the D-level here is not a check: it is a fibre of
     2^174 configurations under 47 equations (§3).
3. **Every candidate split was measured (§3.2).** The ℤ₅ Fourier split above is the best:

   | split (S-level → D-level) | h_D (exact) | C(D \| S) | **D_def** | D-fibre |
   |---|---|---|---|---|
   | ℤ₅ Fourier: w → mixed (the chain of the prompt) | 172.7–175.0 | 209.8–210.3 | **35.0–36.9** | 2^174 per (y, w) |
   | ℤ₂ along an involution of Q, given y (the 523 route; Kumar's U/W split) | 69.8 | 121.4–121.5 | **51.5–51.8** | 2^69.8 |
   | ℤ₂ along an involution, inside (y, w) | 14.6–16.8 | 91.7–92.1 | 74.9–77.4 | 2^15–2^17 |
   | negation pairs (symmetric S, antisymmetric D) | 84.0 | 275.6 (no equation involves S alone [P]) | 191.6 | 2^84 |

   - The Gauss-period form of the equations (norm conditions in the sextic field) is an exact restatement, not a
     split (§1.3).
   - A two-level Wagner tree over w-parents (§3.3) still has 210 bits to match in its last join.
4. **Q3: no-go.** The cost formula is that of 523 (§3.3): log₂(cost per core) ≈ D_def + log₂(cost per S-solution +
   D-check). At 523 that is 26.1 + 18.1 = 44.2.
   - **The binding number for the 523 route is D_def = 51.7 > 46.** A core needs 2^51.7 S-solutions. Even if every
     S-solution and every D-check were free, the route would exceed the 2^46 budget. At 523's per-S-solution costs it
     is about 2^70.
   - **For the chain y → w → G the binding number is the D-level fibre.** A core needs 2^35.9 pairs (y, w).
     - Producing them costs 2^44.9 pair visits at the measured 2^8.96 visits (0.22 ms) per w.
     - But each pair then needs a 47-equation, 210-bit search of a 2^174 fibre. The best route found, an inner ℤ₂
       S-level join, costs about 2^80 pair visits per pair. That is about 2^116 per core.
   - No design reaches 2^46. Q1 is feasible for as many (y, w) as wanted (2^98 per y), and the calibration holds
     (item 5). Only the cost criterion fails. Even with free S-solutions and free D-checks it fails by 2^5.7. At
     523's per-S-solution costs it fails by 2^24, and at the costs measured or estimated here by 2^41 or more. No
     production plan is written.
5. **Calibration of the estimator on this host [M vs E].**
   - **G/ℤ₅ level, exact join counts** (`logs/h863/g5_calib.jsonl`):
     - stage α: 411 hits against 613 predicted over 30 rounds, −0.58 bit;
     - stage β: 803 hits against 1042 predicted on 12 α-solutions, −0.38 bit;
     - so the 65-bit level is overestimated by about 1 bit.
   - **Mixed level, exact hit fractions on 4.2M uniform (y, w)-fibre samples** (`logs/h863/tail_calib.jsonl`): within
     0.0–0.3 bit of the Gaussian up to 16.6 bits, and Poisson-consistent at 20 bits.
   - **Chain rule.** C(all | y) − C(w | y) = 210.2 against 209.8–210.3 measured directly in real (y, w) fibres.
   - **Planted recovery** at both stages: 6 of 6 (`logs/h863/g5_planted.jsonl`).
   - **Not done.** A comparison of cores found against expected at the analogues 519 and 583 was out of reach: their
     lifts cost 2^57 or more per core by the same chain (§3.4).
6. **Transfer (§4).** The same chain at the open orders of the prompt's §2 has D_def between 37.9 and 89.9 (ℤ_m
   split) and between 53.9 and 150 (ℤ₂ split). Every one is worse than 863.
   - Two of the hosts tried are empty [P]: 1223 with a multiplier of order 15, and 1319 with order 23. In both, the
     H-level lag-0 equation has no solution in the y-alphabet.
   - The analogue 519 (4n = 2076, a known order) sits in the 523 regime, with D_def ≈ 22 for both splits.
   - **The large-prime compression route for the k4 family is effectively closed at these orders** [E].
7. **The S14 sibling is empty [M] (§3.5).**
   - **Why it mattered.** The S14 host `fbf4f510b2566552` is the S7 host plus τ = (q₁, q₂, u, t) → (q₂, q₁, −u, t),
     with r = 77, so its cores are S7 cores. It is the only k4 structure at 863 in the 523 regime: the U/W scan of
     `report/b3_push_answer.md` §8.5 puts it at 2^−24.6 completions per U-solution, against 2^−61 for S7 r140 and
     2^−38 for S6 r155. Given a synthetic y it has 2^42.4 cores per y [E], and its ℤ₂ split along the τ-fixed involution has h_D = 34.8 and
     C(D|S) = 57.5, so D_def = 22.7 (523: 36.4, 62.5, 26.1).
   - **But stage 1 has no solution.** All 78 Q-level solutions were decided (up to symmetry: 21 row patterns, a
     superset of the census's 9). For 47 the Q-sums are unreachable in the y-alphabet within the energy; the other 31
     are empty by an exhaustive join of the complete energy-capped block lists. A second method, sharing no code,
     agrees (§3.5). Planted recovery: 6/6 and 4/4.
   - **Consequence.** This confirms Codex's unconfirmed CP-SAT INFEASIBLE for the fibres it tested and extends it to
     the whole host. One more positive estimate, 2^45.1 for the host, decided exactly as zero.
8. **Next computation (§5).** Stop the 863 S7 lift; the S14 slice is closed.
   - If the route is to be kept at all, screen the census's alive rich hosts by their **split deficits** rather than
     by E_Q, at about a minute per host with `theory/h863_analog.py` generalised to their groups.
   - For any host that lands in the 523 regime (D_def ≲ 26, D-fibre ≲ 2^40), decide its stage 1 exactly first. That
     takes minutes, as for S14.
   - Only then build a staged solver, and validate it on the 519 analogue.

---

## 1. Set-up and exact algebra [P]

### 1.1 The host and the levels

- **The family and the host.**
  - Family: k4, t = 3. G = ℤ₂² × ℤ₅ × ℤ₄₃, K = ℤ₅ × ℤ₄₃, Q = ℤ₂².
  - Target: Σᵢ PAFᵢ = 3440 at 0, −12 on K \ 0, +4 off K.
  - S: ×4 on ℤ₄₃ (order 7). Class reps of ⟨4⟩ in ℤ₄₃^* are 1, 2, 3, 6, 7, 9, with negation pairs (1, 2), (3, 6),
    (7, 9).
  - Per block there are 20 cells h = (q, u) ∈ H = ℤ₂² × ℤ₅, each holding a singleton e(h) and six class signs
    c_j(h).
- **The Codex y.** It has rows (4, 6, 6, 58).
  - 64 cells have |y| = 1, 13 have |y| = 13 and 3 have |y| = 15.
  - y = e + 7S fixes e (y mod 14) and the row sums S(h) ∈ {−2, 0, 2} of the class signs. So a y-fibre cell is a row of
    6 signs with 3 or 4 of one sign: 20 or 15 rows. The y-fibre has 2^339.11 points [M].
- **Characters.** A character of G is χ = (ψ, ϑ, φ). The equations are indexed by the S-orbits of characters up to
  conjugation, 72 classes in all:
  - φ = 0: the H-level, 12 classes, fixed by y;
  - ϑ = 0, φ ≠ 0: the G/ℤ₅ level, 12;
  - ϑ ≠ 0, φ ≠ 0: the mixed level, 48.
- **Rank count [P, M].** Parseval adds one linear relation once y is fixed. So the remaining rank is 59 = 12 + 47,
  which the lattice reduction confirms: 59, 12 and 47.
- **The G/ℤ₅ level.** Put w_i(q, j) = Σ_u c_ij(q, u). This is the column-sum vector of the 5 × 6 table of block i,
  cell q. Its row sums S_i(q, u) are fixed by y, and every entry is odd. There are 4132–4332 distinct w per (i, q).
  The level's targets are 3392 at 0, −60 at (0, t ≠ 0) and 20 at q ≠ 0.

### 1.2 The ℤ₂² character split of the G/ℤ₅ level [P]

- **The split.** For ψ ∈ Q̂ put z_{i,ψ}(t) = Σ_q ψ(q) z_i(q, t), a ℤ₄₃-sequence constant on the six classes. Then
  Σ_{q_a} ψ(q_a) PAF_z(q_a, t₀) = PAF_{z_ψ}(t₀).
- **The targets.** Transformed the same way, the G/ℤ₅ equations become, for each ψ:
  - Σᵢ PAF(z_{i,ψ}) = 3452·δ for ψ = 0;
  - 3372 at 0 and −80 off 0 for ψ ≠ 0;
  - three equations per ψ, at the lags 1, 3 and 7. The lag-0 equation follows from y, since
    Σ_{t₀} PAF = (row sum)².
- **The two stages.** With the cells of Q paired as {00, 01} | {10, 11}, ψ = (a₁, a₂) gives
  z_ψ = A_{a₂} + (−1)^{a₁} B_{a₂}. Here A₀ = w(00) + w(01), B₀ = w(10) + w(11), A₁ = w(00) − w(01) and
  B₁ = w(10) − w(11).
  - **Stage α** (ψ00, ψ10) depends only on the half sums (A₀, B₀). It is the level ℤ₂ × ℤ₄₃ = G/(ℤ₅ × ⟨(0,1)⟩).
  - **Stage β** (ψ01, ψ11) depends only on the half differences. Given (A₀, B₀), a choice of w(00) fixes w(01), and
    likewise for w(10).
  - The other two pairings give the other two ℤ₂ quotients.
- **Bits [E]** (table-weighted sampling, 40,000 samples; the stage values are logged in `logs/h863/g5_calib.jsonl`):
  - each ψ alone: 19.3 bits, rank 3;
  - stage α: 35.65 bits, rank 6;
  - stage β given α: 29.5–29.6 bits;
  - all 12 equations: 65.2 bits (the pairs share lattice structure, hence 65 < 4 × 19.3).

### 1.3 Gauss periods and norm conditions [P]

- **The field.** For φ ≠ 0, x̂_i(ψ, ϑ, φ) = Σ_h (ψϑ)(h) [e_i(h) + Σ_j c_ij(h) η_j(φ)]. The η_j are the Gauss periods
  of F = ℚ(ζ₄₃)^⟨4⟩, of degree 6.
  - −1 ∉ ⟨4⟩, so F is a CM field and F⁺ is cubic.
  - ord₄₃(2) = 14 and ⟨2⟩ = ⟨4⟩ × ⟨−1⟩. So 2 splits completely in F⁺ and is inert in F/F⁺.
- **Stage 2 in the 523 form.** Stage 2 is Σᵢ Zᵢ Zᵢ^† = 3452 in O_F[H].
  - **At the G/ℤ₅ level, per ψ:** Σᵢ |βᵢ|² = 863 in O_{F⁺}, with βᵢ = α_{i,ψ}/2 (all coefficients are even).
    Modulo each prime above 2 (residue field 𝔽₄, where conjugation is the Frobenius) this gives |β|² ≡ β³. So an odd
    number of the four βᵢ lie outside that prime. That is a parity condition only, and every solution of §2 meets it
    automatically.
  - **At the mixed level, per (ψ, ϑ):** a sum of four Hermitian norms equals 3452 in the real subfield of F(ζ₅), so 12
    rational equations per ψ.
- **No shortcut.** A single norm equation could be solved by factoring ideals, but a sum of four norms has no such
  shortcut. The Gauss-period form is therefore a restatement, not a split.

### 1.4 The negation split never separates [P]

- **The identity.** Put S = x + x(−·) (symmetric) and D = x − x(−·) (antisymmetric). Then
  PAF_S + PAF_D = 2 PAF_x + 2 Σ_h x(−h) x(−h − g) = 4 PAF_x.
- **Consequence.** Every class equation reads (S-part) + (D-part) = 4F. No equation involves S alone, so there is no
  S-level to solve first. All 275.6 bits of stage 2 fall on the D-level.
- **Numbers.** The D-freedom given S is h_D = 84.0 bits exactly (`theory/h863_negsplit.py`). So D_def = 191.6.

## 2. Q1: the G/ℤ₅ level is solvable [M]

### 2.1 Method (`theory/h863_g5_solve.py`)

- **Lists.** Per (i, q), all column-sum vectors w of the 5 × 6 tables with the y row sums, with exact multiplicities
  (a DP over rows; `h863_lib.cell_lists`).
- **Stage α.** Per block, sample N = 1024 half-sum pairs (A₀, B₀), weighted by tables. Join 2 + 2 over the blocks on
  the 6-dimensional integer signature: the ψ00 and ψ10 PAFs at the lags 1, 3 and 7. Every hit is re-checked exactly.
- **Stage β.** For each α-solution, list per block every split of A₀ and of B₀ into feasible w-pairs: 186–3176
  splits per half. Sample 512–2048 split pairs per block, weighted, and join 2 + 2 on the ψ01 and ψ11 signature.
- **Check.** Every hit is assembled and checked by direct PAF on all 172 lags of ℤ₂² × ℤ₄₃ (`h863_lib.g5_check`).

### 2.2 Results

- **Runs.** The first run found 8 solutions in 1.3 s (16 α-hits from one round; `logs/h863/g5_solve_first_run.json`).
  Four further runs over the three pairings of Q (`logs/h863/g5_solve_*.jsonl`) brought the total to 39 distinct
  solutions (`data/h863_g5_solutions.json`).
- **Independent verification** (`theory/h863_verify_g5.py`, pure Python, no imports from the search): 39 of 39
  pass.
  - **Tables.** Gale–Ryser greedy tables realise every w with the y row sums; their margins are re-checked.
  - **Sequences.** The four sequences of length 860 are S-invariant.
  - **Compressions.** Their ℤ₄₃-compressions equal y on all 80 cells.
  - **Lags.** The summed PAF of their ℤ₅-compressions is 3392 / −60 / 20 on all 172 lags.
- **Throughput** (`theory/h863_g5_rate.py`, one core):
  - nb = 512: 19,264 w in 60 s, 3.1 ms and 2^13.3 pair visits per w;
  - nb = 2048: 430,044 w in 94 s, **0.22 ms and 2^8.96 pair visits per w**.
- **Why CP-SAT failed.** The 1200 s UNKNOWN came from a model with 480 class signs and 1200 product variables. The
  level is easy once the equations are split by the characters of Q.

### 2.3 Calibration of the level count [M vs E]

| split of Q | stage α: found / predicted (10 rounds, N = 1024) | stage β: found / predicted per α-solution (nb = 512) |
|---|---|---|
| q₂ | 133 / 204 (−0.62 bit) | 60/86, 67/87, 70/84, 72/85 (−0.24 … −0.51) |
| q₁ | 130 / 204 (−0.65) | 59/85, 78/87, 76/87, 63/88 (−0.16 … −0.53) |
| q₁ + q₂ | 148 / 205 (−0.47) | 78/87, 80/90, 44/88, 56/88 (−0.15 … −1.0) |

- **Reading.**
  - The Gaussian overestimates by 0.58 bit (α) and 0.38 bit (β), about 1 bit over the 65-bit level.
  - So E[#w | y] ≈ 2^{164.28 − 65.24 − 1.0} ≈ 2^98.
- **Typicality.** The 39 solutions have fibres 2^172.7–2^175.0, against a typical 2^174.8. They are typical w.
- **Planted test** (`theory/h863_g5_planted.py`, `logs/h863/g5_planted.jsonl`).
  - **Set-up.** A random configuration is planted on all 80 cells; its y-fibre and its own G/ℤ₅ PAF values are the
    targets. The planted half sums and splits are forced into the lists.
  - **Result.** 6 of 6 recovered at both stages, each equal to the planted w.
  - **Unforced hits.** These scatter around the Gaussian (5 vs 7.3, 39 vs 24.8, 2 vs 0.08), because planted targets
    are atypical (Mahalanobis 1.2–3.5 against 0.1 for the real target). The test checks correctness, not cost.

## 3. Q2: the lift chain [M sizes, E counts]

### 3.1 Sizes and bits per stage

- **Level bits given y** (`logs/h863/level_bits_given_y.json`, 16,000 y-fibre samples, read at the true target, all
  targets in their lattices):

  | level | rank | bits |
  |---|---|---|
  | G/H = ℤ₄₃ (ψ0, ϑ = 0; Codex's prime projection) | 3 | 19.4 |
  | ℤ₂ × ℤ₄₃ (stage α) | 6 | 35.7 |
  | G/ℤ₅ = ℤ₂² × ℤ₄₃ (w) | 12 | 65.3 |
  | G/Q = ℤ₅ × ℤ₄₃ (ψ0, all ϑ) | 15 | 84.5 |
  | G/⟨z⟩ = ℤ₂ × ℤ₅ × ℤ₄₃ (the S-level of the ℤ₂ split) | 30 | 154.0 |
  | G | 59 | 275.5 |

- **Inside a real (y, w).**
  - **Sampler** (`theory/h863_yw_estimate.py`, `logs/h863/yw_estimate.jsonl`): uniform 5 × 6 tables with both margins
    fixed, by an exact DP over rows (`h863_lib.TableSampler`).
  - **Results** for 8 real w:
    - log₂ fibre(y, w) = 173.4–174.8;
    - the remaining class equations have rank 47 and cost 209.8–210.3 bits (Mahalanobis 0.25–0.70);
    - the target lies in their lattice, so there is no congruence obstruction;
    - log₂ E[cores | y, w] = −35.0, −36.9, −36.2, −36.4, −35.5, −35.2, −35.9 and −36.0.
  - **All 39 w** (exact): log₂ fibre(y, w) = 172.7–175.0 (`logs/h863/yw_fibres_all.json`).
- **Consistency.** E[cores | y] = 2^99.0 × 2^−35.5 = 2^63.5, as computed directly
  (`logs/h863/fibre_estimate_codex_cert.jsonl`).

### 3.2 The candidate splits and their deficits

Measured by `theory/h863_split_deficits.py` (given y, and with `--inner` inside (y, w)) and by
`theory/h863_negsplit.py`. Logs: `logs/h863/split_deficits*.jsonl` and `logs/h863/negsplit.json`.

- **ℤ₅ Fourier split (the chain of the prompt).** The S-level is the ϑ = 0 component, w. The D-level is the ϑ ≠ 0
  components of the table columns.
  - h_D = log₂ fibre(y, w) = 173.9 (median), C(D | S) = 210.0, so **D_def = 35.0–36.9**.
  - The D-level is a 47-equation problem on 2^174 configurations, not a check.
- **ℤ₂ split along z ∈ Q given y** (the 523 route; at this host it is Kumar's U/W split).
  - **S = x + x(· + z).**
    - Per cell pair, the two 6-sign rows are uniform with fixed sums.
    - h_D = Σ [log₂(n₀ n₁) − H(r₀ + r₁)] = **69.8 exactly**, about 1.75 bits per pair of rows.
    - C_S = 153.9–154.2 (rank 30), C(D | S) = 121.4–121.5, so **D_def = 51.5–51.8** for each of the three
      involutions.
  - **Why the deficit is twice 523's.** h_D is 21% of the fibre at both orders. But the D-side carries 44% of a stage
    2 that is 1.9 times larger (275.6 bits against 146).
  - **Agreement.** The U/W scan of `report/b3_push_answer.md` §8.5 gave 2^−61 completions per U-solution at this host,
    without y fixed. That is the same order.
- **ℤ₂ split inside (y, w).**
  - **Exact enumeration.** Every table of both cells of each pair was enumerated (up to 3176² ≈ 10⁷ pairs) to get the
    sum-table distribution.
  - **Numbers.**
    - H(S | y, w) = 156.9–159.5 and h_D = 14.6–16.8, because the column sums pin almost every D sign;
    - C_S = 118.0–118.3 (rank 24) and C(D | S) = 91.7–92.1;
    - **D_def = 74.9–77.4** (12 cases: 4 w × 3 involutions);
    - 2^39–2^41.5 inner S-solutions per (y, w).
- **Negation pairs** (§1.4): there is no S-level. h_D = 84.0 and **D_def = 191.6**.
- **Gauss periods** (§1.3): an exact restatement as sums of four Hermitian norms, with parity conditions only.

### 3.3 Solver designs and cost per expected core

**Cost formula** (523, `report/paley_523_answer.md` §3.6): log₂(pair visits per core) ≈ D_def + log₂(S-step + D-check).
- The S-step of a 2 + 2 join with sides 2^n matching F bits is 2^{F + 1 − n} per solution.
- A D-check by a 2 + 2 join over block D-lists costs about 2^{h_D / 2}.
- At 523: 26.1 + 18.1 = 44.2.

1. **Chain y → w → G** (the prompt's design).
   - **Per core.** 2^35.9 pairs (y, w) are needed [E]. Producing them costs 2^35.9 × 2^8.96 = **2^44.9 pair visits**
     (about 160 core-days) [M rate × E count].
   - **Deciding one (y, w).** It holds 2^−35.9 expected cores. The options:
     - a plain 2 + 2 join over the four block fibres (2^43.5 each) costs 2^87 per (y, w);
     - the best route found is the inner ℤ₂ split: an S-level join with sides 2^79.4–2^79.9, yielding 2^39–2^41.5 inner
       S-solutions, then a D-check of a 2^15–2^17 fibre for each (≈ 2^8). That is **2^80.4 per (y, w)**;
     - further inner levels (G/Q first, then G/⟨z⟩) do not change the 2^35.9 pairs (y, w) needed per core. To meet
       the budget, each pair would have to be decided in about 2^10 pair visits.
   - **Total.** About **2^116 pair visits per core**. The (y, w) production alone uses all but 1.1 bits of the
     2^46 budget.
2. **The 523 route: ℤ₂ S/D given y, with a Wagner tree for S.**
   - **S-step.** The second merge must match F = C_S − C_parent: 118.3 bits with ℤ₂ × ℤ₄₃ parents, or 69.5 with G/Q
     parents. 523 had 54.4.
     - With practical sides (n ≈ 38, as 523's K = 32 lists), that is 2^32–2^81 per S-solution.
     - The D-check of a 2^69.8 fibre costs about 2^35.
   - **Total.** 2^51.7 × (2^32 + 2^35) ≈ **2^87**. Even at 523's costs per S-solution it is 2^70.
   - **Floor.** 2^51.7 with every S-solution and every D-check free.
3. **A two-level Wagner tree over w-parents.**
   - **Design.** The first merge is on the 12 w-equations; the second joins the mixed equations (F = 210) over unions
     of parent-pair fibre products (2^87 each).
   - **Sides.** A middle value of the 2-block w-signature (about 59 bits) has about 2^{82.1 − 59} = 2^23 parent pairs.
     So n ≤ 110.
   - **Cost.** 2^{210 + 1 − 110} ≈ **2^101 per core**, with sides that cannot be stored.
4. **The four other candidates of the prompt** reduce to these. The ℤ₂² character split is the ℤ₂ split applied once
   or twice; the last ℤ₂ stage keeps D_def ≥ 51.7. The ℤ₅ Fourier treatment is design 1. Negation pairs and Gauss
   periods give no S-level (§3.2).

### 3.4 Calibration and analogues

- **The estimator on the mixed level** (`theory/h863_tail_sample.py`, `theory/h863_tail_calib.py`,
  `logs/h863/tail_calib.jsonl`). 4,194,304 uniform configurations of the fibre of w₀. The exact fraction hitting the
  target on nested subsets of the lag classes not fixed by (y, w), against the Gaussian fitted to the same samples:

  | lag classes (q, u, t) | hits | exact (bits) | Gaussian (bits) |
  |---|---|---|---|
  | (0,0,1) | 209,110 | 4.33 | 4.33 |
  | + (0,0,3) | 12,109 | 8.44 | 8.43 |
  | + (0,0,7), + (0,1,0) | 1,257 | 11.70 | 11.71 |
  | + (0,1,1) | 49 | 16.39 | 16.59 |
  | (3,2,9) | 137,461 | 4.93 | 4.94 |
  | + (0,1,2) | 4,696 | 9.80 | 9.83 |
  | + (2,2,6) | 154 | 14.73 | 14.76 |
  | + (1,2,2) | 3 | 20.42 | 19.68 |

  - Two more random orders agree to within 0.25 bit up to 14.8 bits, with 3 and 6 hits at about 19.7 bits.
  - The rank drop from adding (0,0,7) is the exact linear relation fixed by (y, w).
  - Together with §2.3 (exact join counts at 30–36 bits, within 1 bit), the deficits above are good to about 1 bit.
  - A "go" through the 523 route would need D_def about 6 bits lower with free S-solutions and D-checks, or 24 bits
    lower at 523's costs. Design 1 fails on its cost per (y, w), not on its deficit.
- **The chain rule** C(mixed | y, w) = C(all | y) − C(w | y) gives 210.2, against 209.8–210.3 measured directly.
  The synthetic-y version of the general module (`theory/h863_analog.py`, read at the mean) reproduces every 863 number
  to within 0.5 bit. This validates the transfer method of §4.
- **The analogues 519 and 583.** Both are in the table of §4.
  - **519** (G = ℤ₂² × ℤ₃ × ℤ₄₃, ×4; 4n = 2076 is known): ℤ₃ split D_def = 21.7, (y, w) fibre 2^70.7; ℤ₂ split
    D_def = 22.3, h_D = 42.0. That is the 523 regime.
  - **583** (ℤ₂² × ℤ₅ × ℤ₂₉, ×7 of order 7, the same subgroup as the prompt's ⟨16⟩): 22.5 and 30.0.
  - **No core-count comparison.** It was not run. Even at 519, the chain costs 2^21.7 × 2^35 ≈ 2^57 per core with
    2 + 2 lift checks, and a 523-style ℤ₂ solver for this group shape does not exist in the repository. The estimator
    is instead calibrated on the 863 level counts and tails above, and by the earlier calibrations of the same
    estimator (`report/paley_523_answer.md` §3.4–3.5; `report/b3_push_answer.md` §6.1).

### 3.5 The S14 sibling: in the 523 regime, but stage 1 is empty [M]

- **The host.** `fbf4f510b2566552` is the S7 host plus τ = (q₁, q₂, u, t) → (q₂, q₁, −u, t) (Codex's `--kind s14`;
  generator diag-swap, ×4, ×4).
  - It has r = 77 orbits per block, e = 41 class equations and log₂ E_Q = 45.1 (census). Every S14 core is an S7 core.
  - The τ-orbits of H are (00,0), (11,0), {(00,±1)}, {(00,±2)}, {(11,±1)}, {(11,±2)} and {(01,u), (10,−u)} for
    u = 0…4: 11 per block.
- **Deficits given a synthetic τ-invariant y** (`theory/h863_s14_estimate.py`, `logs/h863/s14_estimate.jsonl`, read at
  the mean; two replicates):
  - log₂ fibre(y) = 185.6–186.4, C(all | y) = 143.6 (rank 32), so log₂ E[cores | y] = 42.0–42.8;
  - the ℤ₂ split along the τ-fixed involution z = (1,1): C_S = 86.0–86.2 (rank 18), h_D = 34.5–35.1 exactly (4 orbits
    are z-stable), C(D | S) = 57.4–57.6, so **D_def = 22.5–22.9**;
  - this matches the U/W scan of `report/b3_push_answer.md` §8.5, which gave 2^−24.6 completions per U-solution here;
  - it is the 523 regime (h_D 36.4, C(D|S) 62.5, D_def 26.1). The same scan puts the other 863 k4 hosts at 2^−61
    (S7 r140) and 2^−38 (S6 r155). So S14 is the only 863 k4 structure where a 523-style route would be in range,
    and its stage 1 decides everything.
- **Stage 1 decided exactly** (`theory/h863_s14_stage1_exact.py`, `logs/h863/s14_stage1_exact.jsonl`).
  - **Q-level.** The S14 Q-level solutions are enumerated exactly: per block (a, b, c) = (z(00), z(01) = z(10),
    z(11)), with Σ(a² + 2b² + c²) = 872, Σ 2(ab + bc) = 860 and Σ(2ac + 2b²) = 860. Up to block order, per-block sign
    and the translation by (1,1), there are **78**, in 21 row patterns. The census listed 9 of them.
  - **Per Q-level solution.** Each block's τ-invariant y with those Q-sums is listed completely, under the only
    bound every solution satisfies: the block energy is at most 2856 − Σ_{j≠i} min_j. The four lists (up to 37,000
    entries) are joined 2 + 2 on the exact 20-lag H-level vector.
  - **Result: all 78 are empty.**
    - **47** have Q-sums that no y reaches. The alphabet ±1, ±13, ±15, ±27, ±29, ±41, ±43 has gaps (for example
      a = 7 is not f₀ + 2p₁ + 2p₂), or the minimal energies exceed 2856.
    - **31** are empty by the exhaustive join, the largest in 94 s.
  - **Second method** (`theory/h863_s14_stage1_check2.py`, `logs/h863/s14_stage1_check2.jsonl`). It shares no code
    with the first: a recursive per-orbit enumeration with exact reachability pruning, a ψ-form signature (the Z₅ PAFs
    of the ψ-twisted sequences, 9 values) and a CP-SAT table-constraint model. Its result: all 78 decided, 47 empty by the
    Q-sums or the energy and 31 INFEASIBLE, the slowest in 13 s. On every case its list sizes and minimum energies
    equal those of the first method.
  - **Planted recovery:** 6/6 for the first method (`logs/h863/s14_stage1_planted.jsonl`, with all symmetric images
    re-checked) and 4/4 for the second.
  - **Codex's CP-SAT.** Codex reported INFEASIBLE for two fibres (one with full domains, 91 s), unconfirmed. Both
    methods here confirm it and extend it to every Q-level solution. My own direct CP-SAT model of these fibres
    (`theory/h863_s14_stage1.py`) returned UNKNOWN at 20 s, which decides nothing.
- **Consequence.** The S14 host has no core [M]: every core has an S14-invariant stage-1 compression, and there is
  none. The k4 structures left at 863 are the rich S7 and S6 hosts. The S6 r155 host (2^135.2) was not measured
  here. Its U/W split is at 2^−38 with W-lists of 2^36 per block (b3 §8.5), outside the 523 regime.

## 4. Transfer to the other orders [E; P where marked]

**Method.** `theory/h863_transfer_run.py` uses `theory/h863_analog.py`, logs to `logs/h863/analog_transfer.jsonl`.
- Hosts: H = ℤ₂² × ℤ₅ or ℤ₂² × ℤ₇, the large prime ℤ_p, and a multiplier of odd order f on ℤ_p, so that −1 ∉ S and
  the classes pair under negation, as at 863. Where possible e = (p − 1)/f is near 6.
- y is synthetic: the H-level lag-0 energy, random positions. The estimates are read at the mean, with two
  replicates each (they agree within 0.4 bit).
- H(w | y) is exact except where marked (Gaussian entropy, for (m + 1)^e > 2·10⁶).
- `theory/h863_lag0_check.py` (`logs/h863/lag0_check.jsonl`) tests the H-level lag-0 equation exactly.

| order | H | p, f, e | orbits | log₂ E[cores/y] | log₂ E[w/y] | log₂ fibre(y,w) | C(mixed), rank | **D_def ℤ_m chain** | ℤ₂ split: h_D, **D_def** |
|---|---|---|---|---|---|---|---|---|---|
| 519 (analogue, 2076 known) | ℤ₂² × ℤ₃ | 43, 7, 6 | 84 | 55.0 | 76.7 | 70.7 | 92.4, 23 | **21.7** | 42.0, **22.3** |
| 583 (analogue) | ℤ₂² × ℤ₅ | 29, 7, 4 | 100 | 36.4 | 59.0 | 102.3 | 124.8, 31 | **22.5** | 40.7, **30.0** |
| **863 S7** | ℤ₂² × ℤ₅ | 43, 7, 6 | 140 | 63.5 | 98.8 | 174.6 | 209.9, 47 | **35.4** | 69.6, **51.5** |
| 1223 | ℤ₂² × ℤ₅ | 61, 15, 4 | 100 | no y [P] | | | | | |
| 1223 | ℤ₂² × ℤ₅ | 61, 5, 12 | 260 | 162.8 | 250.8* | 379.4* | 467.4, 95 | **88.0** | 168.8, **106.5** |
| 1423 | ℤ₂² × ℤ₅ | 71, 7, 10 | 220 | 115.5 | 205.4* | 303.4* | 393.3, 79 | **89.9** | 135.5, **95.7** |
| 1583 | ℤ₂² × ℤ₅ | 79, 13, 6 | 140 | 45.7 | 95.8 | 177.6 | 227.6, 47 | **50.1** | 71.5, **60.7** |
| 1783 | ℤ₂² × ℤ₅ | 89, 11, 8 | 180 | 71.2 | 134.8 | 254.2 | 317.7, 63 | **63.5** | 103.2, **82.8** |
| 2063 | ℤ₂² × ℤ₅ | 103, 17, 6 | 140 | 37.1 | 94.5 | 178.2 | 235.5, 47 | **57.3** | 71.9, **64.8** |
| 1039 | ℤ₂² × ℤ₇ | 37, 9, 4 | 140 | 28.6 | 66.5 | 172.1 | 210.1, 47 | **37.9** | 57.9, **53.9** |
| 1039 | ℤ₂² × ℤ₇ | 37, 3, 12 | 364 | 229.2 | 287.7* | 639.7* | 698.1, 143 | **58.4** | 230.7, **150.0** |
| 1319 | ℤ₂² × ℤ₇ | 47, 23, 2 | 84 | no y [P] | | | | | |
| 1487 | ℤ₂² × ℤ₇ | 53, 13, 4 | 140 | 22.1 | 65.9 | 174.2 | 218.0, 47 | **43.8** | 59.0, **56.9** |
| 1879 | ℤ₂² × ℤ₇ | 67, 11, 6 | 196 | 46.0 | 110.6 | 293.2 | 357.9, 71 | **64.7** | 99.3, **94.3** |
| 2215 | ℤ₂² × ℤ₇ | 79, 13, 6 | 196 | 39.6 | 109.8 | 294.3 | 364.5, 71 | **70.2** | 100.0, **97.2** |

\* H(w | y) by the Gaussian entropy.

- **Formula.** D_def(ℤ_m chain) = C(mixed | y, w) − log₂ fibre(y, w).
  - The mixed equations number 2(m − 1)e − 1 at about 4.0–5.0 bits each. The bits per equation grow by about
    ½ log₂ f per equation as f grows: at (m, e) = (5, 6), D_def is 35.4, 50.1 and 57.3 for f = 7, 13 and 17.
  - The fibre(y, w) is 16 doubly fixed m × e tables.
  - The deficit grows with every parameter. **863 is the most favourable open case, and every other order is 2.5 to
    55 bits worse.**
- **The empty hosts [P].** A y requires Σ_cells (y² − 1) = 4|G| − 12(p − 1) − 16m, with
  y² − 1 ∈ {f²S² ± 2fS}.
  - 1223, f = 15: the values are 0, 840, 960, 3480 and 3720, and no sum of them reaches 4080.
  - 1319, f = 23 (e = 2): the values are 0, 2024 and 2208, and no sum reaches 4600.
  - Both hosts are empty. This is the obstruction type of 739 in `report/paley_523_answer.md` §7.
- **The cores per y.** E[cores/y] at 1487 (2^22) and 1039 f = 9 (2^29) is in the range where zero counts are common
  (`report/b3_push_answer.md` §6.1).

## 5. What would change the verdict, and the next computation

- **The requirement.** A staged route needs a final split with **D_def ≲ 26 and a D-fibre ≲ 2^40**, together with an
  S-level whose last merge matches F ≲ 55 bits. At 863 the deficits are 35.9 (fibre 2^174), 51.7 (2^70) and 76
  (2^16), and the S-level merges match 69–118 bits.
  - D_def scales with the size of stage 2: C₂ = 149 bits at the 519 analogue against 275.6 at 863. A host qualifies
    only if C₂ is at most about 150 bits **and** E[cores/y] is large.
  - The S14 sibling met this requirement (C₂ = 143.6 bits, 2^42 cores per y, D_def 22.7), and an exact stage-1 gate
    killed it (§3.5). Otherwise the census found no such structure at the open orders: every structure with
    E_Q ≥ 2^40 has r ≥ 77 (`report/b3_push_answer.md` §8).
- **Next computation (cheap).** Rank the census's alive rich hosts (`report/b3_push_answer.md` §8.1) by the minimum
  over their splits of (D_def, log₂ D-fibre), instead of by E_Q.
  - Generalise `theory/h863_analog.py` and `theory/h863_s14_estimate.py` to their groups. They already handle any
    ℤ₂² × ℤ_m × ℤ_p with a multiplier on ℤ_p, and sub-slices such as S14, at about a minute per host.
  - For a host that lands in the 523 regime, decide its stage 1 exactly first, as in §3.5. That gate is cheap
    (minutes), and at 863 it killed the one candidate.
  - Only then build a solver, and validate it first at the 519 analogue (2076, known), where both splits have
    D_def ≈ 22.
- **Not recommended.** More effort on the 863 S7 lift, or other y at 863. A new y moves deficits by at most about
  0.5 bit (523 §5): the Codex y and synthetic y agree to within 0.5 bit here. The S14 slice is empty (§3.5).

## 6. Reproduce

Python 3 with numpy; ortools for `theory/h863_s14_stage1.py` and the earlier `theory/h863_g5_gate.py`. Run from the
repository root.

```bash
python3 theory/h863_g5_solve.py --na 1024 --nb 512 --seed 1 --want 8 --out /tmp/g5.json
```

```bash
python3 theory/h863_verify_g5.py data/h863_g5_solutions.json
```

```bash
python3 theory/h863_g5_calib.py --rounds 10 --out /tmp/g5_calib.jsonl
```

```bash
python3 theory/h863_g5_planted.py --trials 6
```

```bash
python3 theory/h863_g5_rate.py 60 1024 2048
```

```bash
python3 theory/h863_yw_estimate.py --n 8 --N 16000
```

```bash
python3 theory/h863_split_deficits.py --N 16000
```

```bash
python3 theory/h863_split_deficits.py --inner 0,1,2 --N 12000
```

```bash
python3 theory/h863_negsplit.py
```

```bash
python3 theory/h863_tail_sample.py 0 1048576 1 /tmp/V_1.npy
```

```bash
python3 theory/h863_tail_calib.py /tmp/tail.jsonl /tmp/V_1.npy
```

```bash
python3 theory/h863_analog.py 5 43 4 --y data/h3452_stage1.json --N 6000
```

```bash
python3 theory/h863_transfer_run.py 8000
```

```bash
python3 theory/h863_lag0_check.py
```

```bash
python3 theory/h863_s14_estimate.py 12000 2
```

```bash
python3 theory/h863_s14_stage1_exact.py 0 78 /tmp/s14.jsonl
```

```bash
python3 theory/h863_s14_stage1_planted.py 6 3
```

```bash
python3 theory/h863_s14_stage1_check2.py 0 78 600 /tmp/s14c2.jsonl
```

The tail calibration used four 1M-sample runs (seeds 1–4) of `h863_tail_sample.py`. Measured compute for the whole
session is about 1.5 core-hours.
