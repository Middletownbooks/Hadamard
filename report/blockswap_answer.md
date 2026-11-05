# H(3452), v = 863: the block-exchanging ansatz (A, B, C, τC) has no stage 1. It is killed at Gate 1.

Date: 2026-11-05. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_blockswap_863.md`.
Scripts: `theory/blockswap_*.py`. Logs: `logs/blockswap/`.

**Credit.**
- The ansatz (x₀, x₁, x₂, x₃) = (A, B, C, τC) with τA = A, τB = B, its gates and its cost formula are by the
  **OpenAI model** (`report/openai_deficit_answer.md`, Proposal 1; branch `codex/order-3452-compression`, commit
  5b26dae).
- The k4 set-up at this order (G = ℤ₂² × ℤ₅ × ℤ₄₃, the S7 host, the ℤ₄₃-compressed certificate and the border
  `logs/h3452/border.json`) is by **OpenAI Codex**.
- The bordered construction is **Kumar's form of Theorem 7 of [EKLO]**.
- The method of Gate 1 follows the one that decided the S14 sibling (`theory/h863_s14_stage1_exact.py`,
  `report/h863_lift_answer.md` §3.5). All `theory/blockswap_*.py` files are new. The decision scripts share no code
  with `theory/openai_deficit_*.py`, nor with each other.

**Tags.**
- **[P]** proved here;
- **[M]** exact computation (script named);
- **[E]** estimate (the lattice Gaussian `theory/q523_gauss.gauss_bits`, read at the true target).

UNKNOWN from a solver decides nothing. No core and no H(3452) was found.

---

## 0. Bottom line

1. **Gate reached: Gate 1. Verdict: kill [M].** The ansatz has **no H-level (stage-1) solution at v = 863**.
   - All 40 Q tuples that survive the joint Q-row-sum / energy gate were decided **exactly empty** by an exhaustive
     A + B against C join on the 20-lag H-level vector (`theory/blockswap_stage1_exact.py`).
   - A second exhaustive method that shares no code agrees on every tuple, with the same list sizes and minimum
     energies (`theory/blockswap_stage1_check2.py`). It uses the Q-character form of the equations, exact bit-packed
     keys and a different enumeration.
   - The other 8 of the 48 Q tuples are empty because A's Q row sums cannot be reached in the alphabet.
   - Every core of the ansatz compresses to a stage-1 solution. So **the slice (A, B, C, τC) of the S7 host holds no
     core**, hence no H(3452) by this route.
   - **Binding number: 0 stage-1 solutions** (40/40 survivors empty).
   - The synthetic estimates (2^32.9 cores per y, conditional cost 2^45.1) described a y that does not exist. Like
     S14 (2^45.1), this is one more positive estimate at an open order decided exactly as zero.
2. **Gate 0 passed [M].** My own code (`theory/blockswap_gate0.py`) confirms:
   - the identity R_{τC}(g) = R_C(τg), on random sequences over G, by integer PAFs on all 860 lags;
   - that τ is an automorphism of order 2 that commutes with S and maps K to K, and that the target is τ-invariant;
   - **48 Q tuples and 40 joint survivors at 863, and 24 and 0 at 519.** These are element-for-element the same sets
     as the codex script.
3. **Correctness controls [M].**
   - Planted tests: **14/14 recovered by both methods.** 8 are random plants. 6 are drawn from the real lists of three
     survivor classes, including the largest (|C| = 3.7·10⁷), so the decision runs on exactly the real lists.
   - Both methods return the same number of hits on every plant.
   - Every hit, planted or not, is re-checked in pure Python on the four explicit length-20 arrays.
4. **Why it is empty.**
   - **[M] A one-coordinate proof for 3 of the 5 classes (24 of the 40 tuples)** (`theory/blockswap_single_lag.py`).
     No energy-consistent triple of the complete lists reaches the single target value −516 at lag (0, 1) (classes 2
     and 5) or 172 at lag (1, 1) (class 4).
   - **[E] The other two classes** (the two largest) need the full join. There the energy-consistent triples (2^38.7
     and 2^29.8 per tuple) meet five independent H-level equations, which together cost 38.6–38.8 bits. That is
     2^0.06 and 2^−9.0 expected solution *configurations* per tuple.
   - Solutions come in orbits of at least 5 configurations within a tuple [P], and every planted hit set had a
     multiple of 40. So under one solution orbit was expected even in the largest class (§2.4).
   - The Gaussian is not calibrated cell by cell here: three single-lag cells are exactly unreachable where it predicts
     67–101 hits in 20,000 samples. The decision rests only on the exact joins.
5. **Gates 2–4 were not run.** Without a y there is nothing to border or cost, and no solver specification is written
   (the prompt's "stop at the first kill").
6. **Gate 5, transfer screen [M]** (`theory/blockswap_transfer.py`; §3).
   - **Exactly dead at the joint Q/energy gate:** 1583 and 2063 (the OpenAI report's priority transfer hosts), and
     1879. Also dead: 1223 with f = 15 and 1319, whose lag-0 equation has no solution (as known), and the analogue
     519.
   - **Survivors of the necessary gates only:** the analogue 583 (16), 1223 with f = 5 (56), 1423 (80), 1783 (8),
     1039 with f = 9 (48) and f = 3 (80), 1487 (8) and 2215 (16).
   - 863 had 40 such survivors, and all were empty at stage 1.
7. **Next computation (§4).**
   - Close this ansatz at 863. Do not build its solver.
   - If block exchange is to be kept, decide stage 1 exactly first, at the open hosts with the fewest survivors:
     1783 and 1487 (8 each), then 2215. The analogue 583 (16) can serve as a control. The Gate 1 code needs only the
     alphabet, m and the targets generalised.
   - Before any cost estimate, run the exact single-lag test and read the Gaussian expected count of stage-1
     solutions per survivor class (§2.4). They take minutes and seconds.
8. **Compute used:** about one hour of wall-clock time on 4 cores, against a 12-hour box. Gate 1 took about
   14 CPU-minutes for both exhaustive methods. The CP-SAT supplement took about 50 CPU-minutes, nine 300-second
   UNKNOWNs among them.

---

## 1. Gate 0: independent re-derivation [P, M]

**Script:** `theory/blockswap_gate0.py`, which works from the definitions only. **Logs:** `logs/blockswap/gate0_m5.json`,
`gate0_m3.json`, `gate0_compare.json`. Run time 11 s (m = 5) and 3 s (m = 3).

- **The identity [P, M].**
  - Proof: (τC) has PAF Σ_h C(τh) C(τ(h + g)) = Σ_{h'} C(h') C(h' + τg) = R_C(τg), since τ is an automorphism.
  - Check: three random ±1 sequences over G (|G| = 860, all lags), with explicit integer PAFs, at m = 5 and m = 3.
  - Also checked exhaustively: τ is an involutive automorphism of G, it commutes with S = ×4 on ℤ₄₃, it maps K onto K,
    and the target F is τ-invariant.
- **Q-level tuples [M].**
  - The Q-compression z_i(q) = Σ_{k∈K} x_i(q, k) is odd, since |K| = 215.
  - The target is the integer Q-level PAF Σ_i PAF_{z_i}(a) = 872 at a = 0 and 860 at a ≠ 0. In PSD form this is
    (3452, 12, 12, 12): the prompt's "3452·δ" is the trivial-character value of the PSD form.
  - A and B: z(01) = z(10). The fourth block's z is C's with 01 and 10 exchanged.
  - The tuples were enumerated directly in the cell sums, with no Walsh transform.
  - **Normalisation:** every block's total is made positive (no tuple has a block of total 0), and the pair {A, B} is
    unordered. C ↔ τC is *not* identified, as in the codex count; 30 survivors remain modulo it.
  - **Results:** 48 tuples at 863 and 24 at 519.
- **Joint gate [M].** Each block's reachable excess energies Σ w(y² − 1) under its Q row sums come from exhaustive
  enumeration of all rows (numpy), not from a DP. They are joined exactly against the excess 2856 (1512 at 519).
  - Row weights: a τ-fixed row (q = 00, 11) is f₀ + 2p₁ + 2p₂ in both sum and energy. The mirrored A/B row (q = 01)
    and every C row have unweighted sums and energy counted twice.
  - Result: **40/48 survive at 863, 0/24 at 519.**
- **Agreement with the codex scripts [M].** Both the 48-tuple set and the 40-survivor set are equal as sets to those
  of `theory/openai_deficit_symmetry.py` (`logs/blockswap/gate0_compare.json`).
  - One bug was found and fixed **in my code** before the comparison. The first version weighted the row sums of the
    mirrored and C rows by 2, which emptied the gate.
- **Lag-0 subset-sum only:** reachable at both orders. It has 42 cells (four of weight 1, thirty-eight of weight 2) at
  863 and 26 at 519.
- **The 8 non-survivors** are the class A ∈ {(5, 7, 7, 7), (7, 7, 7, 5)}, B = (13, 13, 13, 13), C a permutation of
  (1, 1, 1, 3).
  - A τ-fixed row of A must sum to 7 = f₀ + 2(p₁ + p₂), which no letters reach [P]. Since p₁ + p₂ is even,
    f₀ ≡ 7 (mod 4), so f₀ ∈ {−1, 15, −13, 27, −29, 43, −41}. Then p₁ + p₂ would have to be one of ±4, ±10, ±18 or 24,
    and no two letters sum to any of these. This is the gap quoted in the S14 decision ("a = 7 is not f₀ + 2p₁ + 2p₂",
    `report/h863_lift_answer.md` §3.5).
  - The primary Gate 1 code confirms "EMPTY (Q sums)" on all 8.
- **Symmetry classes [M].** Each block can be translated separately without changing its PAF:
  - C by any element of Q; τC then moves by the τ-image;
  - A and B by the τ-fixed (1,1);
  - and C ↔ τC.
  
  These maps cut the 40 survivors into **5 classes of 8**. The Gate 1 list sizes are constant on each class, a
  consistency check. All 40 were decided anyway.

## 2. Gate 1: exact stage 1 [M]

### 2.1 Primary method (`theory/blockswap_stage1_exact.py`)

For each survivor (zA, zB, zC):

- **A and B.** Every τ-invariant y (11 orbit values in ±{1, 13, 15, 27, 29, 41, 43}) with the block's Q row sums is
  listed, under the bound every solution satisfies: excess ≤ 2856 − min(other blocks). The list is built as a product
  of the q = 00 row (f₀, p₁, p₂), the q = 01 row (5 free cells, mirrored to q = 10 as y(10, u) = y(01, −u)) and the
  q = 11 row.
- **C.** Every y on the 20 cells with C's row sums and 2·excess ≤ 2856 − min(A) − min(B) is enumerated in compiled
  code. C's contribution is **PAF_C + PAF_C ∘ τ**; τC is never an independent list.
- **The join.** C's 20-lag contributions are keyed by a linear 64-bit hash, so equal vectors have equal keys and nothing
  is lost. The distinct keys go into an open-addressing table, and all |A|·|B| pairs are streamed against it.
  - Every key match is re-checked on the integer vectors.
  - Every hit is re-checked again in pure Python on the explicit arrays (A, B, C, τC): all 20 lags, the alphabet,
    τ-invariance and the Q sums.
- **Completeness [P].** Take any core of the ansatz. Its H-level compression satisfies the 20 H-level equations (target
  2936 / −516 / 172, a consequence of the full target), has letters in the alphabet, and has τ-invariant A and B.
  After block negations and A ↔ B its Q tuple is one of the 48 of Gate 0. The 8 non-survivors have no admissible A,
  and on each of the 40 survivors the lists contain every admissible block. So an empty join for all 40 tuples means
  that the ansatz has no core.

**Results** (`logs/blockswap/stage1_exact_w{0..3}.jsonl`, summary `logs/blockswap/stage1_summary.json`).

| class (zA, zB, zC), 8 tuples each | min excess A, B, C | \|A\| | \|B\| | \|C\| | distinct C keys | energy-consistent triples | result | s per tuple |
|---|---|---|---|---|---|---|---|---|
| (1,1,1,1), (13,15,15,15), (1,1,1,3) | 0, 728, 0 | 14,040 | 35,400 | 37,470,000 | 1,478,008 | 2^38.7 | **empty**, 0 key matches | 40–41 |
| (1,3,3,3), (13,13,13,13), (3,5,5,5) | 0, 672, 0 | 3,510 | 37,000 | 13,035 | 697 | 2^25.5 | **empty**, 0 | 1.1–1.4 |
| (3,3,3,3), (13,13,13,15), (3,3,3,5) | 0, 728, 0 | 1,585 | 36,080 | 563,375 | 25,117 | 2^29.8 | **empty**, 0 | 1.0–1.5 |
| (9,9,9,9), (9,11,11,11), (3,3,3,5) | 672, 728, 0 | 255 | 1,795 | 86,625 | 4,009 | 2^24.5 | **empty**, 0 | 0.1–0.4 |
| (9,9,9,11), (11,11,11,11), (1,1,1,3) | 728, 784, 0 | 255 | 1,635 | 1,350,000 | 56,852 | 2^29.1 | **empty**, 0 | 0.8–1.2 |

- Total: 347 CPU-seconds for all 40 tuples, against a target of 4 CPU-hours.
- No key match occurred at all, not even a hash collision.

### 2.2 Second method (`theory/blockswap_stage1_check2.py`), sharing no code

- **Formulation [P].** For ψ ∈ Q^, write x_ψ(u) = Σ_q ψ(q) x(q, u) and P(s) for the ℤ₅ PAF at lags 0, 1, 2.
  - Then Σ_a ψ(a) PAF_H(a, b) = PAF_{ℤ₅}(x_ψ)(b), and (τC)_ψ(u) = C_{ψ∘τ}(−u).
  - So the H-level is exactly nine equations:
    - P(A₀₀) + P(B₀₀) + 2P(C₀₀) = (3452, 0, 0);
    - P(A₁₁) + P(B₁₁) + 2P(C₁₁) = (2764, −688, −688);
    - P(A₀₁) + P(B₀₁) + P(C₀₁) + P(C₁₀) = (2764, −688, −688).
- **Keys.** Exact bit-packed integers, 13 bits per value with values 0–7 in two words and value 8 compared directly. No
  hashing. The C side is lexsorted and searched by binary search.
  - C entries with a value outside ±4096 get a sentinel. This is safe [P]: in any solution every C value lies within
    ±3452, because the lag-0 values are ≥ 0 and sum to the target, and |P(b)| ≤ P(0).
- **Enumeration.** A and B by recursive per-orbit search, with exact pruning of the remaining row sums (range and
  parity) and of the energy. C rows by depth-first search with pruning, then a compiled loop over row quadruples.
- **Re-check.** Each hit is checked on the explicit arrays against the 20-lag definition, by code in this file.
- **Results** (`logs/blockswap/check2_w{0,1}.jsonl`):
  - **all 40 tuples EMPTY**;
  - list sizes and minimum energies equal to the primary method's on every tuple;
  - 498 CPU-seconds in total.

### 2.3 Planted tests and the CP-SAT supplement

- **Random plants** (`theory/blockswap_stage1_exact.py planted 8 863`, `logs/blockswap/stage1_planted.jsonl`).
  - The plant is a random ansatz y: mostly ±1, with 1–3 letters of size 13 or 15 per block. The target is its 20-lag
    vector; the Q tuple and energy are its own.
  - Results: 8/8 recovered by the primary method, and 8/8 by the second (`logs/blockswap/check2_planted.jsonl`), with
    equal hit counts (40–320 per plant, all verified).
- **Full-size plants on real lists** (`planted_real`, `logs/blockswap/stage1_planted_real.jsonl`).
  - Two plants each in survivor classes 1, 3 and 5. C is drawn uniformly from the real capped C list, and (A, B) from
    the real lists to complete the excess 2856 exactly.
  - Each target therefore has lag 0 equal to 2936, and the decision runs on exactly the lists of the real run.
  - Results: **6/6 recovered by both methods**, with equal hit counts (80–640). This includes the largest class
    (|C| = 3.7·10⁷), recovered in 34–38 s by the primary method and 52–55 s by the second.
- **Planted design note.** These plants are recovered by complete enumeration, so there is no privileged bucket to force
  them into. They check correctness, not cost.
- **CP-SAT supplement** (`theory/blockswap_stage1_cpsat.py`, `logs/blockswap/cpsat_real.jsonl`,
  `logs/blockswap/cpsat_planted.jsonl`). This is a third, direct model with 42 alphabet variables, the row sums, 342
  exact AddMultiplicationEquality products and all 20 lags; one worker, seed 863 + i, 300 s.
  - It solves all 8 random planted targets (OPTIMAL in 0.1–71 s), and every returned solution re-checks on all 20 lags.
  - On the 6 full-size plants drawn from real lists it solves only 2 (11.6 s and 124 s). It returns **UNKNOWN at 300 s on
    the other 4**, including both plants of the largest class. So an UNKNOWN here is uninformative even when a solution
    is known to exist.
  - Each of the 10 solutions it found lies in the primary method's complete hit set for the same target
    (`theory/blockswap_cpsat_membership.py`, `logs/blockswap/cpsat_membership.json`).
  - On the real tuples, one representative per class (survivors 0, 8, 16, 24 and 32), it returns **UNKNOWN at 300 s on
    all 5**. **That decides nothing.** The decision rests on the two exhaustive methods.

### 2.4 Why the slice is empty [M for classes 2, 4, 5; E for classes 1, 3]

**Single-lag obstructions [M]** (`theory/blockswap_single_lag.py`, `logs/blockswap/single_lag.json`, 9 CPU-minutes).

- **Method.** For each class representative and each lag g ≠ 0:
  - take every block of the complete energy-capped lists of Gate 1, with C contributing PAF_C(g) + PAF_C(τg);
  - form the exact set of pairs (excess, lag-g value) per block;
  - join the three sets at total excess 2856, and ask whether the single value T(g) is reached.
- This is a relaxation of Gate 1, so an unreachable lag proves the tuple empty without the 20-lag join.

| class | lags g = (a, b) at which T(g) is unreachable |
|---|---|
| 1 (1,1,1,1), (13,15,15,15), (1,1,1,3) | none |
| 2 (1,3,3,3), (13,13,13,13), (3,5,5,5) | (0, 1), (0, 2), (0, 3), (0, 4) |
| 3 (3,3,3,3), (13,13,13,15), (3,3,3,5) | none |
| 4 (9,9,9,9), (9,11,11,11), (3,3,3,5) | every lag with a ≠ 0 except b = 0: (1, 1)–(1, 4), (2, 1)–(2, 4), (3, 1)–(3, 4) |
| 5 (9,9,9,11), (11,11,11,11), (1,1,1,3) | (0, 1), (0, 2), (0, 3), (0, 4) |

- The ℤ₅ multipliers and τ map the lists of a tuple to themselves, so lags in one orbit under g ↦ −g, τg and
  u ↦ ku have equal reachable sets. The four lags of classes 2 and 5 form one orbit. The twelve lags of class 4 form
  two: a ∈ {01, 10} (exchanged by τ) and a = 11, each with b ≠ 0.
- **Not a congruence.** The reachable values cover the target's residues mod 8, 16 and 32. The single value is simply
  missing from the reachable set.
- **Consequence.** These 24 tuples have a one-coordinate proof of emptiness, which agrees with both exhaustive joins.

**The two classes that need the full join [E]** (`theory/blockswap_stage1_estimate.py`, 20,000 uniform
energy-consistent triples per class, read at the true target; `logs/blockswap/stage1_estimate.jsonl`).

| class | log₂ energy-consistent triples | H-level bits (rank 5) | Mahalanobis | log₂ E[solution configurations per tuple] |
|---|---|---|---|---|
| 1 | 38.65 | 38.59 | 0.54 | **+0.06** |
| 3 | 29.79 | 38.77 | 0.52 | −8.98 |
| (2, 4, 5: exactly empty by one lag; the Gaussian gives −12.2, −11.5, −9.7) | | | | |

- **Rank.** The rank is 5, not 9. For each ψ class, P(0) + 2P(1) + 2P(2) = (Σ_u x_ψ(u))² is fixed by the Q tuple, and
  the total energy is fixed.
- **Orbits [P].** Two symmetries preserve the ansatz, the Q sums and the target:
  - the translations of C by ℤ₅ (τC moves with it);
  - the multipliers u ↦ ku of ℤ₅, applied to all blocks (they commute with τ and S).
  
  In every class some C row sum is 1 or 3, so no C is fixed by a ℤ₅ translation. Orbits within a tuple therefore have
  at least 5 configurations, and every planted hit set had a multiple of 40. Even in class 1 the expected number of
  solution orbits was below 1: 2^−2.3 with the proved bound 5, and about 2^−5 at the observed 40.
- **Calibration [M against E].** The exact hit fractions on lag subsets, against the Gaussian of the same coordinates:
  - 30 cells in all: 3 single lags and 3 lag pairs or triples, in each of the 5 classes;
  - the 12 single-lag cells with hits agree within −1.6 … +2.2 bits;
  - the 15 pair and triple cells expect 0.4–3.4 hits and have 0–2;
  - 3 single-lag cells have 0 hits against 67–101 predicted. These are the exact obstructions above.
  
  So the Gaussian is good to about 2 bits where the support has no holes. It cannot see a target value that is missing
  from the support without any congruence reason. For this explanation that is enough; the decision rests only on the
  exact joins.
- **Conclusion.** The exact zero is the expected outcome. The ansatz's appeal came from fibre estimates that assumed a
  y exists. The first gate that conditions on a real y is stage 1 itself, and it was never likely to be populated.

## 3. Gate 5: transfer screen [M]

**Script:** `theory/blockswap_transfer.py`. **Log:** `logs/blockswap/transfer.jsonl`.

- **What the screen applies.** For each host of `report/h863_lift_answer.md` §4 (G = ℤ₂² × ℤ_m × ℤ_p, a multiplier of
  odd order f, e = (p − 1)/f classes) it uses the same involution τ(q₁, q₂, u, t) = (q₂, q₁, −u, t). The letters are
  y = ±1 + f·s with s ∈ {−e, −e + 2, …, e}. Then:
  1. Gate 0: the Q tuples, then the joint Q-row-sum / energy gate, by a capped (sum, excess) DP over the cells of each
     row;
  2. the lag-0 subset-sum with the ansatz weights;
  3. for the record, the unrestricted lag-0 test of `theory/h863_lag0_check.py`.
- **Validation.** The DP reproduces 24/0 at 519 and 48/40 at 863.
- **Normalisation.** As in Gate 0; survivors are also given modulo C ↔ τC.

| host | m, p, f, e | excess | lag 0 (ansatz, unrestricted) | Q tuples | **joint survivors** (mod C ↔ τC) |
|---|---|---|---|---|---|
| 519 (analogue; 2076 known) | 3, 43, 7, 6 | 1512 | yes, yes | 24 | **0** |
| 583 (analogue) | 5, 29, 7, 4 | 1904 | yes, yes | 72 | 16 (12) |
| **863 S7** | 5, 43, 7, 6 | 2856 | yes, yes | 48 | 40 (30), **all empty at stage 1 [M]** |
| 1223 | 5, 61, 15, 4 | 4080 | **no, no** [P, known] | 56 | **0** |
| 1223 | 5, 61, 5, 12 | 4080 | yes, yes | 56 | 56 (42) |
| 1423 | 5, 71, 7, 10 | 4760 | yes, yes | 144 | 80 (60) |
| 1583 | 5, 79, 13, 6 | 5304 | yes, yes | 64 | **0** |
| 1783 | 5, 89, 11, 8 | 5984 | yes, yes | 200 | 8 (6) |
| 2063 | 5, 103, 17, 6 | 6936 | yes, yes | 64 | **0** |
| 1039 | 7, 37, 9, 4 | 3600 | yes, yes | 80 | 48 (36) |
| 1039 | 7, 37, 3, 12 | 3600 | yes, yes | 80 | 80 (60) |
| 1319 | 7, 47, 23, 2 | 4600 | **no, no** [P, known] | 40 | **0** |
| 1487 | 7, 53, 13, 4 | 5200 | yes, yes | 96 | 8 (6) |
| 1879 | 7, 67, 11, 6 | 6600 | yes, yes | 136 | **0** |
| 2215 | 7, 79, 13, 6 | 7800 | yes, yes | 184 | 16 (12) |

- **Reading.**
  - "0" is an exact emptiness of the ansatz at that host [M]. Every one of its stage-1 solutions would have a Q tuple
    passing both gates. This rests on one method (the DP), validated against the exhaustive-row Gate 0 at 519 and 863.
  - A positive count passes necessary conditions only. At 863, 40 such survivors were all empty at stage 1.
  - No deeper gate was run (as instructed).
- **The priority hosts are closed.** 1583 and 2063 are the OpenAI report's priority transfer hosts (six classes, as at
  863), and both are exactly dead at the joint gate. So is 1879. 1223 (f = 15) and 1319 already fail lag 0 [P,
  `report/h863_lift_answer.md` §4].
- **Run time:** 13 minutes on one core for all 15 hosts.

## 4. Next computation

1. **Stop the (A, B, C, τC) ansatz at 863.** Its stage 1 is exactly empty [M]. Gates 2–4 have no input, and no
   production solver should be built.
2. **If block exchange is pursued elsewhere,** decide stage 1 exactly at the surviving hosts, starting with the
   smallest stage 1.
   - Candidates: 1783 and 1487 (8 survivors each), then 2215 (16), with the analogue 583 (16, e = 4, five-cell rows)
     as a control.
   - `theory/blockswap_stage1_exact.py` needs only the alphabet, m and the targets generalised.
   - At 863 the whole decision took 6 CPU-minutes. Hosts with m = 7 have 7-cell C rows, so their C lists must be
     split by energy first.
3. **Order of operations for any new ansatz.**
   - First, Gate 0 (seconds).
   - Then the exact single-lag test (§2.4; minutes) and the Gaussian count of stage-1 solutions per survivor class at
     the true target (seconds). Here they would have predicted the kill: 3 of 5 classes are empty by one lag, and the
     other two expect fewer than one solution orbit.
   - Then the exact stage 1 (minutes).
   - Only then the fibre and deficit estimates of cores per y.
   - Two consecutive 863 slices (S14, and now this one) had attractive synthetic fibre estimates (2^42 and 2^33 cores
     per y) above empty stage-1 levels.

## 5. Reproduce

```bash
python3 theory/blockswap_gate0.py 5 3
python3 theory/blockswap_stage1_exact.py run 0 40 logs/blockswap/stage1_exact.jsonl
python3 theory/blockswap_stage1_check2.py run 0 40 logs/blockswap/check2.jsonl
python3 theory/blockswap_stage1_exact.py planted 8 863 logs/blockswap/stage1_planted.jsonl
python3 theory/blockswap_stage1_exact.py planted_real 2 864 logs/blockswap/stage1_planted_real.jsonl 0
python3 theory/blockswap_stage1_check2.py planted logs/blockswap/stage1_planted.jsonl logs/blockswap/check2_planted.jsonl
python3 theory/blockswap_stage1_summary.py
python3 theory/blockswap_stage1_estimate.py 20000 863
python3 theory/blockswap_single_lag.py
python3 theory/blockswap_stage1_cpsat.py real 0,8,16,24,32 logs/blockswap/cpsat_real.jsonl
python3 theory/blockswap_stage1_cpsat.py planted logs/blockswap/stage1_planted.jsonl logs/blockswap/cpsat_planted.jsonl
python3 theory/blockswap_stage1_cpsat.py planted logs/blockswap/stage1_planted_real.jsonl logs/blockswap/cpsat_planted.jsonl
python3 theory/blockswap_cpsat_membership.py
python3 theory/blockswap_transfer.py
```

The production runs split the survivor indices over four workers (`stage1_exact_w0..3`, `check2_w0..1`), and ran
`planted_real` with seeds 864, 865 and 866 on survivors 0, 16 and 32.
