# Prompt: twisted (block-permuting) symmetries in bordered b = 1 cores for the open orders (self-contained, maximum reasoning)

*For a maximum-reasoning Claude Code session on `middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`.*

- **Tags:** [M] verified by exact computation in the repo (name the script), [P] proved (give the proof),
  [E] Gaussian estimate (name the estimator and its calibration), [L] literature recollection (check it; arXiv and
  OEIS are blocked from these containers).
- **Branch:** other sessions push to it. Fetch and merge before every push. Never kill processes you did not start,
  and never `pkill -f` with a pattern that matches your own shell; kill by PID.

## 0. Goal and deliverables

Construct a Hadamard matrix H(4v) for an open v. The route is a **bordered Goethals–Seidel core (b = 1)** in a
**twisted slice**: a group automorphism u of the core group G that also permutes (and possibly negates) the four
blocks.

**Open v** (H(4v) unknown):
- below 3000: 537, 571, 599, 643, 659, 669, 719, 739;
- 3000–4000: 751, 789, 823, 839, 859, 863, 907, 947, 955, 971, 991.

**Already built here, not on the list:**

| v | H | construction |
|---|---|---|
| 523, 573, 647 | H(2092), H(2292), H(2588) | b = 1 |
| 721, 917 | H(2884), H(3668) | Cooper–Wallis |
| 933 | H(3732) | b = 0 |
| 919 | H(3676) | b = 0 with a 3-cycle twist (§1.1) |
| 883 | H(3532) | b = 0, untwisted d = 21 multiplier slice, flat-pool k-tree, the newest (`report/b0_883_answer.md`) |

Deliverables, in order of value:
- **(A)** A verified H(4v) for an open v.
  - H Hᵀ = 4v·I, checked by an independent script that shares no code with the search (patterns:
    `scripts/verify_h2292.py`, `scripts/verify_h2588.py`, `scripts/verify_h3676.py`).
  - A bundle `pub/H<4v>/`: the core sequences on G, `border.json`, the normalized matrix, and a standalone
    `verify.py` (copy the layout of `pub/H2292/` or `pub/H2588/`).
- **(B)** A twisted b = 1 estimator, validated against exact counts on small same-shape slices (§3, stage A).
- **(C)** A screen of every open v: every twisted slice in range, each with exact row-sum feasibility, a calibrated
  log₂E, and the structure of the search it implies (number of independent lists, list sizes).
- **(D)** For the best slices: a search design, a positive control, cost and production.
- **(E)** Exact obstructions found on the way.

**Where to write.**
- Report: `report/b1_twist_answer.md`.
- Scripts: `scripts/b1tw_*.py`. Logs: `logs/b1tw/`.
- Commit and push as results arrive; containers restart. Commit trailers follow the repository's existing commits;
  no pull requests.

## 1. Why this route

### 1.1 The precedent: H(3676), b = 0 with a 3-cycle twist (`report/b0_919_answer.md`)

- **Ansatz.** U = ⟨8⟩ has order 51 in F₉₁₉^*. A is U-invariant. B is invariant under U³, and C(x) = B(8⁻¹x),
  D(x) = B(8⁻²x). So the generator fixes A and cycles B → C → D.
- **Equations.** One per class of ⟨U, −1⟩ (9 classes): PAF_A(c) + PAF_B(c) + PAF_B(u⁻¹c) + PAF_B(u⁻²c) = 0.
- **Size.** 19 orbits; `scripts/gs_twist.py` gives log₂E = 24.4 and n = 49.6 bits.
- **Search.** A 2-list search: the complete A-table (48,620 entries) against random B, with T_B computed as
  orbit-level quadratic forms βᵀQ_cβ (`scripts/b0_twist3_search.py`). It found two cores in 8 minutes on 3 cores. The
  same search reproduced the exact count 2,688 at the control v = 61, d = 15.
- **Limitation.** At b = 0 the twist needs an arithmetic coincidence: block sums (s_A, s_B, s_B, s_B) with
  s_A² + 3s_B² = 4v, s_A ≡ ±1 (mod d) and s_B ≡ ±1 (mod d/3). Here 4·919 = 1 + 3·35². No other open prime has
  4v − 1 = 3s².

### 1.2 Why b = 1 should suit twists better

1. **Row sums cost nothing.**
   - A b = 1 core has block sums (±2, 0, 0, 0).
   - A twist that permutes B, C, D (3-cycle or transposition) and fixes A is compatible automatically:
     B = C = D = 0.
   - No coincidence like 919's is needed.
2. **Involutions other than −1 exist.**
   - G has order m = v − 1. For the 15 prime open v, v ≡ 3 (mod 4), so m = 2h with h odd and G ≅ Z₂ × (abelian group
     of order h). Aut(G) = Aut(H_h).
   - (Z/h)* has an involution u ≠ −1 exactly when h has two distinct prime factors. That holds for 12 of the 15:
     every one except 719 (h = 359), 839 (h = 419) and 863 (h = 431).
   - **Lemma H** (`report/theory_arm_2000_3000.md`, "Order-2 twists are what would decouple"): with a block swap
     B ↔ C under u, P = (B + C)/2 and R = (B − C)/2 have complementary supports. Given the support, A, D, P, R are four
     independent lists. That is exact 4-list independence for a k-tree, the same trick that made T(103) and T(131)
     findable (`scripts/tmat_kt.py`).
   - At b = 0 and primes v ≡ 3 (mod 4) the only involution is −1. There Lemma H forces the overdetermined propus
     shape or an impossible periodic Golay pair.
   - The composite v have extra room. For 537 (m = 2³·67), 669 (2²·167) and 789 (2²·197), G can be
     Z₂² × Z_h′ or Z₂³ × Z₆₇. Then GL(2,2) ≅ S₃ (or GL(3,2)) supplies automorphisms of order 2 and 3, like the
     GL(2,2) twist of Alpöge's 479 core.
3. **Evidence that solutions sit in twisted sub-slices.**
   - Alpöge's 479 core (b = 3) has a twisted automorphism of order 12 (`report/bge3_construction.md` §3).
   - `report/REPORT.md` lines 63 and 195: "GL(3,2)-twisted automorphisms for the b ≥ 2 groups were not searched".
   - `report/theory_arm_2000_3000.md` §7 item 4: "Add block-twisted b = 1 slices (signed block permutations combined
     with Aut(G)). These are not covered yet."
   - The 599 session (`report/b599_answer.md` §0, item 2): where solutions are scarce, 93% of the exact solutions at
     the same-shape rung v = 215 have every block fixed by a twisted multiplier or an affine map. A Gaussian over
     generic blocks cannot see that structure.
4. **Untwisted b = 1 is costed out.**
   - Cyclic multiplier sweeps: `logs/msweep_record/m*.jsonl` for 2000–3000, `report/next_design_answer.md` §7 for
     750–1000.
   - Noncyclic b = 1 slices: `report/theory_arm_2000_3000.md` §4.5 (573 n ≈ 167; 537 and 739 killed by exact
     quotient screens).
   - Staged Paley and doubly-staged pipelines built 523, 573 and 647, and are exhausted or empty at most remaining
     orders. See `report/slice_screens_20261013.md`, `report/paley_design_answer.md`, `report/paley_571_answer.md`,
     `report/next_design_answer.md`.
   - Twisted slices are new families with their own row-sum, lattice and compression behaviour.

## 2. Definitions (exact)

**b = 1 core.**
- Four ±1 functions A, B, C, D on an abelian group G of order m = v − 1 with
  - Σ_k PAF_k(g) = 4m at g = 0 and −4 at every g ≠ 0;
  - block sums (±2, 0, 0, 0) up to order.
- Equivalently, Σ_k |X̂_k(χ)|² = 4v at every χ ≠ 1, and 4 at χ = 1.
- Library form: `theory/ab_lib.py` (any abelian G, T = 4m at 0 and −4 elsewhere).

**b = 1 array.**
- The core goes into the Goethals–Seidel array with R = reversal (x ↦ −x; the back-identity for cyclic G).
- Add border rows S ⊗ 1ᵀ and columns Tᵀ ⊗ 1, with S, T ∈ H(4).
- The 4 × 4 corner is forced: −S Σᵀ T / 4, with Σ the block sums. Search the S, T pairs until the corner is ±1.
- Examples:
  - `scripts/verify_h2292.py` (cyclic Z₅₇₂ = Z₅₂ × Z₁₁ by CRT);
  - `pub/H2292/border.json`;
  - `theory/nd_doubly.py` (`bordered_gs_cyclic`);
  - for noncyclic G, `theory/control_340.py`, which certifies H(340) from b = 1 over Z₂²×Z₃×Z₇;
  - `report/alpoge_methods_forensics.md` §2 (the authors' `h` builder).

**Twisted slice.**
- Take u ∈ Aut(G) with U = ⟨u⟩ of order d, and g = (π, ε) ∈ B₄, a signed permutation of the blocks with g^d = 1.
  The ansatz is

      X_{π(i)}(u·x) = ε_i X_i(x).

- **Labels.** Each U-orbit carries one of 16 sign labels (the signs of A, B, C, D at the orbit's base point),
  propagated along the orbit by g.
- **Stabilisers.** An orbit of size s < d (a point with a nontrivial stabiliser ⟨u^s⟩) must carry a label fixed by
  g^s. In particular every fixed point of u, including 0, carries a g-fixed label.
- **Equations.** ΣPAF is constant on the classes of ⟨U, −1⟩ acting on G \ 0, so there is one equation, f_c = −4,
  per class.
- **When −1 ∈ U** the classes do not pair as {c, −c}. Imposing the twist then cuts bits but not equations.
- **b = 1 restrictions on g [P]:**
  - The ±2 block must be π-fixed with ε = +1. If ε = −1 its sum would be 0.
  - So g fixes A and acts on {B, C, D}: a 3-cycle (needs 3 | d), a transposition, or sign changes.
  - 4-cycles and double transpositions are impossible: they would force |s_A| = |s_B|.
- **Lemma G** (`report/theory_arm_2000_3000.md` §5, item 2). If a problem has B free bits and n constraint bits,
  and fixing F bits (a "support") leaves independent lists, then each support carries 2^(B − F − n) solutions in
  expectation. A support split pays only if F ≤ B − n.
- **Cost conventions.** n = log₂(#labelled configurations) − log₂E. log₂E is a Gaussian on the class sums, using
  Monte Carlo moments and per-coordinate lattice steps (`scripts/gs_twist.py`, `scripts/b1_multiplier_enum.py`).

## 3. Plan (reorder if you find better)

1. **Stage A: estimator and exact validation** (`scripts/b1tw_estimate.py`).
   - Start from `scripts/gs_twist.py` (b = 0, prime t, cyclic U) and `scripts/b1_multiplier_enum.py` (b = 1,
     untwisted, Z_m). Generalise to:
     - G = Z_m, and the noncyclic abelian groups of order m where available;
     - u of any order d, with orbits that have stabilisers;
     - every g fixing A;
     - the target −4;
     - exact row-sum feasibility (an ILP over orbit contribution vectors, against (±2, 0, 0, 0));
     - the Gaussian.
   - **Validate against exact counts** for every twist type at small v: m = v − 1 ≤ about 60–80, at most about 12
     orbits. Use brute force, or `theory/slicesolve.py`, the exhaustive S-invariant b = 1 slice solver (numba Gray
     code, PSD filter, exact join).
   - **Pay special attention to involution twists.**
     - For T-matrices the twisted Gaussian predicted 2^3–2^4 solutions where the exact count was 0 (a joint parity
       condition; `report/tmat_routes.md` §10.2 and `logs/tmat/twist_d2_calibration_final.txt`).
     - For ±1 blocks the parity differs, so measure it here.
     - Prove any obstruction you find.
   - Report a table of estimate against exact count for each twist type.
2. **Stage B: screen every open v.**
   - Groups: G = Z_{v−1} for all nineteen v, plus the noncyclic abelian groups of order v − 1 (v − 1 with square
     factors: 537, 669, 739, 751, 789, 955, 991).
   - Twists: cyclic U = ⟨u⟩ ≤ Aut(G) up to conjugacy, keeping slices with about 15–150 orbits per block. Every g with
     g^d = 1 that fixes A.
   - Combine u with odd-order untwisted multipliers when they commute.
   - Record per slice: orbits, classes, exact row-sum feasibility, log₂E, n, and the lists the search would use with
     their sizes.
   - **Rank by search cost, not by n alone.** A 3-cycle makes B determine C and D. That leaves a 2-list problem with
     lopsided lists, U-invariant A against U³-invariant B. Example: b = 0 571/d15 has n = 94.6 but needs about 2^54.6
     B per hit (`report/unattempted_2000_3000.md`, update). 919/d51 worked because the slice was tiny.
3. **Stage C: search design for the best slices.**
   - **Transposition (B C) under an involution u ≠ −1:**
     - Fix the support S = {B = C}. Then A, D, P = (B + C)/2 on S, and R = (B − C)/2 on Sᶜ are four independent lists.
       P is u-invariant, R is u-anti-invariant, and R vanishes on Fix(u).
     - Use a k-tree join: `scripts/kt_generic.py` `kt_join4`, `scripts/kt_b1.py` (b = 1 joins with complete pools),
       and the support-split pattern of `scripts/tmat_kt.py`.
     - Check Lemma G before relying on it.
   - **3-cycle on (B C D):** use the 2-list search of `scripts/b0_twist3_search.py` (complete A-table, sampled B,
     orbit-level quadratic forms; about 10⁵–10⁶ B per second per core), adapted to the b = 1 target. It is viable when
     the expected B-samples per hit are at most about 2^32–2^34.
   - **Combinations:** a twist commuting with a Paley or doubly-staged slice (`theory/paley_*`, `theory/nd_doubly.py`,
     `theory/nd_joint_exhaust.py`) adds structure to stage 1 and stage 2. Consider it where the untwisted staged
     slice was alive but too expensive. Examples: 971 ×469, 859 ×289, 947 ×49 (`report/next_design_answer.md` §7).
   - **Small slices (r ≤ about 34 orbits):** decide them exactly with `theory/slicesolve.py`. This is cheap and gives
     certificates.
4. **Stage D: positive control, then production.**
   - Control: a smaller v with the same group shape and twist type, with known solutions or exhaustive counts.
   - Production: 4 workers, restart-safe.
   - On any core: assemble the bordered array, verify independently, commit at once, then bundle.

## 4. Known numbers and pitfalls

- **Twisted Gaussians.**
  - The 3-cycle Gaussian of `gs_twist.py` runs 0.9–1.8 bits low on the family 4v = 1 + 3s² (37/9, 61/15, and 919/51
    measured in production).
  - The involution (d = 2) Gaussian is unreliable (§3, stage A).
  - Calibrate every Gaussian you use against exact counts at a same-shape small case.
- **Flat-pool k-tree.**
  - The method that found new b = 0 cores at 487: `kaggle/kt_join_kaggle.py`, `report/REPORT.md` §4–5. Its per-bucket
    estimate uses 3.42 bits per remaining class, calibrated at 487.
  - Recalibrate it at any new shape: at 883 the measured value is 3.97 bits, and the realistic rate is about 1/64 of
    the logged `exp_hits` (`logs/k883/README.md`). The 883 run confirmed this: its first core arrived at a logged
    cumulative `exp_hits` of 2,011, as 32 tuples (1/63).
  - Solutions arrive in clusters: with reversal-closed pools, each core appears as 2⁴ tuples (x ↦ −x on each block
    keeps every PAF), times 2 when two blocks share a pool. Count expected clusters, not tuples.
- **Do not duplicate active runs.**
  - 883 is built (H(3532), `report/b0_883_answer.md`). Drop it from the screen.
  - 599: b = 5, c = 6 Paley(11) lead and its neighbours (`report/b599_answer.md`).
  - Twisted b = 1 slices of 599 are fair game. Its untwisted b = 1 doubly-staged slice is exactly empty (598).
- **T-matrix and b = 0 routes are mapped.** See `report/tmat_routes.md` §10 and `report/b0_919_answer.md` §4.
- **Compute:** 4 cores, about 15 GB RAM, no GPU.
  - OOM kills are silent: set `ulimit -v` or `RLIMIT_AS` for heavy Python.
  - pip-install `ortools`, `numba`, `sympy` or `scipy` if missing.
  - Kaggle notebooks (`kaggle/README.md`) suit embarrassingly parallel joins.

## 5. Tools (read before writing new code)

| file | what it does |
|---|---|
| `scripts/gs_twist.py` | twisted b = 0 estimator (prime t, cyclic U, all g ∈ B₄), with exact row-sum ILP and brute-force validation |
| `scripts/tmat_twist.py` | twisted T-matrix estimator (the B₄ machinery, conjugacy classes, centralisers) |
| `scripts/b1_multiplier_enum.py`, `scripts/nbits.py` | untwisted b = 1 configurations on Z_m (orbits, row sums (±2,0,0,0), Gaussian n); b = 0 Gaussian with flatness |
| `scripts/b0_twist3_search.py` | 2-list search for a 3-cycle twist with orbit-level quadratic forms; `--exhaustB` for exact counts on controls |
| `scripts/kt_b1.py`, `scripts/kt_generic.py`, `scripts/tmat_kt.py` | b = 1 k-tree over τ-buckets with complete pools; generic 4-list k-tree join; T-matrix support split (P/R lists) |
| `theory/slicesolve.py`, `theory/test_controls.py`, `theory/control_340.py`, `theory/find_controls.py` | exhaustive S-invariant LP/b = 1/b = 0 slice solver and certified controls H(28)…H(340), including noncyclic b = 1 |
| `theory/ab_lib.py`, `theory/ab_stage1.py`, `theory/ab_stage2.py` | staged-system library: any abelian G, b = 0/1/≥ 3 targets, per-axis multiplier slices |
| `theory/nd_doubly.py`, `theory/nd_joint_exhaust.py`, `theory/paley_*.py`, `theory/q523_*.py`, `theory/p647_*.py` | doubly-staged and Paley pipelines and the 523/647 production kernels (read only; copy to adapt) |
| `scripts/verify_h2292.py`, `scripts/verify_h2588.py`, `pub/H2292/`, `pub/H2588/` | b = 1 verification and bundle patterns (border search, forced corner, packed rows, SHA-256) |
| `logs/msweep_record/`, `report/next_design_answer.md` §7, `report/slice_screens_20261013.md` | existing untwisted b = 1 slice screens (to avoid re-screening and to pick twist partners) |

## 6. Rules

- A construction counts only if it is verified exactly. Separate theorems, proofs and estimates.
- Calibrate every Gaussian against exact counts. Lattice, parity and 2-adic obstructions are real; the twisted d = 2
  case has already shown one.
- Write restart-safe runners.
  - Arm a check-in so the session never idles while workers run: a heartbeat `Monitor` re-armed on expiry, plus
    `send_later` safety check-ins.
  - An idle session was reclaimed once and killed production workers.
  - Commit results as they arrive. Logs and solutions belong in the repo; large pools and caches belong in the
    scratchpad.
- Before a long run, state the expected cost and the calibration behind it. After the run, compare the observed rate
  with the prediction.
