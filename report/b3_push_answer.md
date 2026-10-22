# Bordered Goethals–Seidel arrays with a character set T at the open orders below 4000: census, two obstructions, a calibration, and no new matrix

Date: 2026-10-21. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_b3_push.md`.
Scripts: `theory/bt_*.py`. Logs: `logs/bt/`.

**Tags.** [P] proved here (proof in the text); [M] exact computation (script named); [E] estimate (estimator,
calibration and target named).

**Credit.** The construction is Kumar's (`search/NOTES.md` §18–19 in https://github.com/vatskishan03/hadamard,
Apache-2.0) and, for T = all nontrivial characters, [EKLO] Theorem 7. The compression chain (Q-level, nested
compression levels, DP over fibres, 2 + 2 joins) follows Kumar's `chain.py`/`cchain4.py` and our `theory/s1z_*`; the
border reduction follows Kumar's `border3.c`/`border_probe.py` and our `theory/ab_border*.py`. All code under
`theory/bt_*.py` is new and credits this in its headers.

---

## 0. Bottom line

<!-- BOTTOM-LINE -->
1. **(A) No Hadamard matrix at an open order.** No structure reached a core.
   - **Rich structures.** Every structure that survives the exact gates with a large calibrated estimate is "rich
     and flat": 77 or more orbits per block, and a next exact gate of 2^48 pairs or more (§8):
     - 839 (k4, r = 156: log2 E_Q = 165.9);
     - 863 (r = 140: 136.9);
     - 955 (r = 147: 103.6, its first finer level passed with 5.5·10⁷ solutions);
     - 599 (k4, r = 300: 76.7);
     - 971 (r = 112: 73.3, first finer level passed);
     - 751 (k4 r = 121: 73.5; k8 r = 105: 66.2);
     - 669 (ℤ₃ family, r = 90: 73.1, first finer level passed);
     - 719 (k6, r = 112: 66.2; k4, r = 360: 64.1);
     - 643 (k4 on ℤ₄ × ℤ₃₂ × ℤ₅, S of order 8, r = 122: 91.9; census still running);
     - 739 (k4 on ℤ₂ × ℤ₄² × ℤ₂₃, S of order 6, r = 144: 53.6).
   - **Exact gates kill even large estimates.** The 32 richest 971 hosts (r = 106, estimates 2^83–2^98.7) have no
     solution at the level G/N = ℤ₂ × ℤ₄. The 669 ℤ₃ hosts with estimates 2^26–2^85 die at G/N = ℤ₆ or ℤ₉, all
     except S9 r90 (§5).
   - **Structures like the known successes.** 3964, 2636 and 2148 have 46–54 orbits and estimates of 2^30–2^62. The
     census found nothing comparable in that window. The best alive structure with r ≤ 60 is 839 S18 r60 at
     2^11.2. Its cheapest exact level needs 4.9·10¹¹ join pairs for the first of its 20 Q-level solutions alone
     (§8.2).
   - **Local search.** Annealing over orbit signs did not find a core of the 3964 control in 141 restarts (§8.4).
2. **(B) A ranked census** of <!-- NSTRUCT -->81,151<!-- /NSTRUCT --> structures. It covers:
   - k4 at every n ≡ 3 (mod 4). At 739 every 2-part of order 32 is covered except ℤ₂³ × ℤ₄ and ℤ₂⁵ (2^17 and 2^25
     endomorphisms). At 643 the 2-parts of order 128 with at most 16,384 endomorphisms are covered: ℤ₁₂₈, ℤ₂ × ℤ₆₄
     and ℤ₄ × ℤ₃₂ (at most 4,096), and ℤ₂² × ℤ₃₂ and ℤ₈ × ℤ₁₆ (8,192). This part of the census is still running.
   - k6 at 599, 719, 839, 863, 947 and 971.
   - k8 at 599, 719, 751 and 863; at 839 for the 2-parts with at most 16,384 endomorphisms (running; the full enumeration of S over the 2-groups of order 64 did not finish in 8 hours).
   - The ℤ₃ and product families at 669.
   - Split structures of the 1916 shape (§4).
   - Each structure has its exact gate status, a Q-level stratified count read at the true target, and the cost of
     its next exact gate.
   - The border problem is settled for the top k = 4, k4one and ℤ₃ survivors: all 850 Q-level classes tested (86
     structures) have a border (§7).
     - At k = 8, the monomial classes at Q = ℤ₂³ have a border (Alpöge's Singer-type frames).
     - At Q = ℤ₈ (the 751 survivors) the question is open: no Singer-type symmetry exists there [P], and the
       Kronecker frames give no admissible vector.
     - At k = 6 the survivors' classes are open too. The simple frames fail even on Alpöge's 347 compressions,
       whose border exists [M].
3. **(C) Obstructions.**
   - **Two theorems [P]** with proofs (§3), each checked against an independent brute-force count:
     - **Theorem α:** an element of S that fixes K pointwise and has order 3 on Q = ℤ₂² rules out cores. It explains
       several small-order cases where the estimate was positive and the exact count was zero.
     - **Lemma R:** a symmetry transitive on N \ 0 that fixes a complement rules out cores. It covers primitive
       multipliers and Singer cycles.
   - **Coset-Parseval** adapted to T [P].
   - **Families closed** inside the census range [P, M]:
     - k4 at 859 and 907;
     - k4 at 599 and 719 for 8 ≤ r ≤ 160 (the richer structures there have r ≥ 300);
     - k6 at 839 and 947.
   - **Effectively closed** (every alive structure at log2 E_Q < −20, or no structure alive with a finite estimate):
     - k4 at 739 for the 2-parts with at most 2,048 endomorphisms (ℤ₂ × ℤ₄² adds rich hosts at r = 140–144);
     - k4 at 947; k4one and k4split at 751; k6 at 599 and 971; k8 at 599, 719 and 863; the product family at 669.
   - **Exact closure pass [M]** (`logs/bt/closure.sh`, `theory/bt_gates_levels.py`). It runs every exact level of every
     alive structure in these families, each under a time limit, because the cost bounds overestimate badly.
     - It has closed structures at a deeper exact level, mostly at 739. The current counts are in the "deeper
       level" column of §4.1. The 739 survivors are mostly large-orbit structures (r up to 141, estimates near
       2^−300), whose levels pass or exceed the list cap.
     - At 599 k8, all 845 structures that the census had skipped by its rule of thumb were decided exactly, by the
       Q-level gate or by coset-Parseval (`theory/bt_close.py`).
   - **The 669 ℤ₃ family (Kumar's 2148 family) is closed except for two hosts [P, M].**
     - All 239 structures of the census range are dead except 2: 75 by coset-Parseval, 10 by Lemma R, and 152 at an
       exact finer level G/N = ℤ₉ or ℤ₆ (`theory/bt_gates_levels.py`, timed gates, `logs/bt/gates_run_levels.jsonl`).
     - The dead include every host with an estimate from 2^8 to 2^85. The S9 r82 and S12 r64 cases are confirmed by
       an independent CP-SAT checker (§5).
     - The two left are the S9 r90 hosts (2^73.1, on ℤ₂ × ℤ₃² × ℤ₃₇ and ℤ₂ × ℤ₉ × ℤ₃₇). They pass G/N = ℤ₆ with
       1.3·10⁸ weighted solutions, and their ℤ₉ and ℤ₁₈ levels exceed the list caps: rich hosts.
4. **Calibration (§6).** On 561 small-order structures with exact counts from an independent enumerator:
   - the gates never kill a structure that has cores;
   - when cores exist, E_Q is within −3.9 to +0.6 bits of the count;
   - zero counts are common below about 2^30;
   - at the open orders, every exactly decided structure with a positive estimate was a zero, up to 2^85 at 669.
   - Estimates are therefore used only to rank structures that pass every affordable exact gate.
5. **What would reach further (§8.5):** a staged solver that never lists a level, for the rich hosts.
   - Kumar's ℤ₂ U/W split (the 2148 route) is not that solver [E], calibrated on 2148 (2^−14 estimated, 2^−18
     observed). At the hosts with an S-fixed involution it predicts:
     - 2^−13 to 2^−95 completions per U-solution, with W-lists of 2^19 or more per block;
     - or, where completions are commoner, W-lists of 2^27 or more per block.
   - The 839 hosts have no S-fixed involution.
   - A Gauss-period slice along ℤ₄₃ (863) or ℤ₁₉ (839) is untested.
<!-- /BOTTOM-LINE -->

---

## 1. Controls (all reproduced exactly) [M]

`theory/bt_controls.py`, log `logs/bt/controls.log`; `theory/bt_control2148.py`, log `logs/bt/control_2148.log`;
production-path control `logs/bt/search/`; CP-SAT control `logs/bt/xcheck/ctrl2148_z6.log`.

| control | structure | what was reproduced | reference value |
|---|---|---|---|
| 3964 | ℤ₂² × ℤ₂₄₇, S = ⟨(α, ×66), (1, ×35)⟩ (order 27), 46 orbits, 23 equations | Q-level: 1 solution at rows (−50, −26, −2, 28); ℤ₁₉ level (G/ℤ₁₉ = ℤ₂² × ℤ₁₃): **304 solutions**; ℤ₁₃ level inside the planted fibre: **160** combinations; final join of the planted combination: **3 cores** (fibres 348,502 / 78,327 / 227,993 / 42,567; 173 s); border and assembly of every core: **3 of 3 Hadamard** | 304 and 160 (Kumar) |
| 1916 | split: X, W under ⟨(α, 16), (1, 13)⟩ (46 orbits), Y, Z under ⟨(1, 13)⟩ (128 orbits); 79 classes | **16** Q-level solutions at rows (6, 6, 38, 20); the published 479 core is invariant under an order-12 group after a translation | 16 |
| 2636 | ℤ₄² × ℤ₄₁, K = 2ℤ₄² × ℤ₄₁ (non-split), S = ⟨(β, 1), (1, ×10)⟩ (order 15), 54 orbits, 27 equations | **2,954** compressed vectors per sequence at G/ℤ₄₁ = ℤ₄²; border for Kumar's core with our frames, assembled matrix Hadamard; production path: 192 level-A solutions at Kumar's rows (8, 18, 22, 42) | 2,954 (Kumar) |
| 2148 | ℤ₃ × ℤ₁₇₈, K = ℤ₁₇₈, S = ⟨×153⟩ (the eighth powers of GF(89)), T = all characters of ℤ₃, 54 orbits, 27 equations | Kumar's core (`certificate.json`) passes our core check; its Q-compressions get a border from our Paley-H₁₂ frames (frame 0); the assembled matrix of order 2148 is Hadamard (`theory/bt_control2148.py`, `logs/bt/control_2148.log`). The independent CP-SAT level checker finds the level-ℤ₆ system feasible | 2148 exists |

---

## 2. The engine and the tools

- **`theory/bt_engine.py` — structures.** G = ℤ_{m₁} × … × ℤ_{m_k} (any 2-part), K by generators, Q = G/K,
  Γ = 4 Σ_{χ∈T} χ computed from the characters of G trivial on K (or given), S generated by integer matrices (multipliers
  per axis, automorphisms of the 2-part such as α on ℤ₂² and β on ℤ₄², coupled generators), one group for all blocks
  or one per block (split structures). It produces the orbits, the equation classes (meet of the ⟨S_i, −1⟩ orbit
  partitions) with their Γ-targets, the S-invariant subgroups, the compressions G/N with their own classes and targets
  (`Quot`), and the coset-Parseval screen adapted to T (§3.3).
- **`theory/bt_chain.py`, `bt_gates.py` — exact gates.** Q-level lists inside the Parseval box (two-step construction
  for the product family), nested refinement to G/N (product over the coarse cells, per-character energy bounds),
  fibre DP over cells (distinct vectors of a finer level inside a fibre), meet-in-the-middle fibre expansion, 2 + 2
  joins with a linear 64-bit hash (numba, bucketed), final join on all classes.
- **Census drivers.** `bt_census.py` (automorphism groups of 2-groups, subgroups of given index, all S generated by
  ≤ 2 elements up to conjugation), `bt_census_run.py` (k = 4), `bt_census2.py` (k4one, k6, k8, ℤ₃, product),
  `bt_census971.py` (971: tori of GL(2, 11)).
- **Deeper gates.** `bt_ladder.py` (every invariant N ≤ K, coarse to fine), `bt_level1.py` (level-A gate plus
  stratified count), `bt_decide.py` (exhaustive chain), `bt_search.py` (restart-safe production chain with border and
  assembly of every core), `bt_cost.py` (size of the next exact gate).
- **Independent checks (no code shared with the above).** `bt_xcheck_level.py` (CP-SAT on a compressed level),
  `bt_brute_count.py` (exact core count of small structures by full enumeration), `bt_cpsat_level.py`.
- **Screens [P].** `bt_thm_alpha.py` (Theorem α, §3.1), `bt_thm_rigid.py` (Lemma R, §3.2), coset-Parseval in the
  engine (§3.3).
- **Borders.** `bt_border.py` (reduction to the frame problem, Sylvester frames for (q, t) = (4, 3), Paley frames for
  (3, 3), assembly), `bt_border_census.py`.
- **Calibration.** `bt_calib_small.py`, `bt_calib.py`; local search `bt_anneal.py`; feasibility of Kumar's ℤ₂ U/W
  split at the hosts, calibrated on 2148, `bt_uw_scan.py`.

---

## 3. Obstructions [P]

Notation: k = 4 family, Q = G/K of order 4, T = the three nontrivial characters of Q, so Γ(0) = 12 and Γ = −4 off
K; the target is Σᵢ PAFᵢ(g) = −12 on K \ 0 and +4 off K (with 4|G| at 0); n = |G| + 3.

### 3.1 Theorem α

**Theorem α [P].** Let Q = ℤ₂² and |K| > 1. If S contains an element σ that fixes K pointwise and acts with order 3
on Q, the structure has no core.

*Proof.*
1. *Reduction to σ³ = 1.* Write ord σ = 3^a·m′ with 3 ∤ m′ and put τ = σ^{m′}: it fixes K pointwise and still acts
   with order 3 on Q. Then τ³ = 1 + η with η(G) ⊆ K (τ³ is trivial on Q) and η(K) = 0, so η² = 0 and
   (τ³)^j = 1 + jη. As τ³ has 3-power order, 3^b η = 0; as 2G ⊆ K, 2η = 0. So η = 0 and τ³ = 1.
2. *Splitting.* Take g₁ in a nontrivial coset and g₂ = τg₁, g₃ = τ²g₁. Then c = g₁ + g₂ − g₃ ∈ K, τc = c, and
   h_j = g_j − c satisfy 2h_j = 0 and h₁ + h₂ = h₃. So V = {0, h₁, h₂, h₃} ≅ ℤ₂², G = K ⊕ V and τ = id_K ⊕ α.
3. *Shape of a block.* An invariant block is x_K on K and one function w on each of the three other cosets:
   x(k + h_j) = w(k).
4. *The lags in V \ 0.* PAF(h_j) = 2⟨x_K, w⟩ + 2|K|. The target +4 there gives Σᵢ ⟨x_{i,K}, wᵢ⟩ = 2 − 4|K|. Since
   ⟨x_K, w⟩ = −|K| + 2·#{agreements}, the four blocks agree in exactly one point in total. So x_{i,K} = −wᵢ for three
   blocks, and for the fourth x_K = −w + 2εδ_{k₀} with ε = w(k₀).
5. *The characters of K.* For φ ∈ K̂ \ 1 put a = x̂_K(φ) and b = ŵ(φ). The lags in K give
   Σᵢ (|aᵢ|² + 3|bᵢ|²) = 4|G| + 12 = 4n. The lags in h_j + K give Σᵢ (Re āᵢbᵢ + |bᵢ|²) = 0.
6. *The contradiction.* The three blocks with a = −b contribute 0 to the second sum. So the fourth block has
   Re(conj φ(k₀)·b) = 0 for every φ ≠ 1. Then w′ = w(· + k₀) has w′(m) + w′(−m) constant, equal to 2w′(0), so w is
   constant and b = 0 for φ ≠ 1. Then a = 0 on the three blocks and |a|² = 4 on the fourth, so Σᵢ (|aᵢ|² + 3|bᵢ|²) = 4.
   But it should be 4n. ∎

*Checks [M].* `theory/bt_brute_count.py` (independent full enumeration) gives 0 cores for ℤ₂² × ℤ₅ with S = ⟨α⟩
(where the Q-level estimate was 2^12.3) and for the n = 63 structure ⟨(α, 1, ×4)⟩ (σ⁴ = (α, 1, 1); estimates
log2 E_Q = 4.9, log2 E_A = 12.0; the exhaustive chain also gives 0). With K = 0 (n = 7) the statement fails, as it
must: α-invariant cores exist there.

*Use.* The k = 4 census (`theory/bt_census_run.py`, until this session) did not apply it; the post-screen
`theory/bt_thm_alpha.py` removes <!-- ALPHA-COUNTS -->102<!-- /ALPHA-COUNTS --> survivors, all with log2 E_Q ≤ −26 (the Q-level
conditioning already saw most of the rigidity). `theory/bt_census2.py` applied the same screen from the start
(verdict "order-3 element of S trivial on K"); this section supplies its proof.

### 3.2 Lemma R (rigidity)

**Lemma R [P].** Let N ≤ K be S-invariant with |N| ≥ 5. Suppose S (for a split structure, every Sᵢ) contains τ that
acts transitively on N \ 0 and whose fixed points form a complement of N (|Fix τ| = |G|/|N|). Then the structure has
no core. This holds for every family, since the target on K \ 0 is −Γ(0) = −4t < 0.

*Proof.* With G = M ⊕ N and M = Fix τ, every block is constant on m + (N \ 0): call the value c_m, and put
d_m = x(m). For ν₀ ∈ N \ 0,
PAF(ν₀) = Σ_m [(|N| − 2)c_m² + 2c_m d_m] ≥ (|N| − 4)|M| > 0.
But Σᵢ PAFᵢ(ν₀) = −4t. ∎

Transitivity on N \ 0 forces N to be elementary abelian: ℤ_p with a primitive multiplier, or ℤ₁₁² with a Singer
cycle at 971. `theory/bt_thm_rigid.py` screens every survivor (<!-- RIGID-COUNTS -->799 structures killed<!-- /RIGID-COUNTS -->). Brute-force control: ℤ₂² × ℤ₅
with S = ⟨(1, ×2)⟩ has 0 cores.

### 3.3 Coset-Parseval adapted to T

**[P]** (`bt_engine.BStruct.coset_parseval`; the form of `report/b599x4_answer.md` §1–2 for a character set). Let
N be an S-invariant subgroup of prime order p, and Cᵢ the union of the orbits of block i that are not unions of
N-cosets.
- The characters of G nontrivial on N carry the energy |G|·[4n(1 − 1/p) − 4τ_N]. Here τ_N is the number of χ ∈ T
  nontrivial on π(N), so τ_N = 0 when N ≤ K.
- Parseval caps the same energy at |G|·Σᵢ |Cᵢ|(1 − 1/p²) for odd p, and |G|·Σᵢ |Cᵢ| for p = 2.
- So there is no core when Σᵢ |Cᵢ|(1 − 1/p²) < 4n(1 − 1/p) − 4τ_N.

---

## 4. Census [M gates, E counts]

**Method.**
1. **Families.**
   - **k4.** Q of order 4 (ℤ₂² or ℤ₄), T the three nontrivial characters, every n ≡ 3 (mod 4).
   - **k6, k8.** [EKLO] Theorem 7 with k = 6, 8, where k | n − t.
   - **k4one.** T = {1} plus two nontrivial characters.
   - **ℤ₃ (z3all).** Q = ℤ₃, T = all three characters: 669.
   - **prod12.** Q = ℤ₂² × ℤ₃, t = 9: 669.
   - **971.** G = P × ℤ₁₁², with S inside a torus of GL(2, 11): split, Singer, or Z × U (`bt_census971.py`).
   - **k4split.** Split structures of the 1916 shape: S on blocks 1 and 4, S′ ≤ S of index 2 or 3 on blocks 2 and 3
     (`bt_census_split.py`).
     - The bases are the k4 structures with 16 ≤ r ≤ 50 that the screens leave: 2 at 599, 52 at 751, 100 at 839, 118
       at 863 and 232 at 955.
     - The run is slow: the mixed-group estimate takes 20–45 minutes per structure. It covers the bases at 751 and
       part of 839 so far.
     - Every split structure evaluated has log2 E_Q ≤ −16.
2. **Groups.** G = P × O, with P running over the abelian 2-groups of order 2^{v₂(n−t)} whose endomorphisms can be
   listed (≤ 2^20), and K_P over the subgroups of P with the required quotient, up to Aut(P). O is the odd part,
   cyclic factors plus GL(2, p) for a square p².
3. **Symmetry groups.** S runs over the subgroups of Aut(P)_{K_P} × Aut(O) generated by at most two elements, up to
   conjugation by Aut(P)_{K_P}, that preserve Γ ∘ π, with r in the stated range. This includes the couplings of
   automorphisms of Q with multipliers on the odd part (3964, 1916), non-split K (2636) and the ℤ₃ family (2148).
4. **Gates, cheapest first.**
   - The screens [P]: coset-Parseval for T, Theorem α, Lemma R.
   - The Q-level gate [M]: every row pattern, and every Q-compression quadruple meeting Σᵢ PAF = F_K.
   - The Q-level stratified count [E].
   - Deeper exact levels on the survivors, cheapest first by measured cost: `bt_cost.py`, `bt_gates_run.py`,
     `bt_ladder.py`, `bt_level1.py`, and the exhaustive chain `bt_decide.py` for small r.
   - The border problem for every Q-level class of the top survivors (`bt_border_census.py`).

The tables are made by `theory/bt_summary.py` from the logs; the per-structure records are in
`logs/bt/census/*.jsonl` and `logs/bt/census2/*.jsonl`.

### 4.1 Status of every structure, per order and family

<!-- CENSUS-TABLE -->
| n | family | structures | [P] coset-Parseval | [P] Theorem α | [P] Lemma R | prefilter only (4r − 6.5e < −45) | [M] Q-level empty | [M] deeper level / chain empty | undecided (Q-level cap) | alive | alive with log2 E_Q ≥ 0 | best alive: S, r, e, log2 E_Q |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 599 | k4 | 34 | 0 | 10 | 0 | 0 | 18 | 0 | 0 | 6 | 2 | S2, 300, 224, 76.71 |
| 599 | k4one | 12 | 0 | 0 | 0 | 0 | 12 | 0 | 0 | 0 | 0 | — |
| 599 | k6 | 19 | 15 | 0 | 2 | 0 | 0 | 1 | 0 | 1 | 0 | S10, 84, 55, -49.81 |
| 599 | k8 | 2143 | 923 | 0 | 0 | 0 | 1214 | 0 | 0 | 6 | 0 | S9, 104, 53, -18.54 |
| 643 | k4 | 44564 | 40531 | 72 | 106 | 0 | 2376 | 0 | 0 | 1479 | 446 | S8, 160, 86, 116.51 |
| 669 | prod12 | 282 | 0 | 0 | 8 | 0 | 249 | 2 | 10 | 13 | 0 | S20, 57, 32, -14.28 |
| 669 | z3all | 239 | 75 | 0 | 10 | 0 | 0 | 152 | 0 | 2 | 2 | S9, 90, 45, 73.14 |
| 719 | k4 | 17 | 0 | 5 | 0 | 0 | 6 | 0 | 0 | 6 | 1 | S2, 360, 269, 64.09 |
| 719 | k4one | 4 | 0 | 0 | 0 | 0 | 4 | 0 | 0 | 0 | 0 | — |
| 719 | k6 | 59 | 0 | 0 | 4 | 0 | 27 | 1 | 0 | 27 | 5 | S8, 112, 59, 66.18 |
| 719 | k8 | 142 | 38 | 0 | 0 | 0 | 103 | 0 | 0 | 1 | 0 | S14, 90, 89, −∞ |
| 739 | k4 | 6005 | 4080 | 52 | 170 | 0 | 460 | 1160 | 0 | 83 | 2 | S6, 144, 95, 53.64 |
| 751 | k4 | 217 | 0 | 79 | 24 | 0 | 20 | 27 | 0 | 67 | 8 | S8, 121, 63, 73.45 |
| 751 | k4one | 204 | 0 | 0 | 40 | 0 | 40 | 42 | 0 | 82 | 0 | S20, 57, 33, -22.46 |
| 751 | k4split | 54 | 0 | 0 | 8 | 13 | 1 | 0 | 0 | 32 | 0 | S80, 27/36/36/27, 20, -33.19 |
| 751 | k8 | 441 | 134 | 0 | 0 | 0 | 272 | 1 | 0 | 34 | 7 | S10, 105, 54, 66.19 |
| 839 | k4 | 232 | 0 | 52 | 13 | 0 | 89 | 11 | 0 | 67 | 9 | S6, 156, 83, 165.86 |
| 839 | k4one | 204 | 0 | 0 | 16 | 0 | 108 | 16 | 0 | 64 | 4 | S15, 84, 43, 19.59 |
| 839 | k4split | 9 | 0 | 5 | 1 | 2 | 0 | 0 | 0 | 1 | 0 | S90, 16/30/30/16, 15, -16.41 |
| 839 | k6 | 9 | 0 | 0 | 0 | 0 | 9 | 0 | 0 | 0 | 0 | — |
| 839 | k8 | 16644 | 11405 | 0 | 112 | 1857 | 934 | 3 | 0 | 2333 | 83 | S12, 105, 56, 53.02 |
| 859 | k4 | 81 | 39 | 0 | 0 | 0 | 42 | 0 | 0 | 0 | 0 | — |
| 859 | k4one | 106 | 52 | 0 | 0 | 0 | 54 | 0 | 0 | 0 | 0 | — |
| 863 | k4 | 246 | 0 | 44 | 14 | 0 | 134 | 1 | 0 | 53 | 12 | S7, 140, 71, 136.89 |
| 863 | k4one | 224 | 0 | 0 | 20 | 0 | 168 | 0 | 0 | 36 | 8 | S14, 77, 41, 47.85 |
| 863 | k6 | 60 | 0 | 0 | 10 | 0 | 13 | 5 | 0 | 32 | 1 | S15, 90, 45, 3.28 |
| 863 | k8 | 33 | 5 | 0 | 0 | 0 | 27 | 0 | 0 | 1 | 0 | S14, 108, 107, −∞ |
| 907 | k4 | 384 | 175 | 0 | 0 | 0 | 206 | 3 | 0 | 0 | 0 | — |
| 907 | k4one | 580 | 316 | 0 | 0 | 0 | 264 | 0 | 0 | 0 | 0 | — |
| 947 | k4 | 672 | 340 | 35 | 0 | 0 | 294 | 0 | 0 | 3 | 0 | S12, 120, 89, -113.84 |
| 947 | k4one | 1785 | 1105 | 0 | 0 | 0 | 680 | 0 | 0 | 0 | 0 | — |
| 947 | k6 | 18 | 0 | 0 | 0 | 0 | 18 | 0 | 0 | 0 | 0 | — |
| 955 | k4 | 1321 | 485 | 43 | 79 | 0 | 168 | 31 | 0 | 515 | 62 | S8, 147, 76, 103.58 |
| 955 | k4one | 1636 | 648 | 0 | 106 | 0 | 358 | 34 | 0 | 490 | 0 | S24, 63, 34, -3.7 |
| 971 | k4 | 2440 | 978 | 299 | 50 | 0 | 857 | 48 | 0 | 208 | 16 | S10, 112, 57, 73.3 |
| 971 | k6 | 31 | 0 | 0 | 0 | 0 | 25 | 3 | 0 | 3 | 0 | S12, 118, 95, -368.36 |
<!-- /CENSUS-TABLE -->

"alive" means the structure passed every gate run on it. It is not a claim that cores exist: §6 shows that zero
counts are common among alive structures with E_Q below about 2^30.

### 4.2 The best alive structures

<!-- TOP-TABLE -->
| n | family | G = P × O | S order | r | e | log2 E_Q | exact levels passed (N: weighted solutions) | cheapest unrun level (log2 pairs) | border (classes with a border) |
|---|---|---|---|---|---|---|---|---|---|
| 599 | k4 | [4] × [149] | 2 | 300 | 224 | 76.71 | — | — | — |
| 599 | k4 | [2, 2] × [149] | 2 | 299 | 224 | 15.9 | — | — | — |
| 599 | k4 | [2, 2] × [149] | 2 | 300 | 299 | -201.1 | — | — | — |
| 599 | k6 | [2] × [11, 27] | 10 | 84 | 55 | -49.81 | — | N9: 42.3 | open (k = 6, §7) |
| 599 | k8 | [2, 2, 4] × [37] | 9 | 104 | 53 | -18.54 | N37: 5941493760 | N2: 137.1 | — |
| 599 | k8 | [2, 2, 2, 2] × [37] | 9 | 104 | 55 | -32.91 | N37: 16318464 | N2: 137.1 | — |
| 599 | k8 | [2, 8] × [37] | 12 | 104 | 54 | -161.94 | N37: 447086592 | N2: 137.1 | — |
| 643 | k4 | [2, 2, 32] × [5] | 8 | 160 | 86 | 116.51 | — | — | — |
| 643 | k4 | [2, 2, 32] × [5] | 8 | 160 | 86 | 115.18 | — | — | — |
| 643 | k4 | [4, 32] × [5] | 8 | 124 | 65 | 95.4 | N80: 9534720, N40: 141103911168, N32: 66710016 | N16: 61.0 | — |
| 669 | prod12 | [2, 2] × [3, 5, 11] | 20 | 57 | 32 | -14.28 | — | N11: 40.1 | — |
| 669 | prod12 | [2, 2] × [3, 5, 11] | 10 | 93 | 54 | -17.79 | — | N5: 56.9 | — |
| 669 | prod12 | [2, 2] × [3, 5, 11] | 20 | 63 | 38 | -32.53 | — | N11: 52.4 | — |
| 669 | z3all | [2] × [3, 3, 37] | 9 | 90 | 45 | 73.14 | N111: 130553856 | N6: 52.9 | 1/1 |
| 669 | z3all | [2] × [9, 37] | 9 | 90 | 45 | 73.12 | N111: 130553856 | N6: 52.9 | 1/1 |
| 719 | k4 | [4] × [179] | 2 | 360 | 269 | 64.09 | — | — | — |
| 719 | k4 | [2, 2] × [179] | 2 | 359 | 269 | -11.55 | — | — | — |
| 719 | k4 | [2, 2] × [179] | 2 | 360 | 359 | -279.63 | — | — | — |
| 719 | k6 | [2] × [3, 7, 17] | 8 | 112 | 59 | 66.18 | — | N7: 58.3 | open (k = 6, §7) |
| 719 | k6 | [2] × [3, 7, 17] | 8 | 108 | 59 | 35.3 | — | N7: 51.8 | open (k = 6, §7) |
| 719 | k6 | [2] × [3, 7, 17] | 12 | 84 | 47 | 25.6 | — | N17: 50.4 | open (k = 6, §7) |
| 719 | k8 | [2, 2, 2] × [89] | 14 | 90 | 89 | −∞ | — | — | — |
| 739 | k4 | [2, 4, 4] × [23] | 6 | 144 | 95 | 53.64 | N92: 202515200 | N46: 45.5 | — |
| 739 | k4 | [2, 4, 4] × [23] | 6 | 140 | 93 | 37.37 | N92: 3211008, N46: 579863904 | N23: 53.3 | — |
| 739 | k4 | [4, 8] × [23] | 8 | 118 | 82 | -31.41 | N92: 2880, N46: 142652160 | N23: 42.7 | — |
| 751 | k4 | [4] × [11, 17] | 8 | 121 | 63 | 73.45 | — | N11: 48.3 | — |
| 751 | k4 | [2, 2] × [11, 17] | 8 | 121 | 65 | 66.02 | — | N11: 48.3 | — |
| 751 | k4 | [2, 2] × [11, 17] | 8 | 111 | 65 | 49.36 | — | N11: 48.3 | — |
| 751 | k4one | [4] × [11, 17] | 20 | 57 | 33 | -22.46 | — | — | — |
| 751 | k4one | [4] × [11, 17] | 20 | 60 | 34 | -24.49 | — | — | — |
| 751 | k4one | [4] × [11, 17] | 40 | 36 | 20 | -25.14 | — | — | — |
| 751 | k4split | [4] × [11, 17] | 80 | 27/36/36/27 | 20 | -33.19 | — | — | — |
| 751 | k4split | [2, 2] × [11, 17] | 80 | 27/33/33/27 | 21 | -41.77 | — | — | — |
| 751 | k4split | [4] × [11, 17] | 40 | 32/60/60/32 | 34 | -51.56 | — | — | — |
| 751 | k8 | [8] × [3, 31] | 10 | 105 | 54 | 66.19 | — | N31: 68.8 | open (Q = ℤ₈, §7) |
| 751 | k8 | [8] × [3, 31] | 10 | 105 | 53 | 61.06 | — | N31: 68.8 | open (Q = ℤ₈, §7) |
| 751 | k8 | [8] × [3, 31] | 10 | 91 | 48 | 48.64 | — | N31: 51.2 | open (Q = ℤ₈, §7) |
| 839 | k4 | [2, 2] × [11, 19] | 6 | 156 | 83 | 165.86 | — | N19: 78.4 | 20/20 |
| 839 | k4 | [2, 2] × [11, 19] | 6 | 154 | 83 | 157.24 | — | N11: 81.0 | 20/20 |
| 839 | k4 | [2, 2] × [11, 19] | 9 | 110 | 55 | 105.13 | — | N11: 51.5 | 20/20 |
| 839 | k4one | [2, 2] × [11, 19] | 15 | 84 | 43 | 19.59 | — | N19: 51.7 | — |
| 839 | k4one | [2, 2] × [11, 19] | 15 | 84 | 43 | 19.56 | — | N19: 51.8 | — |
| 839 | k4one | [2, 2] × [11, 19] | 15 | 84 | 43 | 19.31 | — | N19: 51.8 | — |
| 839 | k4split | [2, 2] × [11, 19] | 90 | 16/30/30/16 | 15 | -16.41 | — | — | — |
| 839 | k8 | [8, 8] × [13] | 12 | 105 | 56 | 53.02 | — | — | — |
| 839 | k8 | [8, 8] × [13] | 12 | 105 | 56 | 52.98 | — | — | — |
| 839 | k8 | [2, 2, 16] × [13] | 12 | 104 | 58 | 44.48 | — | N52: 45.4 | — |
| 863 | k4 | [2, 2] × [5, 43] | 7 | 140 | 71 | 136.89 | — | N5: 106.7 | 9/9 |
| 863 | k4 | [2, 2] × [5, 43] | 6 | 155 | 85 | 135.2 | — | N43: 98.0 | 21/21 |
| 863 | k4 | [4] × [5, 43] | 7 | 140 | 70 | 133.35 | — | N5: 106.7 | 9/9 |
| 863 | k4one | [2, 2] × [5, 43] | 14 | 77 | 41 | 47.85 | — | N43: 55.6 | 8/8 |
| 863 | k4one | [2, 2] × [5, 43] | 14 | 77 | 41 | 41.45 | — | N43: 55.9 | 8/8 |
| 863 | k4one | [2, 2] × [5, 43] | 14 | 77 | 41 | 41.45 | — | N43: 55.9 | 8/8 |
| 863 | k6 | [2] × [3, 11, 13] | 15 | 90 | 45 | 3.28 | — | N13: 60.9 | open (k = 6, §7) |
| 863 | k6 | [2] × [3, 11, 13] | 12 | 120 | 71 | -8.65 | — | N11: 106.2 | open (k = 6, §7) |
| 863 | k6 | [2] × [3, 11, 13] | 20 | 66 | 37 | -11.31 | — | — | open (k = 6, §7) |
| 863 | k8 | [2, 2, 2] × [107] | 14 | 108 | 107 | −∞ | — | — | — |
| 947 | k4 | [4, 4] × [59] | 12 | 120 | 89 | -113.84 | N59: 955024128 | — | — |
| 947 | k4 | [4, 4] × [59] | 12 | 120 | 119 | -213.35 | N59: 4007424 | — | — |
| 947 | k4 | [2, 2, 2, 2] × [59] | 12 | 120 | 119 | -213.35 | N59: 4007424 | — | — |
| 955 | k4 | [8] × [7, 17] | 8 | 147 | 76 | 103.58 | N119: 252861382656 | N14: 50.2 | 15/15 |
| 955 | k4 | [2, 4] × [7, 17] | 8 | 140 | 75 | 101.0 | N119: 2211840 | N14: 48.1 | 3/3 |
| 955 | k4 | [8] × [7, 17] | 8 | 147 | 77 | 96.99 | N119: 1104150528 | N14: 50.2 | 15/15 |
| 955 | k4one | [8] × [7, 17] | 24 | 63 | 34 | -3.7 | N34: 4054272 | N119: 33.8 | — |
| 955 | k4one | [8] × [7, 17] | 24 | 61 | 35 | -9.58 | N119: 685056000, N34: 768 | N14: 50.0 | — |
| 955 | k4one | [8] × [7, 17] | 16 | 90 | 53 | -10.26 | — | — | — |
| 971 | k4 | [8] × [11, 11] | 10 | 112 | 57 | 73.3 | N22: 786432 | N121: 49.5 | — |
| 971 | k4 | [8] × [11, 11] | 10 | 112 | 57 | 73.3 | N22: 786432 | N121: 49.5 | — |
| 971 | k4 | [8] × [11, 11] | 10 | 112 | 57 | 73.23 | N22: 786432 | N121: 49.5 | — |
| 971 | k6 | [2] × [3, 7, 23] | 12 | 118 | 95 | -368.36 | N23: 4632576 | — | open (k = 6, §7) |
| 971 | k6 | [2] × [3, 7, 23] | 12 | 120 | 119 | -430.58 | — | — | open (k = 6, §7) |
| 971 | k6 | [2] × [3, 7, 23] | 22 | 72 | 39 | −∞ | — | — | open (k = 6, §7) |
<!-- /TOP-TABLE -->

### 4.3 Notes per order

<!-- NOTES-AUTO -->
- **599.** k4: 34 structures, 6 alive, best S2 r = 300 at log2 E_Q = 76.7; k4one: 12 structures, all dead; k6: 19 structures, 1 alive, best S10 r = 84 at log2 E_Q = -49.8; k8: 2143 structures, 6 alive, best S9 r = 104 at log2 E_Q = -18.5.
- **643.** k4: 44564 structures, 1479 alive, best S8 r = 160 at log2 E_Q = 116.5.
- **669.** prod12: 282 structures, 13 alive, best S20 r = 57 at log2 E_Q = -14.3; z3all: 239 structures, 2 alive, best S9 r = 90 at log2 E_Q = 73.1.
- **719.** k4: 17 structures, 6 alive, best S2 r = 360 at log2 E_Q = 64.1; k4one: 4 structures, all dead; k6: 59 structures, 27 alive, best S8 r = 112 at log2 E_Q = 66.2; k8: 142 structures, 1 alive, best S14 r = 90 at log2 E_Q = −∞.
- **739.** k4: 6005 structures, 83 alive, best S6 r = 144 at log2 E_Q = 53.6.
- **751.** k4: 217 structures, 67 alive, best S8 r = 121 at log2 E_Q = 73.5; k4one: 204 structures, 82 alive, best S20 r = 57 at log2 E_Q = -22.5; k4split: 54 structures, 32 alive, best S80 r = 27/36/36/27 at log2 E_Q = -33.2, 13 prefiltered; k8: 441 structures, 34 alive, best S10 r = 105 at log2 E_Q = 66.2.
- **839.** k4: 232 structures, 67 alive, best S6 r = 156 at log2 E_Q = 165.9; k4one: 204 structures, 64 alive, best S15 r = 84 at log2 E_Q = 19.6; k4split: 9 structures, 1 alive, best S90 r = 16/30/30/16 at log2 E_Q = -16.4, 2 prefiltered; k6: 9 structures, all dead; k8: 16644 structures, 2333 alive, best S12 r = 105 at log2 E_Q = 53.0, 1857 prefiltered.
- **859.** k4: 81 structures, all dead; k4one: 106 structures, all dead.
- **863.** k4: 246 structures, 53 alive, best S7 r = 140 at log2 E_Q = 136.9; k4one: 224 structures, 36 alive, best S14 r = 77 at log2 E_Q = 47.9; k6: 60 structures, 32 alive, best S15 r = 90 at log2 E_Q = 3.3; k8: 33 structures, 1 alive, best S14 r = 108 at log2 E_Q = −∞.
- **907.** k4: 384 structures, all dead; k4one: 580 structures, all dead.
- **947.** k4: 672 structures, 3 alive, best S12 r = 120 at log2 E_Q = -113.8; k4one: 1785 structures, all dead; k6: 18 structures, all dead.
- **955.** k4: 1321 structures, 515 alive, best S8 r = 147 at log2 E_Q = 103.6; k4one: 1636 structures, 490 alive, best S24 r = 63 at log2 E_Q = -3.7.
- **971.** k4: 2440 structures, 208 alive, best S10 r = 112 at log2 E_Q = 73.3; k6: 31 structures, 3 alive, best S12 r = 118 at log2 E_Q = -368.4.
<!-- /NOTES-AUTO -->

<!-- NOTES -->
**Structural remarks** (handwritten; the counts are in the tables above):
- **599 and 719 (k4): the window is empty.** n − 3 = 4·149 and 4·179, and K = ℤ_p is cyclic of prime order. So no
  invariant subgroup lies strictly between K and 0, and the Q-level is the only exact level before G itself.
  - The multipliers of ℤ_p have orders dividing 148 = 4·37 and 178 = 2·89, and 3 divides neither.
  - Every S with an element of order 3 therefore fixes K pointwise, and Theorem α kills it.
  - The census went up to r = 400 here. Every S that Theorem α leaves has r ≤ 20 (multipliers of order 37, 74,
    148, 89 or 178) or r ≥ 114.
  - All of these with r ≤ 152 die at the Q-level.
  - What is left is S of order 2 or 4 with r = 225–360: rich, with no exact level to run.
- **859 and 907: closed for r ≤ 160.**
  - At 859 (n − 3 = 8·107), the multipliers of ℤ₁₀₇ have orders 1, 2, 53 and 106. An S without the order-53 part
    has at least 54 orbits on ℤ₁₀₇ and so r ≥ 162. Every structure with r ≤ 160 therefore carries the order-53
    multiplier and has r ≤ 24, and all of them die (coset-Parseval or Q-level). The rich ones (r ≥ 162) were not
    enumerated.
  - At 907 (8·113; the multipliers of ℤ₁₁₃ form ℤ₁₆ × ℤ₇) the window is populated: S of order 14–32 with
    r = 46–72. All 384 structures still die: 175 by coset-Parseval, 206 at the Q-level and 3 at a deeper level.
- **643 and 739: large 2-parts** (2⁷ and 2⁵).
  - Coset-Parseval removes most structures: 2,890 of 3,101 at 643.
  - The 643 survivors with the largest estimates (G = ℤ₄ × ℤ₃₂ × ℤ₅, S of order 8, r = 119, log2 E_Q ≈ 48) have a
    long chain of index-2 invariant subgroups. The first exact levels pass: N of order 80 with 1.9·10⁶ weighted
    solutions, and N of order 32 with 2.1·10⁷ (§5.2 and `logs/bt/gates_run.jsonl`).
  - But with 119 orbits the last chain step leaves fibres of up to 2^44 per block. So these are rich hosts too.
  - The census of the 2-parts with more than 4,096 (643) or 1,024 (739) endomorphisms is still running.
- **669 (ℤ₃ family).**
  - The timed gates at G/N = ℤ₉ and ℤ₆ (ℤ₉ first; it usually decides in under a second, while the pair bound of
    `bt_cost.py` predicted up to 2^37 pairs) closed every alive host except the two S9 r90 hosts. These pass
    G/N = ℤ₆ (1.3·10⁸ weighted solutions), and their ℤ₉ level exceeds the list cap (7.8·10⁷ vectors for one
    block).
  - The product family (Q = ℤ₂² × ℤ₃, t = 9) dies almost entirely at the Q-level: 249 of 282 structures.
- **751.**
  - k4: the cheapest finer level of the best host (S8, r = 121, 2^73.5), N = ℤ₁₁, costs 2^48.3 pairs.
  - k8: the hosts with log2 E_Q ≥ 0 all have P = ℤ₈, where the border problem is open (§7). Their cheapest finer
    level costs 2^51 pairs or more.
- **839 and 863: the only window structures.** These are the orders with live structures at r = 60–78 and
  E_Q = 2^11–2^45. §8.3 takes the best of them as far as the tools go.
- **947.** Closed at k6. At k4, 3 structures are alive, all at log2 E_Q ≤ −113: S of order 12 on ℤ₄² or ℤ₂⁴ with
  r = 120.
- **955.** The most live hosts of any order (546). The largest, S8 r147, passes the level of N of order 119 with
  2.5·10¹¹ weighted solutions.
- **971.** S lies in a torus of GL(2, 11).
  - Lemma R removes the structures in which a Singer cycle acts on ℤ₁₁².
  - The split-torus hosts with r = 106 die at G/N of order 8 (N = ℤ₁₁²).
  - Of the r = 112 hosts, half die at an axis subgroup of order 22.
  - The survivors pass that level with 786,432 weighted solutions.

---

## 5. Exact obstructions closing structures and families [M]

Every item names its script; "dead" means no core [M] for that structure.

*Soundness check of the recorded kills.* An exact gate refines every Q-level solution, and the 2 + 2 join stores at
most 10⁷ solutions per row pattern. So a kill counts only if the Q-level was stored in full.
- All 1,488 structures killed by an exact gate (every `gates_run*`, `ladder_run*` and `level1_*` log, checked on
  2026-10-22) have at most 4,806 Q-level solutions.
- The exhaustive chains saw at most 64 Q-level and 384 level-A solutions.
- The gate code also refuses a truncated join (`bt_gates.Gates.q_level`, `bt_search.Run`). The census only needs
  the count, which stays exact.

### 5.1 Families closed inside the census range

The census range: S generated by at most two elements of Aut(P)_{K_P} × Aut(O)_{K_O} that preserve Γ ∘ π, up to
conjugation by Aut(P)_{K_P}; 2-parts P whose endomorphisms can be enumerated (≤ 2^20); r between the bounds in §4.
The prefilter of `bt_census2.py` (4r − 6.5e < −45, an estimate) is undone by `theory/bt_close.py`, which runs the
exact gates on every skipped structure.

<!-- CLOSED-FAMILIES -->
- **599, k4one: closed [P, M].** All 12 structures are dead: 0 coset-Parseval, 0 Theorem α, 0 Lemma R, 12 Q-level, 0 deeper level.
- **599, k6: effectively closed.** 1 of 19 structures pass the gates run, all with log2 E_Q < −20 (best -49.8).
- **719, k4one: closed [P, M].** All 4 structures are dead: 0 coset-Parseval, 0 Theorem α, 0 Lemma R, 4 Q-level, 0 deeper level.
- **719, k8: effectively closed.** 1 of 142 structures pass the gates run, all with log2 E_Q < −20 (best −∞).
- **751, k4one: effectively closed.** 82 of 204 structures pass the gates run, all with log2 E_Q < −20 (best -22.5).
- **751, k4split: effectively closed.** 32 of 54 structures pass the gates run, all with log2 E_Q < −20 (best -33.2); 13 prefiltered.
- **839, k6: closed [P, M].** All 9 structures are dead: 0 coset-Parseval, 0 Theorem α, 0 Lemma R, 9 Q-level, 0 deeper level.
- **859, k4: closed [P, M].** All 81 structures are dead: 39 coset-Parseval, 0 Theorem α, 0 Lemma R, 42 Q-level, 0 deeper level.
- **859, k4one: closed [P, M].** All 106 structures are dead: 52 coset-Parseval, 0 Theorem α, 0 Lemma R, 54 Q-level, 0 deeper level.
- **863, k8: effectively closed.** 1 of 33 structures pass the gates run, all with log2 E_Q < −20 (best −∞).
- **907, k4: closed [P, M].** All 384 structures are dead: 175 coset-Parseval, 0 Theorem α, 0 Lemma R, 206 Q-level, 3 deeper level.
- **907, k4one: closed [P, M].** All 580 structures are dead: 316 coset-Parseval, 0 Theorem α, 0 Lemma R, 264 Q-level, 0 deeper level.
- **947, k4: effectively closed.** 3 of 672 structures pass the gates run, all with log2 E_Q < −20 (best -113.8).
- **947, k4one: closed [P, M].** All 1785 structures are dead: 1105 coset-Parseval, 0 Theorem α, 0 Lemma R, 680 Q-level, 0 deeper level.
- **947, k6: closed [P, M].** All 18 structures are dead: 0 coset-Parseval, 0 Theorem α, 0 Lemma R, 18 Q-level, 0 deeper level.
- **971, k6: effectively closed.** 3 of 31 structures pass the gates run, all with log2 E_Q < −20 (best -368.4).
<!-- /CLOSED-FAMILIES -->

### 5.2 Rich structures killed by an exact gate

| order | family | structure | r, e | log2 E_Q [E] | killed at | check |
|---|---|---|---|---|---|---|
| 669 | ℤ₃ | ℤ₂ × ℤ₉ × ℤ₃₇ and ℤ₂ × ℤ₃² × ℤ₃₇, S of order 9 (4 structures) | 82, 41 | 84.7–84.8 | G/N = ℤ₆ (N of order 111): no solution | `bt_ladder.py`; independent CP-SAT `bt_xcheck_level.py`: INFEASIBLE in 66.9 s (`logs/bt/xcheck/x669_S9_N111.log`) |
| 669 | ℤ₃ | S of order 12 (4 structures) | 64, 35 | 45.4–45.6 | G/N = ℤ₆ and G/N = ℤ₉ (N of orders 111, 74) | `bt_level1.py`; CP-SAT INFEASIBLE at both levels (`x669_S12_N111.log`, `x669_S12_N74.log`, 165.6 s) |
| 669 | ℤ₃ | S of order 12 (3 structures) | 66, 37 | 34.2–34.6 | G/N = ℤ₆ | `bt_ladder.py` |
| 669 | ℤ₃ | S of order 18 (4 structures) | 46, 25 | 26.2–26.4 | G/N of order 9 (N of order 74) | `bt_gates_run.py` |
| 669 | ℤ₃ | S of order 18 (3 structures) | 60, 31 | 8.2–9.4 | N of order 74 | `bt_gates_run.py` |
| 669 | ℤ₃ | S of order 12 (2 structures) | 72, 39 | 32.1–32.2 | G/N = ℤ₆ (N of order 111), after passing N of order 6 with 6.9·10⁹ weighted solutions | `bt_gates_run.py` (`logs/bt/gates_run_extra.jsonl`) |
| 971 | k4 | ℤ₂ × ℤ₄ × ℤ₁₁², S of order 10 (split torus), 32 structures | 106, 56–59 | 83.3–98.7 | G/N = ℤ₂ × ℤ₄ (N = ℤ₁₁², order 121): no solution, about 3 s each | `bt_gates_run.py` (`logs/bt/gates_run_971b.log`) |
| 971 | k4 | S of order 10, r = 112, 16 structures | 112, 57–59 | 52.5–75.2 | N of order 22 (one of the two axis subgroups) | `bt_gates_run.py` |
| 839 | k4 | ℤ₂² × ℤ₁₁ × ℤ₁₉, S = ⟨(α, ×3, ×4)⟩ (order 45) | 30, 15 | 5.2 | exhaustive chain (160 level-A solutions, 0 cores); the ℤ₁₉-level gate is empty on its own | `bt_decide.py`, `bt_gates_run.py` |
| 839 | k4 | S of order 45 | 36, 19 | −8.5 | N of order 19 | `bt_gates_run.py` |
| 839 | k4 | nine more structures, r = 18–33 | | ≤ 5 | exhaustive chain | `bt_decide.py` (`logs/bt/decide_run1.log`) |
| 751 | k8 | ℤ₈ × ℤ₃ × ℤ₃₁, S of order 15 | 72 | — | Q-level | `bt_census2.py` |
| 863 | k4 | ℤ₂² × ℤ₅ × ℤ₄₃, S = ⟨(α, ×u, ×v)⟩ coupled through ℤ₄₃ (×9, ×38, ×10; orders 21, 42) | 20–50 | — | Q-level | `bt_census_run.py` (`logs/bt/cand_handpicked.json`) |

**All exact kills of structures with log2 E_Q ≥ 0** (generated from the logs by `theory/bt_summary.py`; levels
named by the order of N; "chain" = the exhaustive chain):

<!-- MORE-KILLS -->
| n | family | G = P × O | S order | r | structures | log2 E_Q [E] | killed at | script |
|---|---|---|---|---|---|---|---|---|
| 971 | k4 | [2, 4] × [11, 11] | 10 | 106 | 24 | 83.2–98.7 | N of order 121 (G/N of order 8) | `bt_gates_run.py` |
| 971 | k4 | [2, 2, 2] × [11, 11] | 10 | 106 | 8 | 90.8–91.6 | N of order 121 (G/N of order 8) | `bt_gates_run.py` |
| 669 | z3all | [2] × [9, 37] | 9 | 82 | 1 | 84.8 | N of order 111 (G/N of order 6) | `bt_level1.py` |
| 669 | z3all | [2] × [9, 37] | 9 | 82 | 1 | 84.8 | N of order 111 (G/N of order 6) | `bt_ladder.py` |
| 669 | z3all | [2] × [3, 3, 37] | 9 | 82 | 2 | 84.7–84.7 | N of order 111 (G/N of order 6) | `bt_ladder.py` |
| 971 | k4 | [2, 4] × [11, 11] | 10 | 112 | 8 | 74.3–75.2 | N of order 22 (G/N of order 44) | `bt_gates_run.py` |
| 971 | k4 | [2, 2, 2] × [11, 11] | 10 | 112 | 8 | 52.5–53.1 | N of order 22 (G/N of order 44) | `bt_gates_run.py` |
| 669 | z3all | [2] × [9, 37] | 12 | 64 | 1 | 45.6 | N of order 111 (G/N of order 6) | `bt_level1.py` |
| 669 | z3all | [2] × [3, 3, 37] | 12 | 64 | 2 | 45.4–45.4 | N of order 111 (G/N of order 6) | `bt_ladder.py` |
| 669 | z3all | [2] × [9, 37] | 12 | 64 | 1 | 45.4 | N of order 111 (G/N of order 6) | `bt_ladder.py` |
| 669 | z3all | [2] × [3, 3, 37] | 12 | 66 | 3 | 34.2–34.6 | N of order 111 (G/N of order 6) | `bt_ladder.py` |
| 669 | z3all | [2] × [9, 37] | 12 | 72 | 1 | 32.2 | N of order 111 (G/N of order 6) | `bt_gates_run.py` |
| 669 | z3all | [2] × [3, 3, 37] | 12 | 72 | 1 | 32.1 | N of order 111 (G/N of order 6) | `bt_gates_run.py` |
| 669 | z3all | [2] × [9, 37] | 18 | 46 | 2 | 26.2–26.4 | N of order 74 (G/N of order 9) | `bt_gates_run.py` |
| 669 | z3all | [2] × [3, 3, 37] | 18 | 46 | 2 | 26.3–26.4 | N of order 74 (G/N of order 9) | `bt_gates_run.py` |
| 863 | k4 | [2, 2] × [5, 43] | 12 | 91 | 1 | 24.0 | N of order 43 (G/N of order 20) | `bt_gates_run.py` |
| 839 | k8 | [2, 2, 16] × [13] | 12 | 82 | 1 | 16.9 | N of order 8 (G/N of order 104) | `bt_gates_run.py` |
| 669 | z3all | [2] × [3, 3, 37] | 18 | 54 | 1 | 9.8 | N of order 111 (G/N of order 6) | `bt_gates_run.py` |
| 669 | z3all | [2] × [9, 37] | 18 | 54 | 1 | 9.6 | N of order 111 (G/N of order 6) | `bt_gates_run.py` |
| 669 | z3all | [2] × [3, 3, 37] | 18 | 60 | 3 | 8.2–9.3 | N of order 74 (G/N of order 9) | `bt_gates_run.py` |
| 839 | k4 | [2, 2] × [11, 19] | 45 | 30 | 1 | 5.2 | chain | `bt_decide.py (exhaustive chain)` |
<!-- /MORE-KILLS -->

### 5.3 What the screens [P] remove

<!-- SCREEN-COUNTS -->
Over 81151 structures: coset-Parseval 61344, Theorem α 696, Lemma R 793 (the screens [P]); Q-level empty 9250, deeper level or chain empty 1541 [M]; prefilter only 1872; alive 5645.
<!-- /SCREEN-COUNTS -->


---

## 6. Calibration of the stratified counts

**Estimators.**
- **log2 E_Q, the Q-level stratified count.** It sums, over the exact Q-level solutions (every row pattern, weighted
  by the number of ordered signed row vectors), the product of the four Q-fibre sizes. Each term is multiplied by the
  lattice-Gaussian probability of the remaining class equations, `theory/q523_gauss.gauss_bits`. That probability uses
  blocks sampled uniformly in their Q-fibres and is **read at the true target** (no reading at the mean). The script is
  `theory/bt_census_run.evaluate` / `bt_census2.evaluate`.
- **log2 E_A** is the same sum taken over the exact solutions of a finer level G/N (`theory/bt_level1.py`).
- Both estimate the number of ordered quadruples of S-invariant sequences, over all signs and orders, that satisfy
  every class equation.

### 6.1 Against exact counts at small orders [M counts]

`theory/bt_calib_small.py` evaluates every k = 4 structure of the census enumeration with r ≤ 18 at the orders n ≡ 3
(mod 4), n ≤ 79, as far as it got. It counts their cores exactly with `theory/bt_brute_count.py`, which shares no code
with the census: full enumeration of the 2^r orbit signs per block, then a 2 + 2 join on the distinct PAF vectors with
two independent 64-bit hashes, which always agreed. The tables are made by `theory/bt_calib_table.py` from
`logs/bt/calib/calib_small.jsonl`.

<!-- CALIB-TABLE -->
(n = 7 left out: K = 0 there, so the Q-level solutions are the cores and E_Q is exact; 6 structures.)
structures: 561 at orders [11, 15, 19, 23, 27, 31]; killed by the gates: 466 (coset-Parseval 312, Theorem alpha 43, Q-level 111); killed but with cores: 0; hash disagreements: 0

| log2 E_Q | structures passing the gates | with cores | log2 N − log2 E_Q (median; min, max) over those with cores |
|---|---|---|---|
| none (−∞) | 11 | 0 | — |
| < 0 | 2 | 0 | — |
| [0, 10) | 8 | 0 | — |
| [10, 20) | 36 | 6 | -0.49; -0.80, +0.03 |
| [20, 30) | 24 | 10 | -0.27; -3.29, +0.56 |
| [30, 40) | 12 | 12 | -0.81; -3.90, -0.14 |
| [40, 50) | 2 | 2 | -0.28; -0.35, -0.28 |

S = 1 (no symmetry): n = 11 [8]x[]: log2 E_Q 21.76, log2 N 21.58, n = 11 [2, 4]x[]: log2 E_Q 22.25, log2 N 22.17, n = 11 [2, 4]x[]: log2 E_Q 20.42, log2 N 18.58, n = 11 [2, 2, 2]x[]: log2 E_Q 19.61, log2 N None, n = 15 [4]x[3]: log2 E_Q 31.6, log2 N 31.46, n = 15 [2, 2]x[3]: log2 E_Q 32.41, log2 N 32.22, n = 19 [16]x[]: log2 E_Q 39.94, log2 N 39.74, n = 19 [2, 8]x[]: log2 E_Q 40.97, log2 N 40.69, n = 19 [2, 8]x[]: log2 E_Q 38.38, log2 N 37.68, n = 19 [2, 8]x[]: log2 E_Q 38.32, log2 N 37.51, n = 19 [4, 4]x[]: log2 E_Q 38.37, log2 N 37.49, n = 19 [4, 4]x[]: log2 E_Q 40.92, log2 N 40.57, n = 19 [2, 2, 4]x[]: log2 E_Q 37.8, log2 N 36.3, n = 19 [2, 2, 4]x[]: log2 E_Q 37.79, log2 N 36.27, n = 19 [2, 2, 4]x[]: log2 E_Q 35.2, log2 N 33.3, n = 19 [2, 2, 2, 2]x[]: log2 E_Q 31.48, log2 N 27.58
<!-- /CALIB-TABLE -->

Reading:
1. **The gates are sound.** No structure that coset-Parseval, Theorem α or the Q-level gate kills has a core.
2. **When cores exist, E_Q is close.** Over the structures with cores, log2 N − log2 E_Q lies between −3.9 and
   +0.6 bits, with per-bin medians between −0.8 and −0.3, so E_Q slightly overestimates.
   - With S = 1 the overestimate is 0.1–0.3 bits for cyclic 2-parts and the split ones at n = 15.
   - It grows with the elementary abelian part of P: up to 1.9 bits for ℤ₂ × ℤ₂ × ℤ₄, and 3.9 bits for ℤ₂⁴
     (n = 19).
   - ℤ₂³ at n = 11 has no core at all, against E_Q = 2^19.6.
3. **Zero counts are common below about 2^30.** Many structures with E_Q between 2^10 and 2^30 have no core at all.
   Above 2^30 every structure seen so far had cores. A positive E_Q is necessary, not sufficient: the zeros are
   arithmetic (the target is not reachable), which the Gaussian cannot see. The n = 63 cases of
   `theory/bt_calib.py` (Q-level estimates 2^4.9–2^19.0, exhaustive chain 0, `logs/bt/calib/calib_k4.jsonl`) are the
   same phenomenon. Three of them are Theorem α (with level-A estimates 2^12.0); the three with S of order 4 and
   r = 19–20 have no known reason. One of the α cases was also confirmed by the brute-force counter.

### 6.2 At the known successes [M levels, E counts]

| control | r, e | log2 E_Q | exact data |
|---|---|---|---|
| 3964 (Kumar) | 46, 23 | 30.0 | 2 Q-level solutions (two sorted row patterns); 304 ℤ₁₉-level solutions; in the planted ℤ₁₉ solution 160 ℤ₁₃ combinations, of which the one run (the planted one) holds 3 cores. A crude extrapolation (3 per combination, every branch alike) gives 2^26–2^27 ordered cores. |
| 2636 (Kumar) | 54, 27 | 53.1 | 4 Q-level solutions; 192 level-A solutions at Kumar's row pattern |
| 2148 (Kumar) | 54, 27 | 62.1 | 27 Q-level solutions (one row pattern, (0, 2, 2, 2)); Kumar: 1 completion of 284,292 U-solutions, at his predicted rate |

So the successes sit at r = 46–54 with E_Q = 2^30–2^62. The 3964 estimate is within 3–4 bits of the crude
extrapolation, which comes from the planted branch and so is biased upward.

### 6.3 At the open orders [M]

Every exactly decided open-order structure with a positive estimate is a zero:
- 839, S45 r30: estimate 2^5.2;
- 669 ℤ₃ hosts:
  - S18 r46: 2^26;
  - S12 r64: 2^45;
  - S9 r82: 2^85;
- and the others in §5.2.

The 669 ℤ₃ zeros sit far above the small-order threshold of §6.1, at the level G/N = ℤ₆. That level is the first
place where the "two spikes" of the ℤ₃ family (Σ r² = 12) meet the order-2 factor of K. Estimates are therefore
used here only to rank structures that pass every affordable exact gate, never to declare a structure live.


---

## 7. The border problem

- **Reduction [P].** This is Kumar §18, as implemented in `theory/bt_border.py`. Let the rows of R̂ span V_T (the
  4t-dimensional space of block vectors with Fourier support in T). Then Ĉ = −Mc R̂ᵀ X/(4t), and Ĉ Ĉᵀ = I₄ ⊗ Γ_Q
  holds automatically, because Mc preserves V_T and acts there with Mc Mcᵀ = 4t. So a border exists for a frame R̂
  iff S_R = {ε ∈ {±1}^{4t} : Mc R̂ᵀ ε ∈ {±4t}^{4q}} contains 4t mutually orthogonal vectors (the columns of X).
- **Frames.** For (q, t) = (4, 3) the frames are Sylvester frames H₄ ⊗ H₄[1:] with the columns of each block permuted.
  For (3, 3) they are Paley H₁₂ with signed column permutations.
- **Status [M]** (`theory/bt_border_census.py`, `logs/bt/border_census.jsonl`). A border depends only on the
  Q-compressions, so checking every Q-level class of a structure settles the border problem for every core it could
  have.
  - The census covered the top 25 survivors per order and family by E_Q: <!-- BORDER-COUNTS -->86 structures (669 z3all: 25, 839 k4: 11, 863 k4: 17, 863 k4one: 8, 955 k4: 25), 850 Q-level classes, 850 with a border<!-- /BORDER-COUNTS -->.
  - Every class has a border. This matches the earlier 2048 of 2048 at Q = ℤ₂²
    (`report/alpoge_design_answer.md` §3) and Kumar's 864 of 864 at Q = ℤ₃.
  - The controls 3964, 2636 and 2148 were bordered and assembled into Hadamard matrices (§1).
- **k = 8 (q = 8, t = 7)** (`theory/bt_border_k8.py`, `logs/bt/border_k8.jsonl`; S_R by a meet in the middle over the
  four blocks, `bt_border.admissible_balanced`, about 0.6 s per frame).
  - *Reduction [P].* The Q-level classes of the 751 survivors are monomial: z_k = c_k + a_k δ_{p_k} with
    {|a_k|} = {2, 2, 2, 4}. Translating, negating or permuting the blocks of a core gives a core. So a monomial class
    can be completed (after such a change of the core) iff the normal form (2, 2, 2, 4)·δ₀ has a border, and the
    answer depends only on Q.
  - *Kronecker frames fail.* The frame H₄ ⊗ H₈[1:] has no admissible vector at Q = ℤ₂³ [P]. The proof is in the
    script header: flatness would need the Walsh coefficients of the blocks to equal (2/7)·Cᵀe with C = H₄ Q_Sᵀ, and
    Cᵀe never has entries in {0, ±14, ±28}. Its block-permuted variants have none either [M]: 0 admissible vectors in
    40 frames at ℤ₂³ and ℤ₈, and in 300 frames at the two 751 classes (`logs/bt/border_k8_751.log`).
  - *Q = ℤ₂³: every monomial class has a border [M].* The c = 8 Singer-type template of `theory/ab_singer.py`
    (Alpöge's 1436 form, `logs/ab/singer/c8_control/`) is used in our layout as the frame R̂, with one block negated.
    Its admissible set has exactly 28 vectors, and they form the corner. The script checks X Xᵀ = 28 I, Ĉ ∈ {±1},
    Ĉ Ĉᵀ = I₄ ⊗ Γ_Q and X Ĉᵀ + R̂ Mcᵀ = 0.
  - *Q = ℤ₈ (the 751 survivors with log2 E_Q ≥ 0): open.* The template frames have no admissible vector, with either
    frame and every block sign, at the normal form and at both 751 classes [M].
    - On V_T the ℤ₈ array is Q_S ⊗ Π₊ + Q_Sᵀ ⊗ Π₋, where Π₊ and Π₋ project onto the balanced functions that are even
      or odd under negation. At ℤ₂³ it is Q_S ⊗ I.
    - No frame symmetry relates the two arrays [P]. Such a symmetry would make the product of the four blocks
      around every rectangle of the block array trivial. But a rectangle through one diagonal block (a translation)
      and three off-diagonal blocks (reflections g ↦ p − g) gives a reflection, which is never the identity on ℤ₈.
    - The same argument excludes a symmetry of order 7, for every placement of the spikes [P]. Its components would
      have to commute with that reflection, and no involution of 8 points commutes with a 7-cycle. So the
      Singer-type ansatz, which needs an order-7 symmetry, is not available.
    - The diagonal symmetries of the ℤ₈ array are the 96 permutations commuting with negation, against 8! at ℤ₂³
      [M].
    - A frame search adapted to that group was not built.
- **k = 6 (q = 6, t = 5)** (`theory/bt_border_k6.py`, `logs/bt/border_k6.jsonl`).
  - *Control [M].* Alpöge's 347 template has non-monomial compressions on Q = ℤ₂ × ℤ₃.
    - It is reproduced in its own layout (`ab_assemble.forced_corner`).
    - Either of its frames, used as R̂, gives |S_R| = 20 = 4t, which is the corner. This holds in its layout and in
      ours (Kumar's), after a CRT relabelling of ℤ₆.
  - *The Lemma-1 frames do not find borders that exist.* For the same compressions, H₂ ⊗ Paley H₁₂ (every choice
    of the row s) and 2,000 block-permuted variants give at most 6 admissible vectors in either layout.
  - *The k = 6 survivors: open.* We tested 12 structures at 599, 719, 863 and 971, with 6 Q-level classes each,
    against the 347 frames (every block sign) and 2,000 permuted variants. The largest admissible set has 6
    vectors, so no border was found.
  - Their classes are not monomial, so the Singer-type search does not apply either.
- **k4one (Q of order 4, T = {1} and two nontrivial characters, t = 3)** (`theory/bt_border_k4one.py`,
  `logs/bt/border_k4one.jsonl`) [M].
  - On one block, V_T is spanned by three orthogonal ±1 vectors c_k: the two T-characters and 1 for Q = ℤ₂²; for
    Q = ℤ₄, 1, (1, 1, −1, −1) and (1, −1, −1, 1). So a frame row is ±c_k on every block.
  - The frames H₄ ⊗ (c_k), with the same k in every block, give no border for any class of the 863 hosts.
  - The frames that let each block permute the c_k (6³ · 2⁴ candidates) give a border for every class tested:
    64 of 64, across the 8 hosts with log2 E_Q ≥ 20, including the two with Q = ℤ₄. Each border passes the exact
    checks.
- **The border problem is not what stops these structures.**
  - At k = 4, in the ℤ₃ family and in k4one, every tested class has a border.
  - At k = 8 the classes at ℤ₂³ have one. The ℤ₈ classes, and the k = 6 classes, are open: borders of this size
    exist (1388, 1436), but only template-type frames find them.
  - The cores of those structures are out of reach anyway: their next exact gates cost 2^50 pairs or more (§8.1).


---

## 8. Costs, the rich hosts, and why nothing reached production

**Unit of cost.** The exact 2 + 2 join at a level G/N costs L_a·L_d + L_b·L_c pair operations, summed over the
Q-level solutions; L_i is the size of block i's refined list (`theory/bt_cost.py` gives a bound, the product of the
per-cell list sizes before the character filter). Measured throughput of the bucketed join (`bt_chain.join4b`) on
this machine under load is about 3·10⁶ pairs per second per core. So 2^34 pairs is about an hour, 2^40 about four
days, and 2^50 about 10⁴ years.

### 8.1 The hosts with the largest estimates

<!-- RICH-TABLE -->
| n | family | S, r, e | log2 E_Q | exact levels passed | capped levels | cheapest unrun level (log2 pairs) |
|---|---|---|---|---|---|---|
| 599 | k4 | S2, 300, 224 | 76.71 | — | — | — |
| 643 | k4 | S8, 160, 86 | 116.51 | — | — | — |
| 643 | k4 | S8, 160, 86 | 115.18 | — | — | — |
| 643 | k4 | S8, 124, 65 | 95.4 | N80: 9534720, N40: 141103911168, N32: 66710016 | — | N16: 61.0 |
| 669 | z3all | S9, 90, 45 | 73.14 | N111: 130553856 | N74, N37, N6, N3, N2 | N6: 52.9 |
| 669 | z3all | S9, 90, 45 | 73.12 | N111: 130553856 | N74, N37, N6, N3, N2 | N6: 52.9 |
| 719 | k4 | S2, 360, 269 | 64.09 | — | — | — |
| 719 | k6 | S8, 112, 59 | 66.18 | — | — | N7: 58.3 |
| 739 | k4 | S6, 144, 95 | 53.64 | N92: 202515200 | — | N46: 45.5 |
| 751 | k4 | S8, 121, 63 | 73.45 | — | — | N11: 48.3 |
| 751 | k4 | S8, 121, 65 | 66.02 | — | — | N11: 48.3 |
| 751 | k4 | S8, 111, 65 | 49.36 | — | — | N11: 48.3 |
| 751 | k8 | S10, 105, 54 | 66.19 | — | N31, N3 | N31: 68.8 |
| 751 | k8 | S10, 105, 53 | 61.06 | — | N31, N3 | N31: 68.8 |
| 751 | k8 | S10, 91, 48 | 48.64 | — | N31, N3 | N31: 51.2 |
| 839 | k4 | S6, 156, 83 | 165.86 | — | N19, N11 | N19: 78.4 |
| 839 | k4 | S6, 154, 83 | 157.24 | — | N19, N11 | N11: 81.0 |
| 839 | k4 | S9, 110, 55 | 105.13 | — | N19, N11 | N11: 51.5 |
| 839 | k8 | S12, 105, 56 | 53.02 | — | — | — |
| 839 | k8 | S12, 105, 56 | 52.98 | — | — | — |
| 839 | k8 | S12, 104, 58 | 44.48 | — | — | N52: 45.4 |
| 863 | k4 | S7, 140, 71 | 136.89 | — | N43, N5 | N5: 106.7 |
| 863 | k4 | S6, 155, 85 | 135.2 | — | N43, N5 | N43: 98.0 |
| 863 | k4 | S7, 140, 70 | 133.35 | — | N43, N5 | N5: 106.7 |
| 863 | k4one | S14, 77, 41 | 47.85 | — | — | N43: 55.6 |
| 863 | k4one | S14, 77, 41 | 41.45 | — | — | N43: 55.9 |
| 863 | k4one | S14, 77, 41 | 41.45 | — | — | N43: 55.9 |
| 955 | k4 | S8, 147, 76 | 103.58 | N119: 252861382656 | N34, N17, N14, N7, N2 | N14: 50.2 |
| 955 | k4 | S8, 140, 75 | 101.0 | N119: 2211840 | N34, N17, N14, N7, N2 | N14: 48.1 |
| 955 | k4 | S8, 147, 77 | 96.99 | N119: 1104150528 | N34, N17, N14, N7, N2 | N14: 50.2 |
| 971 | k4 | S10, 112, 57 | 73.3 | N22: 786432 | — | N121: 49.5 |
| 971 | k4 | S10, 112, 57 | 73.3 | N22: 786432 | — | N121: 49.5 |
| 971 | k4 | S10, 112, 57 | 73.23 | N22: 786432 | — | N121: 49.5 |
<!-- /RICH-TABLE -->

Every structure with E_Q ≥ 2^40 at an open order has r ≥ 77. Its cheapest unrun exact gate costs 2^50 pairs or more,
or its last gate passed with so many solutions (5.5·10⁷ at 955, 2.5·10⁹ at 839) that the next level would have to
run inside each of them.

### 8.2 The nearest analogues of the successes

The three successes reproduced in §1 have r = 46–54 and E_Q = 2^30–2^62, and their chains end in fibres of about 10⁵
vectors (3964). At the open orders, the best live structures in and near that window are:

<!-- WINDOW-TABLE -->
| n | family | best alive with r ≤ 60: log2 E_Q (S, r) | best alive with r ≤ 80 | best alive overall |
|---|---|---|---|---|
| 599 | k4 | — | — | 76.71 (S2, 300) |
| 599 | k4one | — | — | — |
| 599 | k6 | — | — | -49.81 (S10, 84) |
| 599 | k8 | -594.85 (S21, 52) | -594.85 (S21, 52) | -18.54 (S9, 104) |
| 643 | k4 | — | -562.3 (S24, 68) | 116.51 (S8, 160) |
| 669 | prod12 | -14.28 (S20, 57) | -14.28 (S20, 57) | -14.28 (S20, 57) |
| 669 | z3all | — | — | 73.14 (S9, 90) |
| 719 | k4 | — | — | 64.09 (S2, 360) |
| 719 | k4one | — | — | — |
| 719 | k6 | -1.93 (S24, 48) | -1.93 (S24, 48) | 66.18 (S8, 112) |
| 719 | k8 | — | — | — |
| 739 | k4 | — | — | 53.64 (S6, 144) |
| 751 | k4 | -15.38 (S20, 57) | -15.38 (S20, 57) | 73.45 (S8, 121) |
| 751 | k4one | -22.46 (S20, 57) | -22.46 (S20, 57) | -22.46 (S20, 57) |
| 751 | k4split | -33.19 (S80, 27/36/36/27) | -33.19 (S80, 27/36/36/27) | -33.19 (S80, 27/36/36/27) |
| 751 | k8 | -571.37 (S42, 44) | -571.37 (S42, 44) | 66.19 (S10, 105) |
| 839 | k4 | 11.17 (S18, 60) | 44.86 (S15, 78) | 165.86 (S6, 156) |
| 839 | k4one | -26.38 (S30, 45) | -26.38 (S30, 45) | 19.59 (S15, 84) |
| 839 | k4split | -16.41 (S90, 16/30/30/16) | -16.41 (S90, 16/30/30/16) | -16.41 (S90, 16/30/30/16) |
| 839 | k6 | — | — | — |
| 839 | k8 | -53.82 (S24, 60) | -46.96 (S16, 80) | 53.02 (S12, 105) |
| 859 | k4 | — | — | — |
| 859 | k4one | — | — | — |
| 863 | k4 | — | 45.11 (S14, 77) | 136.89 (S7, 140) |
| 863 | k4one | — | 47.85 (S14, 77) | 47.85 (S14, 77) |
| 863 | k6 | -14.56 (S30, 60) | -11.31 (S20, 66) | 3.28 (S15, 90) |
| 863 | k8 | — | — | — |
| 907 | k4 | — | — | — |
| 907 | k4one | — | — | — |
| 947 | k4 | — | — | -113.84 (S12, 120) |
| 947 | k4one | — | — | — |
| 947 | k6 | — | — | — |
| 955 | k4 | -10.95 (S24, 60) | -4.22 (S24, 63) | 103.58 (S8, 147) |
| 955 | k4one | -14.51 (S24, 58) | -3.7 (S24, 63) | -3.7 (S24, 63) |
| 971 | k4 | -71.88 (S50, 54) | -16.12 (S25, 72) | 73.3 (S10, 112) |
| 971 | k6 | — | — | -368.36 (S12, 118) |
<!-- /WINDOW-TABLE -->

**Timed exact gates on the window [M]** (`theory/bt_gates_levels.py --tmax 1800 --maxlen 3e7`,
`logs/bt/gates_run_window.jsonl`). Every alive structure with r ≤ 80 and log2 E_Q ≥ −5 (11 structures) was run at
every invariant level inside K. None died; each level passed or hit a cap:
- 719 k6, S24 r48 (2^−1.9): G/N₁₇ passes with 24,576 weighted solutions; G/N₇ exceeds 30 minutes.
- 839 S18 r60 (2^11.2): the G/N₁₉ list exceeds 3·10⁷ vectors. G/N₁₁ needs 4.9·10¹¹ join pairs for the first of its
  20 Q-level solutions alone.
- 955 k4 and k4one, S24 r63 (2^−4.2 and 2^−3.7): G/N₃₄ passes with 4.1·10⁶. G/N₁₁₉ exceeds 30 minutes, G/N₁₄ needs
  2.7·10¹⁰ pairs, and the other levels exceed the list cap.
- 863 S14 r77 (six structures, 2^27–2^48): G/N₄₃ exceeds the list cap or 30 minutes; the G/N₅ lists have 8·10¹⁰
  vectors.
- 839 S15 r78 (2^44.9): G/N₁₉ passes with 2.5·10⁹; G/N₁₁ exceeds the list cap.

### 8.3 The best candidate measured to the end of what the tools can do: 839, S15, r = 78

- **The structure.** ℤ₂² × ℤ₁₁ × ℤ₁₉, S of order 15, 78 orbits per block, 39 equations; log2 E_Q = 44.9.
- **Level A (ℤ₁₉, G/N = ℤ₂² × ℤ₁₁, 6 coordinates) [M].** 2.46·10⁹ weighted solutions, total multiplicity 2^267.3,
  and log2 E_A = 45.4 [E]. So each level-A solution carries about 2^14 cores in about 2^236 quadruples.
- **Level B (ℤ₁₁) inside one level-A fibre.** The fibre DP over the six level-A cells exceeds 3·10⁶ states and the
  memory of the machine at the first cells (`theory/bt_chain.fibre_dp`, probe on the first level-A solution). Going
  straight from Q to level B costs 2^116 pairs.
- **No other route.** The group has no further invariant subgroup inside K.
- **Same picture at 863** (S14, r = 77, log2 E_Q = 45.1). Its only finer levels, N43 and N5, cost 2^52.8 and 2^79.6
  pairs from Q.

### 8.4 Local search does not substitute

- `theory/bt_anneal.py` is simulated annealing over orbit signs. Its energy is the full PAF residual, scored
  incrementally with exact bookkeeping checked against a direct recomputation.
- It never reached a core of the 3964 control, a structure with at least 3 known cores and log2 E_Q = 30:
  - 108 restarts of 2·10⁶ flips and 33 restarts of 2·10⁷ flips;
  - best energy 11,472, that is one or two classes off by 16;
  - logs `logs/bt/anneal/`.
- Rich open-order structures are larger by every measure, so annealing was not run on them.

### 8.5 What it would take

The structures that the gates leave alive with large estimates are all "rich and flat": many orbits, no deep chain
of invariant subgroups, and level lists far beyond 10⁷. A construction at these orders through this route needs one
of the following:
- a structure with r ≈ 50–60 and a calibrated E_Q above 2^30. The census found none: the best such structure is
  839, S18 r60, at 2^11.2, and below 2^30 zero counts are common (§6.1);
- or a staged solver for the rich hosts that never lists a level, such as a Gauss-period slice along a large prime
  axis (ℤ₄₃ at 863, ℤ₁₉ at 839) with a twisted S/D split, as in our 571 and 523 stage-2 solvers. None was built
  here. §8.1 gives each host's measured gate status and next-gate cost as the starting point.

**Kumar's ℤ₂ U/W split (§19, the route that made 2148) does not reach these hosts [E]** (`theory/bt_uw_scan.py`,
`logs/bt/uw_scan.jsonl`).
- *The split.* Take an S-fixed involution z. U = (x + x(· + z))/2 is the compression to G/⟨z⟩, and W carries one free
  sign on each zero cell of U. The class equations split into U-equations (the characters trivial on z) and
  W-equations.
- *Zero points [P].* The lag-0 equation of the compression fixes the number of zero points of U over the four blocks:
  |G| + Γ(z + K)/4, which is n when z ∈ K. So a U-solution leaves W-lists of about 2^zc sign vectors per block, with zc
  ≈ r/4 zero cells. Its exact W-join costs about 2^(2zc) pairs before the coset-sum constraints.
- *Completion rate.* The scan estimates the completions per U-solution as 2^(4zc) · P(W | U). P(W | U) comes from the
  lattice Gaussian of the class sums, taken on all equations and on the U-equations alone.
- *Calibration on 2148.* The scan gives zc = 12.3 (Kumar: 12 zero cells per block, about 200 W-vectors after the
  constraints) and 2^−14 completions per U-solution. Kumar observed 2^−18 (one in 284,292). The estimate is 4 bits
  optimistic.

At our hosts:

| hosts | zc (W-list size per block, log2) | completions per U-solution [E] |
|---|---|---|
| 863 k4 r140; 955 r147; 971 r112; 669 S9 r90 | 31; 31; 25.5; 19 | 2^−61; 2^−60; 2^−51; 2^−45 |
| 751 k4 r121; 719 k6 r112; 739 r144; 863 k4 r155; 955 r140 | 24; 23.4; 32.8; 36.1; 30.3 | 2^−50 to 2^−57; 2^−70; 2^−95; 2^−38; 2^−31 |
| 751 k8 r105 | 23 | 2^−13 to 2^−19 |
| 643 r122; 599 r300; 719 r360 | 26.7; 75; 90 | 2^−6; 2^8; 2^−1.5 |
| 863 r77; 955 r63; 719 k6 r48 (the window) | 16.6; 11.4; 8.6 | 2^−25; 2^−41; 2^−52 |

- In the first three rows a core needs 2^13 to 2^95 U-solutions, each with an exact W-join of 2^38 pairs or more.
- The hosts with more completions (third row) have W-lists of 2^27 or more per block.
- The 839 hosts (r = 60 to 156) have no S-fixed involution: S permutes the three involutions of ℤ₂² cyclically.


---

## 9. Reproduce

Python 3 with numpy, numba, ortools. Kumar's repository is needed only for the controls; it is cloned read-only at
`/home/user/vatskishan03/hadamard`.

```
# controls (section 1)
python3 theory/bt_controls.py --kumar /home/user/vatskishan03/hadamard      # 3964, 1916, 2636 -> logs/bt/controls.log
python3 theory/bt_control2148.py /home/user/vatskishan03/hadamard/hadamard/2148/certificate.json

# census (section 4); every runner resumes from its output file
python3 theory/bt_census_run.py logs/bt/census/k4_P4_w0.jsonl --ns 839,863,955 --rmin 8 --rmax 160   # k = 4
python3 theory/bt_census_run.py logs/bt/census/k4_P4_w8.jsonl --ns 643 --rmin 8 --rmax 160 --cap 4096   # 643: 2-parts with <= 4,096 endomorphisms
python3 theory/bt_census_run.py logs/bt/census/k4_P4_w4.jsonl --ns 739 --rmin 8 --rmax 160 --cap 16384 --min-cap 1024
python3 theory/bt_census_run.py logs/bt/census/k4_P4_w4.jsonl --ns 643 --rmin 8 --rmax 160 --cap 16384 --min-cap 4096 --part 0/3   # and 1/3, 2/3
python3 theory/bt_census2.py logs/bt/census2/k8_w0.jsonl --fam k8 --ns 751,719,863 --rmin 12 --rmax 110
python3 theory/bt_census2.py logs/bt/census2/z3all_669.jsonl --fam z3all --ns 669
python3 theory/bt_census2.py logs/bt/census2/prod12_669.jsonl --fam prod12 --ns 669 --rmin 12 --rmax 100
sh logs/bt/census2/queue2.sh 0 1                                              # 971, 947, k6, k8 (599, 839), k4one
python3 theory/bt_census2.py logs/bt/census2/k8c_839.jsonl --fam k8 --ns 839 --rmin 12 --rmax 110 --cap 16384 --part 0/2   # and 1/2
sh logs/bt/restart.sh                                                         # every live worker, after a container restart
python3 theory/bt_census_split.py logs/bt/census/split_k4.jsonl logs/bt/census/k4_P4_w0.jsonl --rmin 16 --rmax 50
python3 theory/bt_close.py logs/bt/census2/close.jsonl logs/bt/census2/k8_w0.jsonl     # exact gates for prefiltered

# screens [P] and deeper gates [M]
python3 theory/bt_thm_alpha.py logs/bt/census/thm_alpha.jsonl CENSUS.jsonl[,...]
python3 theory/bt_thm_rigid.py logs/bt/census/thm_rigid.jsonl CENSUS.jsonl[,...]
python3 theory/bt_cost.py logs/bt/cost_next_gate.jsonl CENSUS.jsonl[,...] --top 25 --min-logE -10
python3 theory/bt_gates_run.py logs/bt/cost_next_gate.jsonl logs/bt/gates_run.jsonl CENSUS.jsonl[,...] --cap 34
python3 theory/bt_ladder.py CENSUS.jsonl logs/bt/ladder_run2.jsonl --maxlen 2e6
python3 theory/bt_level1.py CENSUS.jsonl KEY --out logs/bt/level1_KEY.jsonl
python3 theory/bt_decide.py CENSUS.jsonl logs/bt/decide --rmax 45
python3 theory/bt_xcheck_level.py logs/bt/xcheck/x669_S9_N111.json        # independent CP-SAT
python3 theory/bt_gates_levels.py logs/bt/gates_run_levels.jsonl CENSUS.jsonl[,...] --levels 74,111 --fam 669:z3all --tmax 600
sh logs/bt/closure.sh                                                       # timed gates at every level, effectively closed families
python3 theory/bt_rerun_capped.py logs/bt/census2/prod12_669.jsonl --maxlen 6e6 --maxout 1e6   # capped Q-levels

# borders, calibration, local search
python3 theory/bt_border_census.py logs/bt/border_census.jsonl CENSUS.jsonl[,...] --top 25 --min-logE -10
python3 theory/bt_border_k8.py --random 40                                   # k = 8: Kronecker, Singer template, Z2^3 vs Z8
python3 theory/bt_border_k6.py --frames 3000 --survivors 3                   # k = 6: 347 control, Lemma-1 frames, survivors
python3 theory/bt_border_k4one.py --min-logE 20                             # k4one: frames with block-wise character permutations
python3 theory/bt_calib_small.py logs/bt/calib/calib_small.jsonl --ns 7,11,15,19,23,27,31 --rmax 18
python3 theory/bt_brute_count.py RECORD.json                                # independent exact count
python3 theory/bt_anneal.py logs/bt/anneal/ctrl3964 --control 3964 --hours 0.05
python3 theory/bt_uw_scan.py --cert /home/user/vatskishan03/hadamard/hadamard/2148/certificate.json   # U/W split: 2148
python3 theory/bt_uw_scan.py --top 2 --min-logE 40; python3 theory/bt_uw_scan.py --keys @logs/bt/window_keys.txt

# tables of this report
python3 theory/bt_summary.py            # print
python3 theory/bt_summary.py --fill     # rewrite the generated blocks of this file
python3 theory/bt_calib_table.py
```

