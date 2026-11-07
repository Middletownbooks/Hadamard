# Prompt (OpenAI model, maximum reasoning): a genuinely new method for the open Hadamard orders

*Self-contained. You have access to the public repository `https://github.com/middletownbooks/hadamard`, and you may
push only to its branch `codex/order-3452-compression`. Claude Code instances will implement, calibrate and run what
you propose, with 4 cores, about 15 GB of RAM and days of compute. Your job is the idea, its analysis, and a plan
precise enough that they can test it without asking you anything.*

## 0. The framing

**Imagine it is August 2026.** You are the AI that a small team of mathematicians has asked to help find Hadamard
matrices at orders where none is known. That is the position of whatever assistance lay behind the matrices that
Alpöge, Voinov and Reynolds-Haertle announced on 12 August 2026, at the twelve orders then open below 2000. Their
announcement gave the matrices but not the method.

**Start from first principles.** What would you try that the field has not? You now also know:
- what has worked since (§1, from the November 2026 revision of the survey by Egan, Kotsireas, London and Ó Catháin,
  arXiv:2610.27175v2, below "EKLO");
- what has been tried exhaustively in this repository and failed (§2).

**The bar is a new idea, not a re-run.** Re-running an existing search with different parameters, or a refinement
with no new ingredient, does not meet it. An idea counts if it changes what is searched (new arrays, groups, ansätze
or ingredients), or how (new algorithms with a different cost exponent), or replaces search by algebra. It must come
with an honest argument for why it could succeed where the methods in §2 failed.

## 1. What EKLO (v2, 3 November 2026) say about the methods that worked

**Arrays.** Each is over an abelian group G of order v, with group-developed ±1 matrices X, Y, Z, W, and
R = the matrix of x ↦ −x.
- **Theorem 3 (Goethals–Seidel, b = 0):** Σ PAF = 0 off 0 gives H(4v).
- **Theorem 5 (b = 4, Legendre quadruples):** m even, Σ PAF = −4 off 0, row sums (2, 0, 0, 0), and a border of width
  4 give H(4m + 4).
- **Theorem 7 (b = 4(k − 1)):**
  - S is a proper subgroup of index k.
  - Σ_f PAF_f(d) is −b on S \ {0} and +4 off S.
  - A border (K, E, F) is constant on the cosets of S, with K a Hadamard matrix of order b, EEᵀ = 4kI,
    FFᵀ = 4kI − 4(I₄ ⊗ J_k) and MΠEᵀ + ΠFKᵀ = 0.
  - This gives H(4v + b).
- **Character form:** Σ_f |f̂(χ)|² = n = 4v + b for every χ nontrivial on S.

**The twelve 2026 matrices** (EKLO Table 1; all n = 4p with p prime, p ≡ 3 mod 4):

| n | b | G | U (multiplier group), orbits r |
|---|---|---|---|
| 892, 1132, 1244, 1948, 1964 | 0 | ℤ_p, p = 223, 283, 311, 487, 491 | \|U\| = 3, 3, 5, 9, 7; r = 75, 95, 63, 55, 71; no generator symmetric |
| 668, 716 | 4 | ℤ₁₆₆, ℤ₁₇₈ | **no multiplier**: all 4m signs free |
| 1676, 1772 | 4 | ℤ₄₁₈, ℤ₄₄₂ | \|U\| = 15, 12; r = 42, 50 |
| 1916 | 12 | ℤ₂² × ℤ₁₁₉, S = ℤ₁₁₉ | X, W fixed by a group of order 12 (46 orbits); Y, Z by order 4 (128 orbits) |
| 1388 | 20 | ℤ₂ × ℤ₃ × ℤ₅₇, S = ℤ₅₇ | \|U\| = 6, r = 70 |
| 1436 | 28 | ℤ₂² × ℤ₈ × ℤ₁₁, S = 2ℤ₈ × ℤ₁₁ | QR of ℤ₁₁, \|U\| = 5, r = 96 |

**The flatness statistic** (EKLO Table 2).
- **Definition:** Φ_f is the mean of (4|f̂(χ)|²/n − 1)² over the characters nontrivial on S. A random sequence has
  Φ ≈ 1. A typical solution under an exponential model has Φ ≈ 0.6, or about 0.3 if |f̂|² < n/2 everywhere.
- **The b = 0 matrices are extremely flat.** Φ = 0.08–0.14 for every generator, with |f̂(χ)|² < n/2 at every
  character. Under the model, at most 2^−57 (1948) to 2^−141 (1132) of the solutions in the same multiplier slice are
  that flat. EKLO infer that the search was restricted to nearly flat spectra.
- **The b ≥ 4 matrices are typical** (Φ = 0.27–0.93), including 668 and 716, which have no multiplier. There is no
  evidence of either restriction there.

**Matrices found since** (EKLO §4, with Kumar's and this repository's):
- **Theorem 3 with large multiplier groups:** 3732, 4492, 5932, 5020 (|U| = 31, 33, 39, 50); 3676 over ℤ₉₁₉ (also a
  skew one by Schumer); Kumar's 2284 over ℤ₅₇₁.
- **Theorem 5:** this repository's 2092, 2292 and 2588, from Legendre quadruples of lengths 522, 572 and 646 with
  multiplier slices and staged (compressed) solvers.
- **Theorem 7:**
  - 5788 and 8476 over ℤ₂² × F_{q²}, S = F_{q²}, q = 19 and 23;
  - Kumar's 2636 and 3964 (k = 4, over ℤ₂⁴ × ℤ₄₁ and ℤ₂² × ℤ₂₄₇);
  - Kumar's 2148 (border of width 12 over ℤ₃ × ℤ₁₇₈, sum −12 on S \ 0 and 0 off S);
  - Kumar's 3356 = 4·839 (index 8 over ℤ₄³ × ℤ₁₃).
- **Recursive and other constructions:**
  - T-matrices of orders 151 and 163 (multiplier subgroups of orders 5 and 9), and Kumar's 97, 109 and 197;
  - a symmetric BH(86, 4), giving WT(173) and six orders;
  - Kumar's Paley-blocks × T-matrices and Miyamoto-blocks × T-matrices products (thirteen orders);
  - Miyamoto I from 2292 (9172).
- **Kumar's "complex Yamada–Kiyasu" route** settled 4n for n = 643, 859, 955, 1051 and 1699:
  - a quaternary conference matrix of order (n + 1)/2 (zero diagonal, entries ±1, ±i), found as a zero-diagonal
    Szöllősi array over a group of order about (n + 1)/4;
  - then complex Yamada–Kiyasu (n − 2 a prime ≡ 1 mod 4);
  - then Turyn doubling.

  See `hadamard/2572`, `3436`, `3820`, `4204` and `6796` in `https://github.com/vatskishan03/hadamard`, and its
  `FRONTIER.md`. The same route is open, not yet built, at **4n = 5164, 7516, 8332, 8860 and 9532**, which need
  quaternary conference matrices of orders **646, 940, 1042, 1108 and 1192**.

**EKLO Table 3: the smallest open ingredient orders.**
- Golay pairs: 106; TT(n): 48; BS(n + 1, n): n = 44; T-sequences: 97; T-matrices: 127;
- Williamson-type: 71; skew Hadamard: 356; symmetric conference matrices: 146 (118 reported but not public);
- complex Golay pairs: 29.

**EKLO Table 5: the open orders.** It lists 95 orders up to 10000. Since then 2572, 3356, 3436, 3820, 4204, 5912 and
6796 were settled by Kumar, and 6712, 8248 and 8696 shown known. The repository's current list is
`data/open_orders_lt10000.json` (key `open`, 85 orders). The smallest are:
- 4p for p = 599, 719, 739, 751, 863, 907, 947, 971;
- 4·3·223 = 2676.

Every prime in the table is ≡ 3 (mod 4). The shapes are 4p (most), 4·3p, 4·5p, 4·7p, 4·9·191, 4·11·179, 4·13p and
8p.

**Other recent work cited by EKLO:**
- **Williamson, Yacobi and Zinn-Justin (2026):** transformer-guided local search for Goethals–Seidel matrices up to
  order 252.
- **OpenAI (6 October 2026):** a model-generated, Lean-formalised manuscript claiming Ryser's circulant Hadamard
  conjecture.
- **Kotsireas et al.:** compression and "uncompression" problems for Legendre pairs.

## 2. What this repository has already tried (read before proposing; do not repeat)

Everything is in `report/`. These are the main conclusions, with their files.

**What is now known about the 2026 matrices** (`report/alpoge_methods_forensics.md`, `report/alpoge_gap_analysis.md`,
`forensics/`):
- The b = 0 cores have individually flat blocks.
- 668 and 716 hide a **half-shift template**. With m = 2h, σ_k(y) = x_k(y) x_k(y + h) is a two-step function on ℤ_h in
  every block, up to one or two exceptions (log₁₀ P ≈ −71 and −78). This decouples the problem exactly into a
  periodic problem and a negacyclic one on ℤ_h (§4.2 there).
- The originals are in `forensics/originals/`.

**Staged multiplier-slice searches** (compress along a subgroup, solve stage 1 exactly, lift through cyclotomic
classes). They built 2092, 2292 and 2588 (`report/paley_523_answer.md`, `report/paley_571_answer.md`,
`report/paley_647_answer.md`) and 3732 (`report/alpoge_design_answer.md` §4).
- **The cost model works:** log₂(cost per core) ≈ D_def + log₂(cost per stage-1 solution). It reproduces 523 and 933.
- **The budget** is about 2^46.

**Every open order screened** (`report/open_screen_answer.md`, `report/h4892_answer.md`,
`logs/open_screen/followup_20261106/`): the best staged slice per order and family (b = 0, b = 1, k4, ℤ_q,
Cooper–Wallis).
- **Over 20 slices with positive estimates were decided exactly. All are empty.**
- The three whose stage 1 was non-empty (at 4892 and 7572) died at the coarsest stage-2 level, the ℤ_q compression.
  Two fell to an exact energy identity.
- The lesson:
  - slices cheap enough to search are small, and small systems are arithmetically rigid;
  - Gaussian estimates cannot see this;
  - slices rich enough to have solutions are too expensive for exact staged methods.

**The k4 family (Theorem 7, k = 4)** (`report/b3_push_answer.md`, `report/h863_lift_answer.md`,
`report/blockswap_answer.md`, `report/zq_census_answer.md`).
- **A census of 171,001 structures.** Rich hosts survive, e.g.:
  - 863: ℤ₂² × ℤ₅ × ℤ₄₃, |S| = 7, r = 140, E_Q = 2^137;
  - 599: r = 300, 2^77;
  - 751, 971 and 669: about 2^73;
  - 719: 2^66.
- **No affordable exact level exists for them.**
- **At 863**, a ℤ₄₃-compressed solution is easy, and its fibre holds about 2^63.5 expected cores. But every lift split
  has a deficit of 36 bits or more, against about 26 needed.

**Theorems proved here** (some by an earlier OpenAI session; `report/openai_4p3_review.md`):
- Theorem α, Lemma R and coset-Parseval (obstructions);
- the compression normal form for k4 with odd |K|;
- **Theorem B:** one fixed border certificate completes every k4 core with odd |K|;
- Theorem T and Theorem C (no tensor-Paley or three-stratum product identity);
- the elliptic and cyclotomic seed families are closed.

**Also closed:** T-matrix routes at 179, 191 and 223 (`report/tmat_routes.md`), and the 599 lead and the 599 ×4 slice
(`report/b599_answer.md`, `report/b599x4_answer.md`).

**One open tension.** Kumar settled 643, 839 and 955, which our analysis had classed as out of reach rich hosts, by
non-exact or non-search means: heuristic pruning of compressed lists to the heaviest vectors, Galois-ring groups, and
the quaternary-conference route. **Exact-gate pessimism is not the last word on rich hosts.**

## 3. What to think about (suggestions; not limits)

1. **Reverse-engineer the August method.**
   - What procedure would produce b = 0 solutions with Φ ≈ 0.1 at p = 223–491 with |U| = 3–9?
   - What would solve unstructured length-166 and length-178 Legendre quadruples (664–712 free signs) with a hidden
     half-shift template?
   - What does that procedure's cost scale like? Would it reach length 598 (599), 718 (719), 738 or 750, or b = 0 at
     p = 599 with |U| = 13 (r = 47) or |U| = 23?
2. **Spectral-side methods.** Prescribe |f̂(χ)|² profiles, or phases, first, and recover the signs:
   - phase retrieval with rounding;
   - lattice methods;
   - convex relaxations with randomised rounding;
   - factor-graph or belief propagation on the Fourier constraints.

   What is the success probability per trial, and how does it scale?
3. **Heuristic last levels on rich hosts.** Accept incompleteness, as Kumar did successfully. Which rich host at an
   open order has the best expected yield per unit of work under a pruned or heavy-vector search? Calibrate the method
   on 643, 839 or 955.
4. **New arrays or groups:**
   - Theorem 7 with other k and S, or other target functions;
   - non-abelian groups and cocyclic matrices (de Launey–Horadam);
   - Ito-type or weaving constructions;
   - two-circulant cores with richer borders.

   Which open orders does each make newly reachable?
5. **Ingredients first.** Which single new ingredient would close the most open orders through EKLO's Table 4
   constructions? Candidates are a quaternary conference matrix of order 646, 940, 1042, 1108 or 1192; T-matrices of
   order 127 or other p; WT(71); TT(48); BS(45, 44); a symmetric conference matrix of order 146; a skew Hadamard
   matrix of order 356. How would you find it?
6. **Learned or guided search**, in the spirit of Williamson–Yacobi–Zinn-Justin. Where, specifically, would a learned
   proposal distribution beat exact staged joins, and how would you validate that on solved orders before running it
   on open ones?

## 4. What to deliver

Up to **three ranked proposals**. Each needs:
1. **The idea:** a precise statement of what is searched or constructed. Give arrays, groups, ansätze and algorithms
   explicitly enough to implement.
2. **Why it could work where §2 failed:** which assumption of the failed methods it escapes, or which structure of the
   successful ones it exploits.
3. **A cost model** with its assumptions. Give the predicted work per solution at one or more open target orders,
   compared with about 2^46 elementary operations on 4 cores.
4. **A calibration plan on solved orders:**
   - which known order or orders to reproduce blind (e.g. 892, 1244, 668, 1676, 1916, 2092, 3732, or smaller
     analogues);
   - the expected outcome;
   - pass and kill criteria fixed in advance.
5. **The first experiment for a Claude Code instance:** exact steps, inputs, expected run time, and what result
   would justify scaling up.
6. **Small runs you did yourself,** if any, kept under about 10 minutes each, with scripts in
   `theory/openai_fresh_*.py` and logs in `logs/openai_fresh/`.

Also include:
- **"Ideas considered and rejected,"** one line each, with the reason. This includes anything that turns out to be a
  rerun of §2.
- **A reconstruction of the August methods** (§3 item 1), if you reach one, with the evidence from the originals in
  `forensics/originals/`.

Mark claims [P] proved, [M] exact computation (script named) or [E] estimate.
**Never claim a Hadamard matrix** without an exact integer check of H Hᵀ = n·I by code sharing nothing with its
construction.

## 5. Working rules (sessions may time out: checkpoint early and often)

- **Setup:**
  - Run `git fetch origin codex/order-3452-compression claude/focused-planck-j1645f`.
  - Check out `codex/order-3452-compression` and merge `origin/claude/focused-planck-j1645f` into it.
  - If the merge is not clean, abort it and read the Claude branch with `git show origin/claude/focused-planck-j1645f:PATH`.
- **Pushing:** push only to `codex/order-3452-compression`. Never push to `main` or to `claude/focused-planck-j1645f`.
  Open no pull requests.
- **The answer file:** within your first ten minutes, create `report/openai_fresh_answer.md` with the skeleton of §4,
  and push it. Push again after each finished proposal, or every 30 minutes, whichever comes first. If the session
  ends early, the file says what is unfinished.
- **Credit:**
  - EKLO for the survey, Theorem 7 and the flatness analysis;
  - Alpöge, Voinov and Reynolds-Haertle for the twelve matrices;
  - Kumar (`https://github.com/vatskishan03/hadamard`, Apache-2.0) for his constructions and routes;
  - this repository for §2.
