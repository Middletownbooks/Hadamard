# The six undecided Paley-type c = 4 configurations: all empty

Date: 2026-10-17. Branch `claude/focused-planck-j1645f`. Answers the request "Run the six undecided paley type configs":
the c = 4 (b = 3) Paley-type configurations left undecided in `report/alpoge_design_answer.md` §7 (over the
6·10⁷ block-vector cap or out of memory): 839 q = 11 with ⟨8⟩; 955 q = 7, s = 2 and 3; 955 q = 17, s = 4;
991 q = 13, s = 4 with ⟨8⟩ and with ⟨4⟩.
Scripts: `theory/ab_s1_chars.py` (stage 1), `theory/ab_s2_zq_screen.py`, `theory/ab_s2_z2zq_screen.py`,
`theory/ab_s2_qzq_screen.py` (exact stage-2 screens), `theory/ab_s2_s4.py` (the stage-2 Gaussian, for comparison).
Logs: `logs/ab/paley6/`.

**Tags.** [M] exact computation (script named); [P] proved here; [E] Gaussian estimate.

---

## 0. Bottom line

| v | slice (Q = Z₂², q, s, H', U) | stage 1 [M] | stage 2 | status |
|---|---|---|---|---|
| 839 | q = 11, s = 2, H' = Z₁₉, U = ⟨8⟩ (order 6) | **0** (21,500 column-valid blocks, join empty) | – | **empty [M]** |
| 955 | q = 7, s = 2, H' = Z₂ × Z₁₇, U = Z₁₇^* | **0** (the ψ₁ system has no solution) | – | **empty [M]** |
| 955 | q = 7, s = 3, same H', U | **0** (the ψ₁ system has no solution) | – | **empty [M]** |
| 955 | q = 17, s = 4, H' = Z₂ × Z₇, U = Z₇^* | **0** (725,018 blocks over 402,487 ψ₀ quads, join empty) | – | **empty [M]** |
| 991 | q = 13, s = 4, H' = Z₁₉, U = ⟨8⟩ (order 6) | **0** (2,661,224 blocks over 20,760 ψ₀ quads, join empty) | – | **empty [M]** |
| 991 | q = 13, s = 4, H' = Z₁₉, U = ⟨4⟩ (order 9) | **119,858,496** | every y has no lift (exact compression screens) | **empty [M]** |

- **All six are empty [M]: none of these slices contains a core, so none gives an H(4v)** (the border and H(4v)
  steps of the pipeline therefore do not arise).
  - Five have no stage-1 solution at all.
  - **991 at s = 4 with ⟨4⟩** has 119,858,496 stage-1 solutions, all verified exactly (raw file sha256
    `cd199cde…61e1`, a 100,000-row uniform sample committed), and **none of them lifts to a core**. Exact
    necessary conditions for stage 2 — the core identity compressed along subgroups of Q × H' (§3) — kill every one:
    - the four Z₁₃-level conditions (compress along all of Q × H', one per Q-character): 72,243,520;
    - the three joint Z₂ × Z₁₃ conditions (pairs of characters): 44,856,320 more;
    - the joint Q × Z₁₃ condition (compress along H' only): the last 2,758,656.
- **Bonus: the two configurations the earlier report called "dead [E]" are empty too [M]** (§4):
  - 991 q = 13, s = 3, ⟨4⟩ (98,304 y's, earlier "−78 bits per y [E]"): every y fails the Z₁₃-level condition, under
    each of the four Q-characters separately;
  - 751 q = 17, s = 4, ⟨3⟩ (311,296 y's, earlier "−72 bits per y [E]"): the Z₁₇-level conditions kill 253,952, a
    joint Z₂ × Z₁₇ condition the last 57,344.
- With these, **all 61 Paley-type c = 4 configurations with ≤ 16 G'-orbits are decided empty [M]** (the sweep of
  `report/alpoge_design_answer.md` §7 had 53 empty, 2 "dead [E]", 6 undecided).
- **The stage-2 Gaussian is unreliable in this family** (§5): on a 991 ⟨4⟩ y it predicts 2^76.7 lifts at the Z₁₃
  level, where the exact count is 0 (a block's letters fix its Z₁₃ value at 0 so that its lag-0 energy alone
  exceeds the target). The earlier "dead [E]" labels rested on that Gaussian; the exact screens replace it here.
- Every planted control passes (all six configurations): each character system contains the planted solution and
  the full pipeline recovers the planted quadruple.

## 1. Stage 1 by the character decomposition over Q [P, M]

### 1.1 The identity [P]

Setting (`report/alpoge_design_answer.md` §1): G = Q × H' × Z_q, Q = Z₂², the slice group S = ⟨U on H', K on Z_q⟩
(K the index-s subgroup of Z_q^*), H = {0} × H' × Z_q, target T = 4m at 0, −4b on H \ 0, +4 off H (b = 3). Stage 1
lives on G' = Q × H': four S-invariant y_k with letters y = ε₀ + f(ε₁ + … + ε_s) on each G'-orbit and
Σ_k PAF_{G'}(y_k) = T1 (T summed over the Z_q fibres).

For a character ψ of Q put Z_k[ψ](h) = Σ_q ψ(q) y_k(q, h). Since ψ is real and ψ(q)ψ(q + δ) = ψ(δ),

  PAF_{H'}(Z[ψ])(η) = Σ_h Σ_{q,q'} ψ(q)ψ(q') y(q, h) y(q', h + η) = Σ_δ ψ(δ) PAF_{G'}(y)(δ, η),

so Σ_k PAF_{G'}(y_k) = T1 ⟺ for all four ψ: Σ_k PAF_{H'}(Z_k[ψ]) = T^ψ, T^ψ(η) = Σ_δ ψ(δ) T1(δ, η)
(the 4 × 4 Hadamard transform over Q is invertible). On each H'-orbit o the column of letters (y(q, o))_q and the
column of character values (Z[ψ](o))_ψ determine each other.

Weighted sums [P]: ẑ_k(ψ) = Σ_o |o| Z_k[ψ](o) is the Fourier coefficient of y_k at the G'-character (ψ, 1_{H'});
the core identity Σ_k |x̂_k(χ)|² = 4b for χ ≠ 1 trivial on H gives Σ_k ẑ_k(ψ)² = 12 for ψ ≠ 1, and ẑ_k(ψ) is even
(a ±-sum of 4|H'| odd letters), so ẑ_k(ψ) ∈ {0, ±2}. For ψ = 1, ẑ_k is the block's row sum.

### 1.2 Algorithm and completeness [P]

`theory/ab_s1_chars.py`:
1. per ψ: every vector Z[ψ] with energy ≤ T^ψ(0) and allowed weighted sum (ψ = 1: a row sum of an admissible row
   quadruple, taken ≥ 0, which fixes the block sign; ψ ≠ 1: 0, ±2); the 4-block system is solved on the distinct PAF
   signatures (all multiset quads; energy-split hashed join, light pairs tabulated in bounded bands, every quad
   checked exactly);
2. a vector is kept only if its signature occurs in a solution of its own system;
3. blocks: per kept ψ₀ vector, every triple of kept ψ₁, ψ₂, ψ₃ vectors whose (block, orbit) columns are transforms
   of letter columns;
4. per ψ₀ signature quad, a hashed 2+2 join of the four block lists on the ψ₁..ψ₃ PAFs (the ψ₀ part holds by
   construction); every hash match is re-checked exactly;
5. every result is checked again directly on G' (Σ_k PAF_{G'} = T1 at every lag class), canonicalised (block order;
   sign of zero-row-sum blocks) and written as raw int8 rows with a resume checkpoint.

Completeness: every stage-1 solution passes 1–4 (its ψ-vectors meet the bounds, form solutions of their systems,
its blocks are column-valid and its ψ₀ signatures form a quad), and the joins enumerate all candidates; hashes give
no false negatives (exact arithmetic mod 2^64). [P]

### 1.3 Validation [M]

- 839 q = 11, ⟨4⟩ (the earlier decided neighbour): 40,960 solutions, **set-equal** to the earlier independent run
  (`logs/ab/paley6/chars_ctl_839_u4.jsonl.gz`), every one verified.
- 991 q = 13, s = 3, ⟨4⟩ and 751 q = 17, s = 4, ⟨3⟩: 98,304 and 311,296 solutions, the counts of the earlier direct
  method (`logs/ab/paley_c4.jsonl`).
- The character systems: banded and single-band joins give identical quad sets on three configurations; the
  per-character counts for 839 ⟨8⟩ equal those of the first (vector-level) implementation.
- Planted controls (`--plant 1`, a random quadruple of realistic energy with its own target), all six
  configurations: every character system contains the planted quad, and the planted quadruple is recovered
  (`logs/ab/paley6/plant_*.log`; 955 q = 7, s = 3 joined only on the planted ψ₀ quad, whose 1,032,192 solutions
  include it).

## 2. Stage-1 results [M]

Per configuration: H'-orbit weights; per character the block vectors, distinct PAF signatures and signature quads
(ψ₁, ψ₂, ψ₃ coincide by symmetry); the vectors kept (ψ₀; ψ₁ = ψ₂ = ψ₃); the column-valid blocks; the stage-1
solutions. Logs `logs/ab/paley6/s1v2_*.log`, `s1v3_*.log`.

| configuration | weights | ψ₁..₃: vectors / sigs / quads | ψ₀: vectors / sigs / quads | kept | blocks | solutions |
|---|---|---|---|---|---|---|
| 839, q = 11, s = 2, ⟨8⟩ | 1, 6, 6, 6 | 4,047 / 2,024 / 939 | 10,166 / 10,166 / 21,776 | 3,467; 1,384 | 21,500 | **0** |
| 955, q = 7, s = 2 | 1, 16, 1, 16 | 1,639 / 478 / **0** | – | – | – | **0** |
| 955, q = 7, s = 3 | 1, 16, 1, 16 | 1,639 / 478 / **0** | – | – | – | **0** |
| 955, q = 17, s = 4 | 1, 6, 1, 6 | 14,793 / 3,611 / 2,437 | 122,746 / 53,198 / 402,487 | 13,268; 1,783 | 725,018 | **0** |
| 991, q = 13, s = 4, ⟨8⟩ | 1, 6, 6, 6 | 5,613 / 2,807 / 1,357 | 9,792 / 9,792 / 20,760 | 3,152; 1,798 | 2,661,224 | **0** |
| 991, q = 13, s = 4, ⟨4⟩ | 1, 9, 9 | 211 / 46 / 370 | 371 / 138 / 2,726 | 357; 211 | 27,202 | **119,858,496** |
| 991, q = 13, s = 3, ⟨4⟩ (§4) | 1, 9, 9 | 211 / 46 / 370 | 209 / 79 / 1,258 | 191; 211 | 3,704 | 98,304 |
| 751, q = 17, s = 4, ⟨3⟩ (§4) | 1, 5, 5 | 383 / 69 / 1,113 | 2,902 / 865 / 119,415 | 2,699; 383 | 22,868 | 311,296 |

The two q = 7 configurations of 955 die in the first character system: no four ψ₁-vectors (weighted sums 0 or ±2)
have PAFs summing to T^ψ₁ = (3708, −112, −112, −112). Runtimes (one core): 2 s for the q = 7 pair, 7 s for 839,
10.5 min for 955 q = 17, 58 min for 991 ⟨8⟩, about 1 h 45 min for 991 ⟨4⟩ (dominated by writing and checking
1.2·10⁸ solutions).

## 3. Stage 2 of 991 ⟨4⟩: exact compression screens [P, M]

### 3.1 The conditions [P]

For a character ψ of Q put x̄_k^ψ(t) = Σ_{q,h} ψ(q) x_k(q, h, t), a K-invariant vector on Z_q. Then

  Σ_k PAF_{Z_q}(x̄_k^ψ)(τ) = Σ_{δ,η} ψ(δ) Σ_k PAF_G(x_k)(δ, η, τ) = Tq^ψ(τ) := Σ_{δ,η} ψ(δ) T(δ, η, τ),

so each core satisfies four Z_q-level conditions. For a stage-1 y, the value of x̄_k^ψ on a K-orbit of Z_q is a sum
over the G'-orbits j of |j| ψ(q_j) times the sign the chosen fibre pattern puts there; the patterns of the letter
y_kj range over its fibre independently for each j.

### 3.2 Exact count per y [M]

`theory/ab_s2_zq_screen.py`: the distribution of x̄_k^ψ over block k's fibre depends only on the multiset of
per-orbit (signed) contribution sets, the block's ψ-type. Per type an exact DP (exact integer multiplicities), then
the PAF signatures (lag 0 and every class but the last, which the fixed ψ-sum ẑ_k(ψ) of y determines); per distinct
quadruple of types the exact number of lifts by a 2+2 convolution pruned by the lag-0 budget. Zero is exact (no
cancellation); a y with count 0 for some ψ has no lift at all.

### 3.3 Validation [M]

- On the v = 159 analog (q = 13, s = 4, H' = Z₃ with −1; 928 stage-1 y's) the screen's counts equal a direct
  per-y computation (`theory/ab_s2_s4.zq_exact`), 2^64.34 and 2^64.36.
- Planted targets (random lifts x, T' := Σ_k PAF_G(x_k)): the screen never kills y = compress(x) (on the 991 ⟨4⟩
  slice through `theory/ab_s2_zq_screen.run`: 3 logged trials, `validation_planted_zq_991.log`, after 10 unlogged
  ones during development; on v = 159, 6 trials with an independent dict convolution,
  `validation_planted_zq_159_dict.log`).
- The example that started it: for a y whose fourth block has x̄(0) = 68 forced by its letters, the lag-0 energy of
  that block alone is at least 4,672 > Tq(0) = 3,964 (the DP's minimum equals the bound).

### 3.4 Results at the Z₁₃ level [M]

Over all 119,858,496 stage-1 solutions of 991 ⟨4⟩ (811 block types; 1,056 ψ₀ and 36,248 ψ₁, ψ₂, ψ₃ type
quadruples; `logs/ab/paley6/zq/zq_991_q13_s4_u4_full.log`):

| condition | y's with no lift |
|---|---|
| ψ₀ (trivial character) | 66,689,344 |
| ψ₁, ψ₂, ψ₃ (each) | 22,587,200 |
| at least one of the four | **72,243,520** |

47,614,976 y's pass all four.

### 3.5 The joint Z₂ × Z_q levels [P, M]

For an element r of order 2 of Q, compress along ⟨r⟩ × H': x̄_k(a, t) = Σ_{q ∈ a + ⟨r⟩, h} x_k(q, h, t) on
Z₂ × Z_q (Q/⟨r⟩ ≅ Z₂). Summing the core identity over the same subgroup gives [P]

  Σ_k PAF_{Z₂ × Z_q}(x̄_k) = T_r,  T_r(α, τ) = Σ of T over the points of G that project to (α, τ).

This is the pair of characters {ψ₀, ψ_r} of §3.1 taken jointly (x̄(0, ·) ± x̄(1, ·) are the two compressions), so it
is stronger than both. Per block the two halves have independent exact DP distributions (each over the six
G'-orbits of its two Q-cosets); the joint states are the pairs within the lag-0 budget, and their PAF signatures on
Z₂ × Z₁₃ (all lag classes but one) are joined 2+2 with a 64-bit linear hash and energy bands
(`theory/ab_s2_z2zq_screen.py`). A zero count is exact: a true solution always matches its hash. Validation [M]:
planted targets (v = 55 and v = 71, the H' = 1 analogs with the same Z₁₃ and Z₁₇ fibres) are never killed; on
planted cases the counts equal an independent pure-Python exact integer count, for all three r (800 …
197,623,584 lifts; `logs/ab/paley6/zq/validation_*.log`).

Results on the 47,614,976 Z₁₃-survivors of 991 ⟨4⟩, the three r in sequence (each on the previous survivors;
`logs/ab/paley6/zq/z2r*_991_q13_s4_u4.log`):

| r | y's screened | block types | type quadruples | no lift | survive |
|---|---|---|---|---|---|
| 1 | 47,614,976 | 618 | 26,612 | 44,348,416 | 3,266,560 |
| 2 | 3,266,560 | 126 | 2,440 | 507,904 | 2,758,656 |
| 3 | 2,758,656 | 84 | 1,640 | 0 | 2,758,656 |

### 3.6 The joint Q × Z_q level [P, M]

Compressing along H' only, x̄_k(q, t) = Σ_h x_k(q, h, t) on Q × Z_q, gives [P]

  Σ_k PAF_{Q × Z_q}(x̄_k) = T_Q,  T_Q(δ, τ) = Σ of T over the points of G that project to (δ, τ),

the four characters jointly (it contains every level above). Per block the four cosets have independent exact DP
distributions; the joint states are the 4-tuples within the lag-0 budget; their PAF signatures on Q × Z₁₃ (12 lag
classes, one implied) are aggregated by (lag-0 energy, 64-bit linear hash) and joined 2+2 in energy bands, with the
2+2 split that minimises the pair lists (`theory/ab_s2_qzq_screen.py`). Aggregating by energy and hash together
means a hash collision can only add to a count, never move a true solution out of its band, so **a zero count is
exact**. Validation [M]: on planted targets of the H' = 1 analogs v = 55 and v = 71 the counts equal an independent
pointwise exact count (every lift enumerated, PAF on Q × Z_q computed directly; counts 1 … 49,184), and no planted y
is killed (`logs/ab/paley6/zq/validation_qzq_crosscheck_*.log`).

Result [M] on the 2,758,656 survivors of §3.5 (287 block types, 89,152 type quadruples;
`logs/ab/paley6/zq/qzq2_991_q13_s4_u4.log`, per-quadruple counts in `qzq2_991_q13_s4_u4.jsonl.gz`): **every type
quadruple has count 0, so no y lifts to a core.
991 q = 13, s = 4 with ⟨4⟩ is empty.** (The first pass, with hash-only aggregation, gave the same; the reported run
is the one with energy-and-hash aggregation.)

## 4. The two earlier "dead [E]" configurations [M]

Both stage-1 sets were regenerated with the character solver (the earlier counts, 98,304 and 311,296, reproduced)
and screened exactly.
- **991, q = 13, s = 3, ⟨4⟩ (earlier: −78 bits per y [E], "dead").** Each of the four Z₁₃-level conditions kills
  all 98,304 y's (one ψ₀ type quadruple and 24 per nontrivial character, every count 0): **empty [M]**
  (`logs/ab/paley6/zq/zq_ctl_991_q13_s3_u4.log`).
- **751, q = 17, s = 4, ⟨3⟩ (earlier: −72 bits per y [E], "dead").** The Z₁₇-level conditions kill 253,952 of the
  311,296 y's (all by ψ₀; the nontrivial characters kill none). On the remaining 57,344 the joint Z₂ × Z₁₇ level
  with r = 1 has count 0 for all 128 type quadruples: **empty [M]** (`logs/ab/paley6/zq/z2r1_ctl_751_q17_s4_u3.log`).

## 5. The Gaussian against exact counts [E vs M]

The stage-2 Gaussian (`theory/ab_s2_s4.py`: H2 exact, the lattice Gaussian of Σ_k PAF at the lag classes under
uniform lifts) was compared with exact counts at the Z₁₃ level, the one level where both are available per y.
- v = 159 analog (q = 13, s = 4, H' = Z₃ with −1): exact − Gaussian = +0.93 … +1.05 bits on 4 y's (Mahalanobis
  distance 4.2–4.3). In this regime it is calibrated to about a bit.
- 991 ⟨4⟩: on y #5 of the stage-1 file the Gaussian predicts 2^76.7 lifts at the Z₁₃ level (Mahalanobis distance
  2.8), the exact count is **0**. The letters of the fourth block fix its value at 0 ∈ Z₁₃ to x̄(0) = 68, so its lag-0
  energy alone is at least 4,672 > 3,964: a hard support bound that no Gaussian sees.
- 991 ⟨4⟩, full stage 2: the Gaussian gave −37.9 … −43.1 bits per y on six sampled y's (rank 23, Mahalanobis
  distance 6.1–6.5); five of the six have exact count 0 already at the Z₁₃ level, the sixth at the Z₂ × Z₁₃ level
  (r = 2).

So in this family the stage-2 Gaussian is not a usable estimator: its earlier "dead" verdicts (§4) happen to be
right, but for a reason it cannot see. No conclusion here rests on it.

## 6. Reproduce

```
python3 theory/ab_s1_chars.py logs/ab/paley6/cfg_<name>.json logs/ab/paley6/s1v3_<name>.ybin            # stage 1
python3 theory/ab_s1_chars.py logs/ab/paley6/cfg_<name>.json /tmp/p.ybin --plant 1 [--plant-quad-only]    # control
python3 theory/ab_s2_zq_screen.py logs/ab/paley6/cfg_<name>.json logs/ab/paley6/s1v3_<name>.ybin \
        logs/ab/paley6/zq/zq_<name>.jsonl --keep-alive logs/ab/paley6/zq/alive_<name>.ybin                 # Z_q screens
python3 theory/ab_s2_z2zq_screen.py logs/ab/paley6/cfg_<name>.json IN.ybin OUT.jsonl R [--shard i/n] \
        --keep-alive ALIVE.ybin                                                     # Z_2 x Z_q screen, r = R in 1, 2, 3
python3 theory/ab_s2_qzq_screen.py logs/ab/paley6/cfg_<name>.json IN.ybin OUT.jsonl [--shard i/n] \
        --keep-alive ALIVE.ybin                                                     # Q x Z_q screen
python3 logs/ab/paley6/zq/planted_screen_run.py logs/ab/paley6/cfg_991_q13_s4_u4.json 3 21                 # Z_q test
python3 logs/ab/paley6/zq/z2zq_crosscheck_planted.py logs/ab/paley6/calib/cfg_a55_q13_s4.json 4 9         # Z_2 x Z_q
python3 logs/ab/paley6/zq/qzq_crosscheck_planted.py logs/ab/paley6/calib/cfg_a55_q13_s4.json 4 5          # Q x Z_q
```
Runtimes (4 cores): stage 1 of the five empty configurations 2 s – 58 min each; 991 ⟨4⟩ stage 1 about 1 h 45 min
on one core (5.75 GB of raw solutions, regenerable; sha256 in `logs/ab/paley6/s1v3_991_q13_s4_u4.ybin.sha256`); the
Z₁₃ screens 40 min; the Z₂ × Z₁₃ screens about 45 min; the Q × Z₁₃ screen 15 min.

## 7. Files

- Code: `theory/ab_s1_chars.py` (stage 1), `theory/ab_s2_zq_screen.py`, `theory/ab_s2_z2zq_screen.py`,
  `theory/ab_s2_qzq_screen.py` (exact compression screens), `theory/ab_s2_hzq.py` (H' × Z_q level, validated by
  brute force on v = 159, not needed in the end), `theory/ab_s2_s4.py` (stage-2 Gaussian and exact Z_q counts),
  `theory/ab_full_chars.py` (unstaged solver, too heavy for v = 159).
- Logs: `logs/ab/paley6/` (configurations, stage-1 logs, planted controls), `logs/ab/paley6/zq/` (screens,
  validation), `logs/ab/paley6/calib/` (small analogs).
