# The 4892 slice decided exactly (v = 1223, b = 1), and the remaining cheap slices

Date: 2026-11-06. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_4892.md` (commit 2e41933). Scripts: `theory/h4892_*.py`. Logs: `logs/h4892/`. Data: `data/h4892_stage1.json`.

**Tags.** [P] proved; [M] exact computation (script named); [E] estimate. A cap or a timeout is never a zero.

---

## 0. Bottom line

<!-- BOTTOM -->
1. **4892: the slice `["b1", 1223, [], [47, 2], 13, 4, [[0, 23]]]` has no core [M].** No H(4892) comes out of it.
   - **Stage 1 is complete and confirmed** by two enumerators that share no code. Each finds **384 solutions up to
     block order and block sign**, which is **147,456 ordered quadruples**. The two sets agree as sets of explicit
     (y₀, y₁, y₂, y₃) (§1).
   - **Stage 2 is exactly empty for all 384 classes** (§2). Every stage-2 lift already fails at the ℤ₁₃
     compression, before the full 15-class join. Two independent codes show this:
     - the bucketed stage-2 join (`theory/h4892_stage2_run.py`): 0 matched S-buckets for every y;
     - a direct check on explicit length-13 compressed sequences (`theory/h4892_z13check.py`): 0 of 384.
   - Planted controls were recovered in every case: stage 1 4/4 (and method B's own plants, §1.3); the stage-2 join
     4/4 in production (6/6 in the trial run); the ℤ₁₃ check 6/6.
2. **The measured fibre was larger than the estimate.** Each y has block lists of sizes 9216, 20736, 20736 and 576.
   That makes the fibre 2^51.0 per y against the screen's 2^45.9 [E], so the plain 2+2 join would have cost
   2^36 in all. The ℤ₁₃ S-level bucketing made the real cost about 80 s for all 384 y.
3. **The 4892 neighbourhood** (§3) and **the remaining cheap and tier-0 slices** (§4): table in §4.
<!-- /BOTTOM -->

---

## 1. Step 1: stage 1, complete and confirmed

### 1.1 The slice

- G = G′ × ℤ₁₃ with G′ = ℤ₄₇ × ℤ₂, so |G| = 1222 = v − 1.
- S = U × M:
  - U = the squares in ℤ₄₇^* (order 23), acting on axis 0;
  - M = {1, 3, 9} ⊂ ℤ₁₃^*, with cosets C_j = 2^j M.
- The six U-orbits on G′ are {0} × {b}, QR × {b} and QNR × {b}, of sizes 1, 1, 23, 23, 23, 23.
- Letters y = α + 3σ: {±1, ±5, ±7, ±11, ±13}.
- Target: T₁ = 4840 at 0 and −52 at every other lag.
- Index order: (a, b) ↦ 2a + b.

### 1.2 Two enumerations [M]

**(a) Method B** of `theory/open_screen_stage1_exact.py --strong`, with every solution written out by
`theory/h4892_stage1_toolB.py` (`logs/h4892/stage1_toolB_full.json`):
- 718 pruned blocks with y[0] > 0;
- one row pattern, (0, 0, 0, 2);
- **384 solutions**, each unique up to block order and block sign (63 s).

**(b) An independent enumerator**, `theory/h4892_stage1.py`, which shares no code with the stage-1 tool, ab_lib or
ab_b0_stage1:
- The QR set and the orbits are recomputed from scratch.
- PAF_{G′}(y) is evaluated as the quadratic form Σ v_i v_j N_ij(h), with N counted by brute force over G′.
- All 10⁶ letter vectors are listed per block, keeping those with energy ≤ 4840 − 3·94 and r² ≤ 4. That leaves
  **2964 blocks** and 421 distinct PAF vectors.
- An ordered 2+2 join then gives **147,456 ordered solutions**.
- Every solution is re-checked on all 94 lags by the definition of PAF.

**Agreement** (`theory/h4892_stage1_compare.py`, `logs/h4892/stage1_compare.json`):
- The two canonical sets are equal (384 = 384).
- The ordered set of (b) is exactly the expansion of the canonical set by the 24 block orders and 16 block signs:
  384 · 384 = 147,456.
- Every class has four distinct blocks, and there are 28 distinct blocks up to sign.
- Every class has rows (0, 0, 0, ±2).

**Symmetries.**
- The ordered list in `data/h4892_stage1.json` mods out nothing.
- The canonical list mods out only block permutation and block negation. Both preserve the slice and map stage-2
  fibres bijectively onto each other: a permuted or negated core is again a core [P].
- Not used for reduction, though they also preserve the slice:
  - translation by (0, 1) on G′;
  - the multiplier by a non-square on ℤ₄₇, which swaps the QR and QNR orbits;
  - the multipliers of ℤ₁₃, which leave y fixed.
- Step 2 runs all 384 classes.

### 1.3 Planted controls

**Enumerator (b): 4/4** (`logs/h4892/stage1_own.out`).
- Each plant is one row-2 block and three row-0 blocks with total energy ≤ 4840.
- The list is then complete for the planted target, since every block of a planted solution has energy ≤ Tp(0) − 282
  and Σ r² = 4.
- In each case the planted ordered quadruple is among the returned solutions.

**Method B** with `--planted 4` (`logs/h4892/stage1_methodB.jsonl`): **4/4 recovered**, and 384 solutions again.
This completes the run that did not finish within 900 s on 2026-11-06.

`data/h4892_stage1.json` contains:
- the 384 canonical classes as 94-vectors;
- a 56-entry block table;
- the 147,456 ordered solutions as index quadruples;
- the planted controls.

---

## 2. Step 2: exact stage 2 for every stage-1 solution [M]

### 2.1 The lift

- On G = ℤ₄₇ × ℤ₂ × ℤ₁₃ (index 13(2a + b) + z), each block is constant on 30 S-orbits.
- Given y, the unknowns are the class signs c_j(O) with Σ_j c_j(O) = σ(O).
- PAF_G of a block is a quadratic form in its 30 orbit values. It is recorded on one lag per {orbit, −orbit} class:
  15 nonzero classes, each with target −4.
- The orbit and lag bookkeeping (`theory/h4892_stage2.py`) was checked against a direct PAF over G on random
  orbit-constant sequences.

### 2.2 The search (`theory/h4892_stage2_run.py`)

- **The ℤ₄ symmetry [P].** The multiplier z ↦ 2z of ℤ₁₃ fixes y, maps C_j to C_{j+1} and preserves the (constant)
  target. So the largest block (the pivot) is restricted to the least pattern of each orbit, which cuts it from
  20736 to 5184 patterns.
- **The S-level.** Each block list is bucketed by its S-level. This is the PAF of the compression along G′ (to ℤ₁₃)
  on its two nonzero lag classes, a linear image of the 15-vector (`theory/h4892_slevel.py`).
- **The join.** Pairs {a, b} | {c, d} are formed only in S-bucket pairs whose sums can meet. Inside those, a 64-bit
  linear hash of the 15-vector is joined by sort and searchsorted, and every hit is re-checked exactly. A real core
  would then be checked on all 1221 lags.

### 2.3 Result (`logs/h4892/stage2_Gp.jsonl`)

| quantity | value |
|---|---|
| stage-1 classes run | 384 (all) |
| block list sizes per y | {9216, 20736, 20736, 576} in every class (fibre 2^51.0 per y) |
| distinct S-level values per block | 22–511 |
| S-bucket pairs, left / right (max) | 75,721 / 9,788 |
| matched S-bucket pairs | **0 for every y** |
| left / right pairs visited | 0 / 0 |
| cores | **0** |
| time | 80 s in all |

**Planted** (same log, first 4 lines): **4/4 recovered**.
- Each plant uses a random stage-1 class and random class signs, with the target set to the planted 15-vector.
- The z ↦ 2z restriction is complete only for invariant targets, so plants with the pivot on draw the pivot pattern
  among its orbit minima, and the other plants switch the pivot off.
- Recovery means the exact planted quadruple is returned. For example, one plant visited 7.1·10⁶ left and
  1.3·10⁶ right pairs and returned 39 cores of its own target, the planted one among them.

### 2.4 The independent check: the ℤ₁₃ compression (`theory/h4892_z13check.py`, `logs/h4892/z13check.json`)

**Why it suffices [P].** If x₀, …, x₃ is a core, then the compressions w_k(z) = Σ_{g′} x_k(g′, z) satisfy
Σ_k PAF_{ℤ₁₃}(w_k) = −4 · 94 = −376 at every t ≠ 0.

**The check** shares no code with the stage-2 files:
- For each block it builds the exact set of reachable (w(C₀), …, w(C₃)): 79, 1027 or 1105 vectors.
- Each vector is expanded to a length-13 sequence, its 12 nonzero-lag PAF values are computed by the definition, and
  the four blocks are joined 2+2.
- **0 of 384 classes have a quadruple. Planted 6/6.**

So the slice is empty: no lift exists, for any stage-1 solution [M].

**Calibration.** The screen estimated, per y:
- C_S = 24.2 bits;
- H_S = 29.0 bits, so about 2^4.8 S-level solutions per y [E].

The real count is 0 for all 384 y. The estimate was positive and the exact answer is zero, as at every open order
decided before.

---

## 3. Step 4: the 4892 neighbourhood

<!-- NEIGH -->
**Same group, q = 13.** G = (ℤ₄₇ × ℤ₂) × ℤ₁₃ throughout. U ranges over the subgroups of ℤ₄₇^* of order 23 and 46
(order 2 gives 48 U-orbits; it fits 24 ≤ r ≤ 160 only at s = 2). s ranges over {2, 3, 4, 6, 12}.

**Inclusion [P].** If U ⊆ U′ and M ⊆ M′, every (U′, M′)-invariant core is (U, M)-invariant. So slice(U′, s′) is
contained in slice(U, s) when U′ ⊇ U and s′ | s with the same q (M′ ⊇ M).

**Stage-1 methods.**
- **A, B:** `theory/open_screen_stage1_exact.py --strong --planted 2`, with logs in `logs/h4892/part2/4892_*.jsonl`.
- **G:** `theory/h4892_stage1_gen.py d s`, the q = 13 generalisation of the step-1 enumerator, with logs in
  `logs/h4892/stage1_gen_d*_s*.json`.
- **Y:** `theory/h4892_stage1_any.py`, the generic independent enumerator, with logs in
  `logs/h4892/part2/any_4892.jsonl`.

The control (23, 4) is reproduced by G and by Y: 384 classes and 147,456 ordered quadruples.

| U order | s | r | stage 1 | stage 2 | verdict |
|---|---|---|---|---|---|
| 23 | 4 | 30 | 384 (A ✗ mem; B, own; sets equal) | 0 for all 384 (§2) | **empty [M]** |
| 23 | 2 | 18 | 0 (G, planted 2/2) | – | **empty [M]**; also ⊂ (23, 4) [P] |
| 46 | 2 | 12 | 0 (G, planted 2/2) | – | **empty [M]** |
| 46 | 4 | 20 | 0 (G, planted 2/2) | – | **empty [M]**; also ⊂ (23, 4) [P] |
| 23 | 3 | 24 | 0 (A, B, agree; planted 2/2 each, after the odd-s plant fix; G, Y, planted 2/2) | – | **empty [M]** |
| 46 | 3 | 16 | 0 (G, planted 2/2) | – | **empty [M]**; also ⊂ (23, 3) [P] |
| 46 | 6 | 28 | 0 (B, planted 2/2; G, Y, planted 2/2; A ✗ mem) | – | **empty [M]** |
| 23 | 6 | 42 | **58,048 classes** (Y, planted 2/2; G `--canonical`; the two sets are equal, `logs/h4892/stage1_s6_compare.json`; B timed out at 3600 s, and so did A) | **0 at the ℤ₁₃ compression [M]** (below) | **empty [M]** |
| 46 | 12 | 52 | same stage-1 problem as (46, 6) (identical letters and T₁): 0 | – | **empty [M]** |
| 23 | 12 | 78 | same stage-1 problem as (23, 6): 58,048 classes | f = 1, so α is not fixed by y; the energy bound kills none | **undecided**; f = 1 is outside the screen's menu |
| 2 | 2 | 144 | <!-- S2U2 -->queued<!-- /S2U2 --> | – | – |

**The ℤ₁₃ compression decides (23, 6) [M].** The real-y re-read [E] (`theory/h4892_reread.py`,
`logs/h4892/reread.jsonl`, 5 real y) gives:
- E per y ≈ −111, C₂ ≈ 192;
- S-level (L = G′) solutions per y ≈ 2^−37;
- cost ≈ 2^136 per core.

That is far outside the 2^40 rule, so no stage-2 search was run. Two cheap necessary tests at the ℤ₁₃ compression
decide it exactly anyway:
1. **The energy bound** (`theory/h4892_energybound.py` [P + M]).
   - At lag 0, Σ_k E(w_k) = Σ r_k² + 4(q − 1)|G′| = 4516.
   - Here E(w) = w(0)² + f Σ_j W_j². The value w(0) = Σ_O |O| α(O) and the sum Σ_j W_j = Σ_O |O| σ(O) are both fixed
     by y, and every W_j ≡ Σ_O |O| (mod 2).
   - Minimising Σ_j W_j² under these constraints gives a closed-form lower bound per block.
   - **47,232 of the 58,048 classes** exceed 4516.
2. **The exact reachable-set test** (`theory/h4892_zqnp.py`): numpy reachable sets of the compression, a range test
   per lag class, and the window filter.
   - The window filter keeps a block value only if it lies in [T − Σ_others max, T − Σ_others min], iterated to a
     fixed point [P].
   - It kills **all 10,816 remaining classes** (`logs/h4892/zqnp_4892_s6_survivors.json`).
   - The window filter never rejected a planted target (4/4).
   - On the s = 4 control it again gives 0 of 384, with planted 4/4 found by its join.

**q = 47** (G′ = ℤ₁₃ × ℤ₂, s = 23, f = 2):

| U on ℤ₁₃ | r | stage 1 | verdict |
|---|---|---|---|
| order 12 | 96 | 0 (B, planted 2/2; Y, planted 2/2; A ✗ mem) | **empty [M]** |
| order 6 | 144 | <!-- Q47U6 -->running (A, B); Y capped at 3·10⁶ partial blocks<!-- /Q47U6 --> | – |

**The other v = 1223 rows of the screen.**
- The four other b1 keys, with q = 47, s = 2 and U of order 4, 3, 2 or 1, all died at screen A.
- The k4 slices (G′ ∋ ℤ₆₁ or ℤ₅, q = 61 or 5) are A-dead or have C₂ ≥ 165 with negative E per y. The cheapest k4 row
  has cost 2^73.6 [E]. None was run.
- There is no b0 slice.
<!-- /NEIGH -->

## 4. Second part: the remaining cheap and tier-0 slices

<!-- PART2 -->
**Methods.**
- **A and B:** `theory/open_screen_stage1_exact.py KEY --strong --planted 2`, at 3600 s per key.
  - Keys where method A used the whole hour (A runs before B) were rerun with `--methods B` (logs `B_*.jsonl`).
  - Keys where B ran into the tool's 6 GB address-space cap were rerun one at a time with `--aslimit 13`
    (`M_*.jsonl`).
- **Second method:** `theory/h4892_stage1_any.py`, which shares no code with the tool.
  - It is a level-wise numpy enumeration under the energy and row bounds, with the k4 compression normal form where
    it applies (Q = ℤ₄, |K| odd).
  - The join is bucketed by pair energy and keyed on a 64-bit linear hash with exact recheck. Every solution is
    re-checked on all lags by the definition.
  - Plants are drawn by row class, with total energy ≤ T₁(0), so the list is complete for the planted target.
  - Before each run its set-up (orbits, T₁) was compared with the tool's: all equal
    (`logs/h4892/part2/any_setup_check.jsonl`).
- The table is rebuilt by `python3 theory/h4892_part2_table.py`.

**Two tool fixes** (prompt §4, with regression):
- `b_setup` now accepts `[axis, order, generator]` U entries, which fixes the 5892 key. The regression test
  `theory/h4892_test_bsetup.py` compares 1295 screen keys against the previous version: 1223 identical set-ups,
  0 differences, and the 72 previously failing keys are exactly the 3-entry ones. The 5892 orbits also agree with
  `ab_sweep.build`.
- `planted_y` seeded every orbit with ±1. For odd s, σ is odd, so ±1 is not a letter and the plants were not letter
  vectors. Both A and B then "missed" them on 4892 s = 3 (log kept as
  `4892_q13_s3_U23.invalid_plants.jsonl`). The plants now start from the least letters; for even s the random draws
  are unchanged.
- Every slice decided in `logs/open_screen/` has even s, so no earlier result is affected.

<!-- TABLE2 -->
<!-- /TABLE2 -->

**7572 (v = 1893, b1, ℤ₁₁ × ℤ₂² × ℤ₄₃, s = 6, U of order 10) is the one non-empty stage 1 among them.**
- **Stage 1** has **1360 classes**. B and the second method agree as sets of explicit quadruples
  (`logs/h4892/part2/7572_stage1_compare.json`; data `data/h7572_stage1.json`).
- **Re-read at 6 real y [E]:** E per y ≈ −1620 and mahal² ≈ 47, giving a cost of about 2^1650. This is far beyond
  the 2^40 rule.
- **It is decided exactly anyway, at the ℤ₄₃ compression:**
  - `theory/h4892_zqnp.py` finds that, for every one of the 1360 classes, the sum of the four blocks' least
    reachable compression energies exceeds the required Σ_k E(w_k) = 4 + 42 · 4 · 44 = 7396.
  - The closed-form bound of `theory/h4892_energybound.py` gives the same: every class has lower bound 7564 > 7396.
  - The cause is w(0) = Σ_O |O| α(O) ≈ ±42 in every block.
  - So **7572's slice is empty [P + M]**.
  - Small-modulus tests (w mod 3, …, 9; `theory/h4892_zqmod.py`) pass every class. Only the integer energy decides
    it.

**k4 slices with odd |K|** (3004 and 4156): a core would have completed by Theorem B. 3004 is empty by two methods.
<!-- /PART2 -->

## 5. Next computation

<!-- NEXT -->
1. **Gate every staged slice at the ℤ_q compression before costing it.**
   - Every non-empty stage 1 met here has died at the coarsest stage-2 level, the compression along G′ to ℤ_q:
     4892 s = 4 by the ℤ₁₃ join (and by mod 23), 4892 s = 6 by the energy bound plus the window filter, 7572 by the
     energy bound.
   - That level is cheap: reachable sets of ≤ 2·10⁶ vectors per block, and the closed-form bound
     `theory/h4892_energybound.py` needs only y.
   - The screen's Gaussian estimate saw these levels as costly but feasible: S-level solutions per y of 2^3.3 at
     4892 s = 4 and 2^−37 at s = 6. The real counts were 0.
   - **Proposal:** add to `theory/open_screen_stage1_exact.py` a stage "1½", run on every stage-1 solution found:
     `h4892_energybound` first, then `h4892_zqnp` (range test, window filter, join). Treat a slice as alive only if
     some y passes.
2. **A joint (y, w) enumeration for the slices whose stage 1 is too large for B** (5892, 5164, 7916, 4792, 4156,
   4748).
   - The energy identity Σ_k [w_k(0)² + f Σ_j W_{kj}²] = Σ r_k² + 4(q − 1)|G′| is a strong extra row condition, but
     it is useless as a per-block prune while y is dominated by ±1 letters (there the bound is ≈ r²).
   - So the remaining route is a memory-lean method B: a streaming pair join by bucket, written to disk, instead of
     materialising all pair arrays. Alternatively, CP-SAT (`--methods C --strong`) with the ℤ_q-level energy equation
     added as a redundant constraint.
   - Budget: up to 1 h per slice on 4 cores.
3. **The two 4892 neighbours left undecided** are (U₂, s = 2), with 48 U-orbits, and (q = 47, U₆).
   - Both have C₂ > 150 or negative E per y at the screen [E], so even a non-empty stage 1 would not be searched.
   - Running them is cheap only once item 2 exists.
4. **(U₂₃, s = 12) at 4892** (f = 1, outside the screen's menu) has the same 58,048-class stage 1 as s = 6.
   - Its lift needs the f = 1 form of the ℤ₁₃ test, where α is free in every orbit.
   - Carrying both decodings into the reachable set is easy. But with 12 cosets the reachable sets can reach about
     10⁸ vectors per block, so the energy window has to be applied while the sets are built, not after.
<!-- /NEXT -->
