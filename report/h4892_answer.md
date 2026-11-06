# The 4892 slice decided exactly (v = 1223, b = 1), and the remaining cheap slices

Date: 2026-11-06. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_h4892.md` (the prompt of commit
2e41933). Scripts: `theory/h4892_*.py`. Logs: `logs/h4892/`. Data: `data/h4892_stage1.json`.

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
(being filled)
<!-- /NEIGH -->

## 4. Second part: the remaining cheap and tier-0 slices

<!-- PART2 -->
(being filled)
<!-- /PART2 -->

## 5. Next computation

<!-- NEXT -->
(being filled)
<!-- /NEXT -->
