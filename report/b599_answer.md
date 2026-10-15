# H(2396), v = 599: the b = 5, c = 6 Paley(11) lead and the other routes

Date: 2026-10-15. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_599_design.md`.
Scripts: `theory/b599_*.py`. Logs: `logs/b599/`.

**Tags.** [M] exact computation (script named); [P] proved here (proof in the text); [E] Gaussian estimate (estimator
and calibration named); [L] literature recollection, unchecked.

---

## 0. Bottom line

1. **(A) No H(2396) was built.** No route found in this session reaches a core within the budget. §8 ranks the routes.

2. **The lead slice has essentially no stage-1 solutions; the "+48 bits" premise was an artefact [E, calibrated;
   M for the exact searches].**
   - Stage 1 of the lead (G′ = Z₂ × Z₃ × Z₉, −1 on Z₉, letters ±1, ±9, ±11) must have a Q-compression z that solves
     the Q-level exactly: Σ_k |ẑ_k(ψ)|² = 20 at the five nontrivial characters of Q = Z₆ and Σ_k z_k² = 416. Those z
     form 9 classes (420 compression quadruples up to symmetry, §2.1) [M].
   - Conditioned on each z, the y-quadruples are counted exactly and the 15 remaining lag classes are estimated by the
     lattice Gaussian (`theory/b599_flevels.py`). Summed over all 420 quadruples: **2^4.9 ordered stage-1 solutions**.
     The two classes with a known border hold **2^0.3** (347 template class, 56 quadruples) and **2^3.3**
     (monomial/Singer class, 44 quadruples). Divided by the symmetries that fix z (translations of constant blocks,
     the multipliers of Z₉, block signs), that is ≪ 1 distinct solution [E].
   - The whole-slice estimate (+48.1) sampled y with the stage-1 energy only. There Σ_k z_k² has mean 3,588 and
     standard deviation 1,188 (minimum 1,008 in 2,000 samples), against the target 416, and the Q-character powers run
     1,757–6,253 against the target 20 (§2.3). The Gaussian evaluated that far lower tail of a positive quadratic form.
   - **Calibration** (§3): on the same-shape rung v = 215 (Paley(7), G′ = Z₂ × Z₃ × Z₅) the conditional Gaussian
     agrees with exact counts within 0.6 bit where solutions are plentiful (2^17 per z). Where they are scarce
     (Gaussian 2^6.3–6.6 per z) it **undercounts by 1.3–5.5 bits**, because 93% of those solutions have every block
     fixed by a nontrivial symmetry (a twisted multiplier or affine map) that a Gaussian of generic blocks cannot see.
     Same-shape rungs with Paley(11), the lead's letters, are exactly empty: v = 203 and **v = 335 (all 296
     compression quadruples) [M]**.
   - So the lead was searched exactly for symmetric solutions (§4, `theory/b599_symblock.py`, validated on v = 215):
     **every block fixed by an order-3 element: 0 solutions over all 420 quadruples [M]**; order-3 or the 15-dim
     involutions: **0** over the bordered quadruples run so far (§4).
   - With about 0.22 cores per y [E], even a handful of stray y's would not be worth the stage-2 run. **The lead is
     dead as a route [E + M].**

3. **The slices next to the lead split into "cheap stage 2, no stage 1" and "rich stage 1, expensive stage 2" [E].**
   - **H′ = Z₃² with −1 on both axes ("neg23")**: the Z₃² twin of the lead. Stage 1 2^4.9 (same 9 classes), stage 2
     −1.1 to −1.9 bits per y, D_def 21, about 2^43 per core. Same verdict as the lead, but a far larger symmetry
     group for symmetric searches (GL(2,3), §5).
   - **H′ = Z₃² with −1 on one axis ("Z₃²-a", 36 G′-orbits)**: stage 1 is **rich**, 2^45.7 over all classes and
     2^42.4 / 2^42.8 in the two bordered classes. An S/D stage-1 generator costs about 2^21 per y with complete
     σ-lists, about 2^32 with lists of 2^24 [E]. But stage 2 is −9 to −13 bits per y with D_def 28.4–30.5 and C_S
     77–80 (six typical y's): **2^57.9–2^61.1 per core** by the cost model calibrated at 347 and 523, which is
     2^13–2^16 over budget.

4. **(C) Exact obstructions found here [M]** (§7):
   - 598 = 2·13·23, b = 1, the doubly-staged multiplier slices ×3, ×49, ×31, ×25, ×45, ×51: **stage 1B is empty in
     every one** (`theory/nd_joint_exhaust.py`, plants pass). With the Paley(23) ×27 slice already dead (w-screen),
     the b = 1 multiplier slices of 598 with ≤ 72 orbits are all dead or lattice-obstructed (§7).
   - The lead restricted to order-3-symmetric blocks: empty (all 420 compression quadruples).
   - Same-shape c = 6 slices with Paley(11) at v = 203 and v = 335, and v = 131 (Paley(7), n = 3): stage 1 empty.

5. **(D)** The ranked list is §8. Nothing is within budget. The two concrete next computations are (i) a three-level
   stage-2 solver for Z₃²-a, whose level structure could cut the plain 2^58 toward 2^45, and (ii) the remaining
   symmetric families of neg23 and of the lead (the large involution families), which are exact and cheap per
   element but need a leaner join.

---

## 1. Setting and tools

- Lead (prompt §2): v = 599, b = 5, c = 6, h = 99, G = Z₂ × Z₃ × Z₉ × Z₁₁, H = Z₉ × Z₁₁, Q = Z₂ × Z₃,
  Paley(11) on Z₁₁ (squares, f = 5), U = −1 on Z₉. Stage 1 on G′ = Z₂ × Z₃ × Z₉: 30 G′-orbits per block, letters
  α + 5(p + n) ∈ {±1, ±9, ±11}; T1 = 2176 at 0, −220 on H′ ∖ 0, +44 off H′ (`theory/b599_lib.py`, self-test
  against `ab_sweep.build`) [M].
- **Same-shape family** `theory/b599_fam.py`: c = 6, G′ = Z₂ × Z₃ × Z_n (−1 on Z_n), Paley(q): v = 6nq + 5.
  599 = (9, 11); the ladder (3, 11) → 203, (5, 11) → 335, (7, 11) → 467, (3, 7) → 131, (5, 7) → 215,
  (9, 7) → 383, (7, 3) → 131′, (9, 3) → 167. Self-tested against `ab_sweep.build` for every member [M].
- **Q-level classes** `theory/b599_qlevel.py` and, for any slice with U trivial on Q, `theory/b599_cond.py` [M].
- **Conditional stage-1 estimator** `theory/b599_flevels.py` (family), `theory/b599_cond.py` (generic): exact DP
  counts of the y-quadruples with given z and stage-1 energy, uniform sampling, lattice Gaussian
  (`theory/q523_gauss.gauss_bits`) on the remaining lag classes; S-level and exact D-freedom h_D for the S/D budget
  [E, counts exact].
- **Exact counters**: `theory/b599_exact.py` (block lists, naive 2+2 hash join), `theory/b599_calib.py` (block
  lists, the 647 kernel `theory/p647_join.join4` imported unchanged, optional subsampling with an unbiased
  hits/p² estimate). They agree on 5/5 planted targets with non-zero counts (288, 288, 144, 648, 288) and with
  `theory/ab_s1sd.py` at v = 203 [M].
- **Symmetric-block searches** `theory/b599_symblock.py` (family), `theory/b599_symgen.py` (generic, Aut(H′) =
  units of Z_n or GL(2,3)) [M]; validation in §4.

---

## 2. Stage 1 of the lead

### 2.1 The Q-level classes [P, M]

Every coset sum z_k(q) = Σ_{h ∈ H′} y_k(q, h) is odd, so z_k = c_k·1 + 2u_k with c_k odd and u_k ≥ 0 integral,
min u_k = 0. The Q-level (`report/alpoge_design_answer.md` §1.2) is Σ_k |ẑ_k(ψ)|² = 20 for ψ ≠ 1, i.e.
Σ_k |û_k(ψ)|² = 5 [P]. Parseval gives Σ_k ‖u_k − ū_k‖² = 25/6, so each u_k has range ≤ 2 and u_k ∈ {0, 1, 2}⁶.
`theory/b599_qlevel.py` lists every such quadruple up to per-block translation and negation, block order and t → −t:

| class | shapes (u_k as indicator sets / spikes) | profiles | row solutions | Σ 2^log₂E [E] |
|---|---|---|---|---|
| 0 | 0, 0, δ, 2δ (**monomial / Singer border**) | (1,1,1), (4,4,4) | 44 | 2^3.28 |
| 1 | 0, 0, δ, (0,0,1,0,2,1) | | 44 | 2^1.93 |
| 2 | 0, 0, {two in a layer}, (0,0,2,1,1,1) | | 44 | 2^1.27 |
| 3 | 0, 0, {two in a layer}, (0,1,1,1,0,2) | | 44 | 2^1.37 |
| 4 | 0, δ, {two in a layer}, (0,1,2,1,1,1) | | 44 | 2^1.50 |
| 5 | 0, δ, Z₂-line, (0,1,1,2,1,1) | | 44 | 2^1.60 |
| 6 | 0, δ, (0,0,1,0,1,0), (0,1,1,1,1,2) | | 44 | 2^1.62 |
| 7 | 0, {two in a layer}, Z₂-line, {one + two complementary} (**347 template border**) | | 56 | 2^0.35 |
| 8 | 0, {two in a layer}, (0,0,1,0,1,0), (0,0,1,0,1,1) | | 56 | 2^0.29 |
| all | | | 420 | **2^4.89** |

Rows: Σ_k (6c_k + 2Σu_k)² = 2396. The 347 template class reproduces the template transfer of
`report/alpoge_design_answer.md` §3 exactly: 56 row solutions, each fitting the frame (`theory/ab_borderfit.fits`),
with 8 fitting translation variants each [M].
Log: `logs/b599/flevels_599_all.jsonl`, `logs/b599/levels_bordered347.jsonl`, `logs/b599/levels_monomial.jsonl`.

### 2.2 Per quadruple [E]

For each z the exact count H1 of y-quadruples with these compressions and stage-1 energy is 84.0–94.7 bits; the
remaining 15 lag classes (20 classes minus energy and the four Q-level combinations) cost C = 94.3–96.1 bits
(rank 15, Mahalanobis ≤ 0.31). So log₂E(stage 1 | z) runs from +0.17 (best monomial quadruple) to −12.0. The S/D
split of stage 1 gives h_D 9.4–13.8 and D_def 25.5–30.1: each quadruple has 2^16.2–2^26.0 S-solutions and needs
2^25.5–2^30.1 of them per stage-1 solution.

### 2.3 Why the whole-slice estimate said +48 [M]

`ab_gauss.stage1` samples y uniformly with the stage-1 energy (and rows) and fits one Gaussian to all 19 lag classes.
The Q-level is four linear combinations of them. Under that sampler (`ab_gauss.stage1_sampler`, 2,000 samples):
Σ_k Σ_q z_k(q)² has mean 3,588 and standard deviation 1,188 (minimum seen 1,008); the target is 416. The summed
power at a nontrivial Q-character runs 1,757–6,253 on average (min/max over ψ); the target is 20. The target sits far
in the lower tail of a positive quadratic form, where a Gaussian is far too generous. Conditioning on z (§2.2) removes
that tail entirely.

---

## 3. Calibration of the conditional estimator [M against E]

(`theory/b599_exact.py`, `theory/b599_calib.py`, `logs/b599/exact_*.jsonl`, `logs/b599/calib/`)

| rung (n, q) | v | quadruples | conditional Gaussian (sum) | exact |
|---|---|---|---|---|
| (3, 11) | 203 | 180 | none in lattice / energy | **0** (also `ab_s1sd.py`: 0 α-quadruples) |
| (3, 7) | 131 | 104 | 2^4.84 (blocks of 3–9 vectors: degenerate) | **0** |
| (5, 11) | 335 | 296 | 2^−2.3 | **0** |
| (5, 7) | 215 | 3 rich quadruples | 2^17.09, 2^16.99, 2^17.09 | 2^16.46, 2^16.81, 2^17.49 (± 0.2–0.3, p = 1/64 subsample) |
| (5, 7) | 215 | 20 scarce quadruples | 2^6.27–2^6.66 | 2^7.91–2^11.87 (exact, p = 1) |

- **Rich regime:** within 0.6 bit.
- **Scarce regime:** the Gaussian undercounts by 1.3–5.5 bits. The solutions are structured: in one quadruple 3,744
  solutions use only 9 and 6 distinct blocks in positions 2 and 3. Γ = (affine maps of Q) × (Z₅*/±1) acts on single
  blocks and maps stage-1 solutions to stage-1 solutions. **10,560 of the 11,328 solutions of six scarce quadruples
  have every block fixed by a nontrivial element of Γ** (most often ×2 on Z₅ composed with t → −t + s). A sampler of
  generic blocks does not see this family.
- The low-pattern ansatz (each block uses ≤ 3 distinct coset patterns, `theory/b599_twopat.py`) finds none of these
  solutions: the structure is symmetry, not low rank.

Consequence for 599: the conditional Gaussian is trustworthy for generic solutions, and the scarce regime has to be
searched exactly for symmetric families (§4).

---

## 4. Exact searches of the lead's symmetric families [M]

`theory/b599_symblock.py`: Γ = (Q-affine maps) × (Z₉*/±1), 36 elements. Per block, the union over prime-order
elements whose Q-part fixes z_k of their fixed blocks (compression z_k, energy ≤ its cap); then `join4`; every hit
verified by the full G′ PAF. Validation: on v = 215 it returns exactly the 10,560 all-symmetric solutions of §3
(768 + 768 + 768 + 768 + 3,744 + 3,744), with no code shared with the stabiliser analysis.

| element types (order: fixed dimension) | list sizes per block | quadruples | solutions |
|---|---|---|---|
| order 3 (3:10 twisted translations, 3:18 multipliers ×2, ×4) | 64–24,112 | all 420 | **0** |
| order 3 + involutions of dimension 15 (translation by Z₂, e-shifting reflections) | up to 983,402 | bordered (classes 0, 7) so far | **0** (see log) |
| involutions of dimension 20 (reflections t → s − t) | 2^25–2^26 | – | not run (lists too large for a 4-list join) |

Logs: `logs/b599/sym/sym3_599_p*.jsonl`, `logs/b599/sym/sym15_599_p*.jsonl`.

---

## 5. The other c = 6 configurations of 599 [E]

Stage 1 by the conditional estimator (`theory/b599_cond.py`), stage 2 by the typical-y budget
(`theory/ab_typical_levels.py`). Cost per core by the model of `report/paley_523_answer.md` §3
(`report/alpoge_design_answer.md` §11.1): 2^(D_def + C_S − H_S/2) with H_S = H₂ − h_D; it gives 2^43.5 for the lead,
matching the prompt's 2^43.

| slice | G′-orbits | stage 1 (all classes; bordered) | per y | D_def | cost per core |
|---|---|---|---|---|---|
| lead: H′ = Z₉, −1 | 30 | 2^4.9; 2^0.3 + 2^3.3 | −2.2 | 21.5 | 2^43.5 |
| neg23: H′ = Z₃², −1 on both | 30 | 2^4.9 | −1.1 / −1.9 | 20.8 / 21.6 | 2^43 |
| **Z₃²-a: H′ = Z₃², −1 on one axis** | 36 | **2^45.7; 2^42.4 + 2^42.8** | −8.8 … −12.8 (six y) | 28.4 … 30.5 | **2^57.9–2^61.1** |
| H′ = Z₉ or Z₃², −1 on Q's Z₃ | 36 | (U acts on Q: not run) | −13.1 / −13.8 | 31.2 / 31.9 | 2^61–2^62 |
| H′ = Z₉ or Z₃², U = 1 | 54 | – | −9 | 44.5 | 2^84 |
| index 5 on Z₁₁ (f = 2), −1 on Z₉ / U = 1 / −1 on Z₃ and Z₉ | 30–54 | – | −83 / +32 / −213 | 123 / 156 / 135 | ≫ 2^100 |

**Z₃²-a in detail.**
- Stage 1, best bordered quadruples (class 0 and 7): H1 134.2 / 134.5, C 96.7 / 97.0, log₂E +37.5 per quadruple.
- S/D budget of stage 1 (`b599_cond.budget_sd`): C_S 56.0 / 56.3, h_D 17.9 / 16.7, 2^60.3 / 2^61.6 S-solutions per
  quadruple, D_def 22.8 / 24.0. With complete σ-lists (2^29 per block, too large to hold) an S/D generator costs
  2^(D_def + C_S − H_S/2) ≈ 2^21 per y; with σ-lists subsampled to 2^24 per block each S-solution costs
  2^(C_S − H_S/2 + 2·5) ≈ 2^8 and the D-check about 2^9, so about 2^32 per y [E].
- Stage 2, six typical y: H₂ 124–126, C₂ 134.4–136.8, C_S 77.3–79.6, h_D 26.5–29.0, 2^17.3–2^21.7 S-solutions per y
  (`logs/b599/budgets_z3z3a.jsonl`). It needs about 2^11 y per core, and the S-join dominates:
  2^(C_S − H_S/2) ≈ 2^29 per S-solution and 2^28.4–2^30.5 S-solutions per core.
- The single middle level of `ab_stage2.py` does not change the per-S-solution cost. With a second, nested level M,
  the cost per S-solution is 2^(C_S − H_S + H_M/2) + 2^((C_S − C_M) − (H_S − H_M)/2). That helps exactly when M is
  constraint-dense, C_M > H_M/2. The Z₃² factor gives several U-stable candidate levels; measuring (H_M, C_M) on real
  y's is the next computation (§8).

---

## 6. c = 8 configurations of 599 [E]

h = 74 = 2·37, H′ = Z₂, Q = Z₂³ (359 template class), staging prime 37. s = 2, 3, 4 were already dead (screen A or
energy). The others had died of MemoryError in screen A; their typical-y budgets (`logs/b599/budgets_c8.jsonl`):
s = 6: −111.9 per y (D_def 106.8); s = 9: −155.1 (D_def 163.8); s = 12: +194.3 per y but D_def 52.3 with h_D 135.7
(about 2^120 per core). s = 18: the estimator ran out of memory. **All dead [E].**

---

## 7. b = 1 at m = 598 = 2·13·23 [M]

`theory/nd_joint_exhaust.py 598 μ 13 23` (stage 1B along Z₂₃ to Z₂₆, then the joint test), plants first.

| μ | |M| (on Z₁₃, Z₂₃) | orbits | canonical α-quadruples | stage-1B solutions | plants |
|---|---|---|---|---|---|
| 3 | 33 (3, 11) | 30 | 219 | **0** | 4/4 |
| 49 | 66 (6, 11) | 18 | 2 | **0** | 3/3 |
| 31 | 44 (4, 11) | 24 | 2 | **0** | 3/3 |
| 25 | 22 (2, 11) | 42 | 18 | **0** | 3/3 |
| 45 | 12 (12, 2) | 70 | 28 | **0** | 3/3 |
| 51 | 22 (2, 11), −1 ∈ M | 40 | 18 | **0** | 3/3 |

All six are **exactly empty**. Every multiplier subgroup that acts on Z₂₃ with index 2 has letters ±1, ±21, ±23 on the
Z₂₃ fibres, and these die at stage 1B as Paley(23) ×27 dies by the w-screen. The remaining subgroups of Z₅₉₈* in
`logs/msweep_record/m598.jsonl` (×53, ×61, ×5, ×17, ⟨5, 31⟩, ⟨25, 51⟩) have the target outside the lattice of
their whole-slice Gaussian (a parity obstruction, not proved exactly here). ⟨45, 47⟩ and ⟨183, 231⟩ (48 and 72
orbits) have naive estimates −62 and −95. Logs: `logs/b599/joint/`.

---

## 8. Ranked list of the 599 routes (D)

| rank | route | status | next computation |
|---|---|---|---|
| 1 | c = 6, Z₃²-a (H′ = Z₃², −1 on one axis), Paley(11) | stage 1 rich and bordered [E, rich regime calibrated]; stage 2 2^58–2^61 per core [E] | a three-level stage-2 solver: measure (H_M, C_M) for the U-stable subgroups of Z₃ × Z₃² on real y's (an S/D stage-1 generator first, about 2^21 per y) |
| 2 | c = 6, neg23 and the lead | stage 2 2^43 per core; generic stage 1 ≪ 1 [E]; lead order-3-symmetric: empty [M] | the GL(2,3)-symmetric families of neg23 and the dimension-20 involution families of the lead (exact, need a join that streams 2^25-lists) |
| 3 | c = 6, 36-orbit slices with U on Q | stage 2 2^61–2^62 per core [E] | extend `b599_cond.py` to U acting on Q; only useful with rank 1 solved |
| 4 | c = 12 (h = 49, F₄₉), c = 10 (h = 59) | no border known (Singer ansatz exhausted at c = 10) | the corner-flatness search of `report/alpoge_design_answer.md` §3.1 |
| – | c = 8 (s = 6, 9, 12) | dead [E] (−112, −155 per y; D_def 52 for s = 12) | – |
| – | c = 6, index 5 or 10 on Z₁₁ | dead [E] (D_def 123–156) | – |
| – | b = 1, 598 multiplier slices | dead [M] (six stage-1B empty), the rest lattice-obstructed or −62/−95 | – |
| – | b = 0 over Z₅₉₉, T(599) | closed (`report/tmat_routes.md` §10) | – |

---

## 9. Reproduce

```
python3 theory/b599_lib.py                                  # lead self-test against ab_sweep.build
python3 theory/b599_fam.py                                  # family self-tests
python3 theory/b599_qlevel.py --v 599                       # the 9 Q-level classes and their row solutions
python3 theory/b599_flevels.py --n 9 --q 11 --out FILE      # conditional stage-1 estimates, all 420 quadruples
python3 theory/b599_exact.py --n 5 --q 11 --out FILE        # exact per-z stage 1 (v = 335: 0)
python3 theory/b599_calib.py --n 5 --q 7 --top 1 --p 0.015625     # rich v = 215 calibration
python3 logs/b599/calib/run_low.py 0 1                      # scarce v = 215 calibration (exact)
python3 theory/b599_symblock.py --n 9 --q 11 --types 3:10,3:18 --out FILE            # order-3 symmetric: 0
python3 theory/b599_symblock.py --n 9 --q 11 --types 3:10,3:18,2:15 --classes 0,7 --out FILE
python3 theory/b599_cond.py '<cfg json>' --out FILE         # conditional stage 1 for any slice with U trivial on Q
python3 theory/ab_typical_levels.py '<cfg json>' --seed S   # stage-2 budgets
python3 theory/nd_joint_exhaust.py 598 3 13 23 --dir logs/b599/joint --workers 2      # 598 x3: stage 1B empty
```
