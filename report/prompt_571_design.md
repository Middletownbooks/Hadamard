# Prompt: reach H(2284) (v = 571) inside the Paley + multiplier slices (self-contained, for a fresh Claude Code session)

*For a maximum-reasoning Claude Code session on `middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`.
Tags: [M] verified by exact computation in the repo, [E] Gaussian estimate, [L] literature recollection (check it;
paper search may be blocked).*

## 0. Goal and deliverables

Construct a Hadamard matrix of order 2284 = 4·571. 571 is prime; 570 = 2·3·5·19.

The route is the one that just produced **H(2292)** (v = 573, `report/RESULT_2292.md`): a bordered (b = 1)
Goethals–Seidel core in a **Paley(p) + multiplier slice**. For 571 the slices are rich in solutions but out of reach of
our staged pipeline (§3). I want a design that makes some 571 slice reachable at about 2^45 operations or less.

Deliverables, in order of value:
- **(A)** a verified H(2284). Use `paley_b1pipe.finish` or the independent `scripts/verify_h2292.py` pattern.
- **(B)** a method whose cost at 571 is ≤ 2^45 [E] and which reproduces a known core as a positive control: 419 or 573
  (both b = 1, p = 11), or the p = 5 control.
- **(C)** an exact obstruction showing a slice family is empty, so we stop.

Write to `report/paley_571_answer.md`, with scripts in `theory/p571_*.py`. Commit and push.
- Mark every claim [M] (naming the script), [P] (with a proof) or [E].
- Background jobs may be running; do not kill processes you did not start.

## 1. The architecture (exact)

**Bordered GS core (b = 1).** Four ±1 functions x_k on G = G′ × Z_p, of order m = v − 1 = n′p, with:
- Σ_k PAF_k = −4 off 0;
- row sums (±2, 0, 0, 0).

The border (S, T ∈ H(4), with the corner forced) completes it to H(4v); see `theory/paley_b1pipe.py` `finish`.

**Paley-coset ansatz.** On each coset g′ + Z_p:

    x_k = α·δ + u·χ + w·ε,

where χ is the Legendre symbol mod p, δ = [z = 0], ε = 1 − δ, and exactly one of u (twisted) or w (constant) is ±1.

**Multiplier slice.** The labels (α, type, sign) are constant on the orbits of g ↦ μg on cyclic G′ = Z_n′. This is the
orbit slice of the multiplier group (squares mod p) × ⟨μ⟩.

**Exact decomposition** (`report/paley_design_answer.md` §2, verified [M]). Write β = α − w, y = α + (p−1)w (the
Z_p-coset sums), and CC(f, g)(h) = Σ_x f(x) g(x + h), with all sums taken over the four blocks.
- **S1:** Σ_k PAF(y_k) = T_J, where T_J = 4m − 4(p−1) at 0 and −4p elsewhere (b = 1).
- **S2 (real):** Σ_k [PAF(β_k) + p·PAF(u_k)] = 4v·δ.
- **S2 (χ-part):**
  - p ≡ 3 (mod 4): Σ_k [CC(β_k, u_k) − CC(u_k, β_k)] = 0;
  - p ≡ 1 (mod 4): Σ_k [CC(β_k, u_k) + CC(u_k, β_k)] = 0. The symmetric form is new; it was validated by a verified
    H(164) at p = 5 (`data/core_v41_b1_paley_p5_symX_control.txt`).

**Constraints.**
- Norm identity: N_tw + (p−2)²·N₋ + p²·N₊ = T_J(0), where N₋ counts J′ cosets (w = −α), N₊ counts J cosets (w = +α)
  and K = N₋ + N₊.
- The w-screen (Lemma 4) and Lemma 6 (u row sums at the real characters) are in `theory/paley_wscreen.py`.

**Staged pipeline**, which found H(2292):
1. Designed gadgets (w on constant orbits), ranked by predicted stage-1 yield.
2. An exact 2+2 join for α on S1, which gives y.
3. My stage-2 engine `scripts/paley_s2join.py`: an exact 2+2 join for u given y.
   - Per-block u-lists: all signs on the twisted orbits.
   - **Z_d class split:** the Z_d-compressed S2 equations are additive over blocks, so only class pairs completable on
     the other side are joined. At 573 (d = 4) this cut each side from 2^34 to 2^23–2^25.
   - Prime-residue hash passes; signature dedup.

At 573 (p = 11, Z₅₂, ×29; slice log₂E = +23) stage 2 completed **1 of 141** y's, far better than the predicted 2^−12
to 2^−17.

## 2. Tools and data

| file | role |
|---|---|
| `theory/paley_slicepipe.py` (`guided`, `run_gadget`, `gadget_alpha_lists`) | gadget design + exact stage 1 |
| `scripts/run_slice_b1.py` | worker loop: stage 1 → stage 2 → finish |
| `scripts/paley_s2join.py` | exact stage-2 engine; `solve_y(A, W, v, p, n, mu, dcomp=d)`; p ≡ 1 supported |
| `theory/paley_slices.py`, `scripts/slice_sweep.py` | slice log₂E [E] |
| `theory/paley_split_slice.py v p n mu` | stage-1 vs stage-2 split per (K, N₋) [E]; patched for p ≡ 1 |
| `scripts/slice_sa.py` | joint annealing / tabu on orbit labels; **fails the 419 control** (§3) |
| `data/core_v573_b1.txt`, `data/core_v419_b1*.txt` | positive-control cores |

## 3. What is established for 571 (do not redo; challenge if wrong)

**Slice sweep** (b = 1; `logs/slice_sweep/`, `report/unattempted_2000_3000.md`, last section). Orbits are per block:

| p, G′ | μ (order) | orbits | log₂E [E] | constraint bits |
|---|---|---|---|---|
| 3, Z₁₉₀ | 27 (12) | 22 (sizes 1×2, 4×2, 6×6, 12×12) | **+40.6** | 200 |
| 3, Z₁₉₀ | 31 / 61 / 69 | 40 / 30 / 36 | +98 / +40 / +46 | 360 / 298 / 364 |
| 3, Z₁₉₀ | 9 = 99 | 18 | −10.8 | |
| 5, Z₁₁₄ | 11 (6) | 28 (1×2, 2×2, 3×12, 6×12) | **+20.6** | 292 |
| 5, Z₁₁₄ | 31 / 25 / 5 | 24 / 18 / 12 | −7.8 / −10.0 / −60.6 | |
| 19, Z₃₀ | 19 / 11 / 17 / 7 | 18 / 20 / 10 / 12 | −10.6 / −16.3 / −40.5 / −140 | |

**Stage split** (`logs/split_571_p3_mu27_partial_record.log`, `logs/split_571_p5_mu11.log`):
- **p = 3, ×27:** the norm identity forces **N₊ = 189** (N₋ free).
  - Stage 1 log₂E is +37 to +49 per (K, N₋).
  - Stage 2 has only 61–68 twisted u-bits against 107–116 constraint bits, so log₂E per y = **−39 to −53**.
  - Core log₂E is about −12 to +2 per K, for N₋ ≤ 14 (the run was stopped there).
- **p = 5, ×11:** stage 1 is +63 to +96; stage 2 is **−77 to −105 per y**; core is up to +10.8 (K = 104, N₋ = 43).
- **So the solutions exist [E] but are stage-2 starved.** Fixing y first is hopeless, since about 2^45 or more y's
  would be needed.

**Joint local search does not work** [M] (`scripts/slice_sa.py`):
- Orbit-label annealing and best-move tabu were run on the exact orbit-weighted objective; the incremental updates
  were unit-checked.
- Neither solves the 419 slice (best cost about 1,800–2,800), which the staged pipeline solves in minutes.
- Earlier, plain SA on the full reduced system failed at n′ ≥ 12 (`report/unattempted_2000_3000.md`).

**Other open orders** (for context): 647 (p = 17, Z₃₈, ×7, −4.7) gave 0 stage-1 solutions in about 2,200 gadgets;
599 at b = 5 is dead; 523, 537, 669, 719 and 739 at b = 1 fail the exact w-screen; 643 and 659 have no usable
multiplier slice.

## 4. Questions (reorder freely)

1. **Re-split the problem so that no stage is starved.** The constraint count is fixed, so the slack must be
   distributed differently. Candidates:
   - **Level decomposition of stage 2.** Z₁₉₀ = Z₂ × Z₅ × Z₁₉ (×27 acts as 1, ×2, ×8), and Z₁₁₄ = Z₂ × Z₃ × Z₁₉.
     S2 at the characters of each level involves only that level's compressions. Can u be solved level by level
     (characters trivial on Z₁₉ first, then the rest), with α also left partly free per level? The goal is to move
     stage-1 slack into stage 2.
   - **Mixed staging:** fix only part of y (for example its Z₁₉-compression, or the J/J′ gadget and the α-compressions),
     then solve the remaining α and u **jointly** by a join whose lists carry both α and u per block.
     - Per-block list size is the issue. For p = 3 ×27: α has about 22 − (#constant orbits) free orbits, u about 63/4
       bits per block. So per-block (α, u) lists of about 2^30, with class splitting.
     - Estimate the cost, and design the compressions/classes that keep it ≤ 2^45.
   - **Different variable split:** any exact reformulation (e.g. over Z[(1+√−3)/2] for p = 3, Eisenstein integers) in
     which the four-block condition becomes a norm-sum condition with fewer free variables per equation.
2. **p = 3 specifics.** With p = 3 the Paley twist is tiny:
   - a coset is (α, u, −u) twisted or (α, w, w) constant;
   - J cosets are forced to N₊ = 189.
   Is there an equivalent re-encoding of a p = 3 slice as a GS problem on G′ = Z₁₉₀ with ternary or quaternary
   alphabets (each coset is a symbol from {(α, u, −u), (α, w, w)}) that existing compression / PSD methods handle well?
   [L: Đoković–Kotsireas compression for GS/D-optimal searches.]
3. **The 19 factor.** Every 571 slice has a Z₁₉ factor in G′, with μ of order 3 or 6 on it.
   - Is there a second Paley-type twist (cubic or sextic characters of F₁₉, i.e. cyclotomic classes of order 3 or 6,
     rather than the quadratic character) that cuts the variables per block while keeping log₂E > 0?
   - Note that the plain doubly-quadratic-twisted slices were shown empty or negative by the previous design session
     (`report/paley_design_answer.md` §5).
4. **b ≥ 3 for 571.** v + 1 = 572 = 4·11·13. Options: c = 4 (dead in the Paley ansatz), c = 22, 26, 44, 52. Borders
   for c ≥ 10 are unknown. Is a border construction for c = 22 or 26 within reach, and do those Paley slices have
   positive log₂E with a staged-friendly split?
5. **Positive control.** Whatever route you propose must first reproduce a known core: 419 (×7 on Z₃₈) or 573 (×29 on
   Z₅₂, a core in `data/core_v573_b1.txt`), or a smaller p = 3 or p = 5 instance that you show is non-empty.

## 5. Rules

- A construction counts only if it is verified exactly (H Hᵀ = N·I).
- Separate theorems (with citation or proof) from estimates.
- For each route, give the exact objects, their cost in bits, the smallest checkable instance, and the first
  computation to run.
- Gaussian estimates have been wrong in both directions here:
  - they overcount at p ≥ 11 without a slice (rigid gadgets, Lemma 6);
  - they undercounted stage 2 at 573 by about 2^5–2^10.
  Calibrate any estimate you rely on.
