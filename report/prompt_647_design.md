# Prompt: construct H(2588) (v = 647) in the 647 ×353 doubly-staged slice (self-contained, maximum reasoning)

*For a maximum-reasoning Claude Code session on `middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`.*

- **Tags:** [M] verified by exact computation in the repo (name the script), [P] proved (give the proof),
  [E] Gaussian estimate (name the estimator and its calibration), [L] literature recollection (check it; arXiv and
  OEIS are blocked from these containers).
- **Branch:** other sessions push to the same branch. Fetch and merge before every push. Never kill processes you did
  not start, and never `pkill -f` with a pattern that matches your own shell; kill by PID.

## 0. Goal and deliverables

Construct a Hadamard matrix of order **2588 = 4·647**. 647 is prime and the order is open.

This repo has built:
- H(2884) (v = 721, `cw/RESULT_2884.md`);
- H(2292) (v = 573, `report/RESULT_2292.md`);
- **H(2092) (v = 523, 2026-10-13)**, in a quartic-Paley slice, by a calibrated two-level stage-2 solver
  (`report/paley_523_answer.md`).

647 ×353 is the only live slice the last design session found (`report/next_design_answer.md`; read §0, §5, §9.3a,
§9.4 and §11 in full).

Deliverables, in order of value:
- **(A)** A verified H(2588): the core, the border, and H Hᵀ = 2588·I checked by an independent script that shares no
  code with the search (pattern: `scripts/verify_h2092.py`, `scripts/verify_h2292.py`). Also a publication bundle
  `pub/H2588/` in the format of `pub/H2092/`. The border and corner of `pub/H2292/border.json` work unchanged for any
  core with rows (2, 0, 0, 0) in block A, since the corner is −S·Tᵀ/2.
- **(B)** A calibrated cost per core for the cross level [M speeds × E counts, with the count calibrated exactly], and
  a production run if it is ≤ about 2^45 operations per expected core (≤ about 100 core-hours).
- **(C)** An exact obstruction: the cross level empty for every joint solution, or the joint stage 1 empty for every
  y_B orbit. That would prove the slice empty.

Write to `report/paley_647_answer.md`, with new scripts in `theory/p647_*.py` and logs in `logs/p647/`. Commit and push
results as they arrive, because containers restart. Commit trailers follow the repository's existing commits; no pull
requests.

## 1. The construction (exact)

**Bordered GS, b = 1.**
- Four ±1 sequences x_k on G = Z₆₄₆ (m = v − 1) with Σ_k PAF_k(g) = −4 for g ≠ 0.
- Then Σ_k r_k² = 4, so the rows are (2, 0, 0, 0) after reordering and negation. For every nontrivial character χ,
  Σ_k |x̂_k(χ)|² = 4v.
- Completion: the GS array of the circulants with back-identity R, border rows S ⊗ 1ᵀ and columns Tᵀ ⊗ 1
  (S, T ∈ H(4)), and the corner −S·Tᵀ/2 (`pub/H2092/verify.py` `build`).

**The slice ×353** [P] (`report/next_design_answer.md` §5.1):
- m = 646 = 2·17·19.
- 353 ≡ 13 (mod 17), of order 4: index 4 on Z₁₇*, f = 4.
- 353 ≡ 11 (mod 19), of order 3: index 6 on Z₁₉*, f = 3.
- Z₂ is fixed, and M = M₁₇ × M₁₉ is a product. Each x_k is constant on M-orbits; there are 70 per block.

**Character levels** on Z₂ × Z₁₇ × Z₁₉:

| level | characters | content |
|---|---|---|
| level 0 | χ = ψ = 1 | rows, alternating sums |
| stage 1A | χ = 1, ψ ≠ 1 | the Z₃₈ compression along Z₁₇: y_A on Z₃₈ |
| stage 1B | χ ≠ 1, ψ = 1 | the Z₃₄ compression along Z₁₉: y_B on Z₃₄ |
| cross | χ ≠ 1, ψ ≠ 1 | given (y_A, y_B), per block and Z₂-layer, a 5 × 7 ±1 label matrix Λ |

- In the label matrix Λ, the rows are {0, A₁..A₄} (Z₁₇ classes) and the columns {0, B₁..B₆} (Z₁₉ classes). Its
  weighted row sums are y_B and its weighted column sums y_A. So the 4 × 6 cell block has fixed margins.
- The remaining equations split by S/D on Z₆₄₆ = Z₂ × Z₃₂₃:

      Σ_k PAF₃₂₃(σ_k) = (645, −2, …),   Σ_k PAF₃₂₃(δ_k) = (647, 0, …)   at 20 lag classes,

  with σ = (layer₀ + layer₁)/2 and δ the half-difference. That is 12 + 12 new constraints, and
  ‖σ‖² + ‖δ‖² = 323 per block.

## 2. What is established [M unless marked]

- **Stage 1B** has exactly **495 orbits** under the slice's symmetry group: multipliers through Z₃₄*/⟨13⟩, per-block
  shifts by 323, negation and permutation of row-0 blocks. They come from 962 canonical α-quadruples.
- **Joint stage 1 for y_B #0: 20,928 joint solutions (y_A, y_B)**, from 928,080 S-solutions in 140 s. Orbits #1–#494
  are not run: the run died on #1, probably out of memory, since its S-join has more than 10⁶ matches. Up to the
  stabiliser of y_B (at most 6), that is about 3,500 inequivalent joint solutions.
- **Cross level, Gaussian [E, uncalibrated]:** H 79.4 bits, C 84.3 (C_S 49.8, C_D 49.4), Mahalanobis 0.48. So
  **log₂E ≈ −4.9 per joint solution**, about 2^7 inequivalent cores from y_B #0 if it holds.
- **Per-layer lists** for one joint solution: (1200, 584), (740, 1110), (1860, 1110), (710, 936). σ has about 2^17.7
  distinct values per block.
- The cost estimate is about 2^21 S-solutions expected per joint solution, then a cheap D-check, and about 2^40 per
  core [E]. **The S-join work is not measured**; its 2+2 bound is about 2^35.
- **No same-shape positive control exists.** 239 ×115 and 443 ×225 are exactly empty, as is 647 ×183 (all 4 y_B
  orbits fail the joint test).
- **Calibration history:**
  - whole-slice Gaussians of +46 to +102 sat on exactly empty slices about 15 times;
  - given one compression, the Gaussian of the other compression is meaningless, because it takes only a few discrete
    values per block (§9.4);
  - for the cyclic stage 2 of 523 and 573 the estimators were within 1–2 bits;
  - the cross level here has no further compression, so the discrete failure mode does not apply directly. But it
    is unproven.

## 3. Tools (read before writing new code)

| file | what it does |
|---|---|
| `theory/nd_cross.py` | cross level given (y_A, y_B): `Cross.lists` (exact per-layer labellings), `Cross.estimate` (lattice Gaussian), `solve_cross` / `planted_cross` (exact S-then-full join on orbit representatives, 35 per layer, to avoid the 5 GB arrays that killed a first attempt). **Written, not yet validated or run.** |
| `theory/nd_joint_exhaust.py` | symmetry group, canonical α-quadruples, restart-safe stage-1B enumeration and the joint test; `--plant N`, `JointTest.planted_real` |
| `theory/nd_ejoin.py` | exact 4-list energy-split join (numba open addressing, partitioned passes, growing buffer); `--selftest` 30/30 |
| `theory/nd_doubly.py` | `bordered_gs_cyclic` (completion); older, slow joint code |
| `theory/nd_s1_sd.py`, `theory/nd_lag0_screen.py` | exact compressed stage-1 solver; lag-0 screen |
| `theory/q523_wagner.py`, `q523_dcheck_nb.py`, `q523_prod.py` | the H(2092) kernels: radix-partitioned cache-resident join (`mjoin_rad2`), compiled D-check, restart-safe production with checkpoints and `--resume`. Read only; copy into `theory/p647_*` to adapt. |
| `pub/H2092/verify.py`, `scripts/verify_h2092_rebuild.py` | standalone verification patterns |

Loading a joint solution (from `report/next_design_answer.md` §11):

    X = Cross(646, 353, 17, 19)
    ys = enumerate_yB(...)
    rec, sols = X.J.test(ys[0])
    ya = X.J.ya_full(sols[j]['layers'])

## 4. Plan (reorder if you find better; every step's result goes in the report)

1. **Validate the exact cross solver.**
   - `planted_cross(X, yb, ya, seed)` on joint solutions of y_B #0 must recover every plant. Use at least 10 plants
     over different joint solutions.
   - Add a brute-force check on a reduced instance.
2. **Calibrate the S level.**
   - Run `solve_cross(..., s_only=True)` on several real joint solutions. Record the S-join time and work, and the
     exact S-solution count.
   - Compare the count with 2^(H_S − C_S) (H_S = Σ_k log₂ #distinct σ_k ≈ 4 × 17.7; C_S ≈ 49.8 from
     `Cross.estimate`). This is the only available calibration of the cross Gaussian; report the per-solution log₂
     errors.
   - Then measure the D-level success on the S-solutions found, as the 523 session did (`theory/q523_dfine_calib.py`,
     §3.4 of its report), against C_D.
3. **Cost per core** = (joint solutions per core) × (S-join + D-checks per joint solution), using the measured speeds
   and the calibrated counts.
   - If the S-join dominates, port the 523 join kernels (radix-partitioned passes, bitset prefilters).
   - If it is ≤ about 100 core-hours, go to step 4. If not, find the lever: middle levels, list reweighting, or
     choosing joint solutions with larger cross freedom. Cost it before running.
4. **Production.**
   - Deduplicate the 20,928 joint solutions of y_B #0 under the stabiliser of y_B (multipliers u acting only on Z₁₉,
     i.e. u ≡ ±(M mod 34) on Z₃₄; at most 6).
   - Run the exact cross solver over them on 4 workers, restart-safe, with per-joint-solution checkpoints and a tally
     script.
   - Any core: complete with `nd_doubly.bordered_gs_cyclic`, verify with a new independent `scripts/verify_h2588.py`,
     commit at once, then build `pub/H2588/`.
   - Cores arrive as a Poisson process. H(2092) came at about 3 expected cores, with probability about 0.05 of that
     wait. Do not abandon a calibrated run early.
5. **If y_B #0 is exhausted without a core,** resume the joint enumeration over orbits #1–#494
   (`nd_joint_exhaust.py 646 353 17 19 --dir logs/nd/joint --workers 1`). First add a cap on S-solutions or chunk
   the full join (the #1 run died out of memory).
6. **Positive control (strongly preferred before production; required before claiming (B)).**
   - Search small orders for a product slice of the same shape: G = Z₂ × Z_p × Z_q with M = M_p × M_q acting by
     cyclotomic subgroups on both primes and fixing Z₂, i.e. v = 2pq + 1, with v not necessarily prime.
   - Screen candidates exactly (lag 0, joint stage 1), and solve one end to end with the same pipeline to a verified
     Hadamard matrix.
   - If no such control exists below a stated bound, prove or exhibit that, and rely on the step-2 calibration.

## 5. Rules

- A construction counts only if it is verified exactly. Separate theorems, proofs and estimates.
- Calibrate every Gaussian you use against exact counts.
- **Compute:** 4 cores and about 15 GB of RAM. Size join buffers for 15 GB with headroom.
- Install `ortools` and `numba` with pip if they are missing.
- Write restart-safe runners: checkpoints, `--resume`, and the rule that a run adopts another run's checkpoint only
  when told to.
- Arm a watch or check-in so the session does not idle while workers run; an idle session was reclaimed once and
  killed the 523 workers.
- Report status in `report/paley_647_answer.md` §0 as the run proceeds: joint solutions done, S-solutions, D-checks,
  expected cores so far.
