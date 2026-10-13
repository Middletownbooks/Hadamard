# Prompt: design a tractable route to a new Hadamard order other than 2092 (self-contained, maximum reasoning)

*For a maximum-reasoning Claude Code session on `middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`.*

- **Tags:** [M] verified by exact computation in the repo, [P] proved, [E] Gaussian estimate, [L] literature recollection
  (check it; arXiv and OEIS are blocked from these containers).
- **Branch:** another session runs the 523 production on this same branch. Fetch and merge before every push, and
  never touch its files (`theory/q523_*`, `logs/q523/`) or its processes.

## 0. Goal and deliverables

Construct a Hadamard matrix H(4v) for some open v other than 523.
- Open v < 750: 537, 571, 599, 643, 647, 659, 669, 719, 739.
- Open v between 750 and 1000 (user-supplied list, not yet adjusted for recent Alpoge/Schneider work): 751, 789, 823,
  839, 859, 863, 883, 907, 917, 919, 933, 947, 955, 971, 991.
- We built H(2292) (v = 573) and H(2884) (v = 721) in this repo, and validated a pipeline on H(940) (v = 235).

No candidate currently has a costed, validated path; that is why this is a design session. Deliverables, in order of
value:
- **(A)** A verified H(4v): H Hᵀ = 4v·I exactly, checked by an independent script (pattern: `scripts/verify_h2292.py`,
  `scripts/verify_h2092.py`).
- **(B)** For one specific slice, a method costed per solution [E calibrated] at ≤ about 2^45 elementary operations
  (≤ about 100 core-hours), validated end to end on a smaller order of the same shape, together with the exact first
  computation to start it.
- **(C)** Exact obstructions that remove candidates (proofs or INFEASIBLE certificates).
- **(D)** A ranked list of the remaining candidates with exact status and the next computation for each.

Write to `report/next_design_answer.md`, with scripts in `theory/nd_*.py` and logs in `logs/nd/`. Commit and push.
Mark every claim [M] (naming the script), [P] (with proof) or [E].

## 1. Architecture (exact)

**Bordered GS, b = 1.**
- Four ±1 functions x_k on G of order m = v − 1 with Σ_k PAF_k(g) = −4 for g ≠ 0.
- Then Σ_k r_k² = 4, so the rows are (±2, 0, 0, 0). Every nontrivial character χ satisfies Σ_k |x̂_k(χ)|² = 4v.
- Completion: the GS array, border rows and columns S ⊗ 1ᵀ and Tᵀ ⊗ 1 (S, T ∈ H(4)), and the corner −SΣᵀT/4
  (`theory/paley_b1pipe.py`, `scripts/verify_h2292.py`).

**Slices.**
- A multiplier group M = ⟨μ⟩ ≤ Z_m*; each x_k is constant on M-orbits. `logs/msweep_record/m*.jsonl` (v < 750) and
  `logs/msweep34/m*.jsonl` (750–1000) list every cyclic M with its orbit count and a slice Gaussian log₂E [E].
- **Generalized Paley (GP) slice:** M fixes a factor Z_F pointwise and acts on a prime factor Z_q by the index-s
  subgroup; f = (q − 1)/s. Per Z_q-coset, x = α at z = 0 and c_j on the class C_j.
- **Exact two-stage split:**
  - stage 1 is the coset sums y = α + f·Σc_j on Z_F, with Σ_k PAF(y_k) = 4m − 4(q − 1) at 0 and −4q elsewhere;
  - stage 2 is the cyclotomic equations R(h, t) = −4 (`scripts/cyclo_stage2_cpsat.py`, docstring). Equivalently
    Σ_k Z_k Z_k^† = 4v·δ in O_K[Z_F], where K is the degree-s Gauss-period field (`theory/q523_lib.py`, class
    `Cyclo(q, s)`, general in q and s).

## 2. What is established

**The 523 pipeline** (`report/paley_523_answer.md`, read §0, §3 and §5):
- slice Z₁₈ × quartic Z₂₉;
- stage 1 by CP-SAT or by a residue-aware 4-list join (`theory/q523_ymitm.py`);
- stage 2 by the Z₂ split of Z₁₈ into S and D levels on Z₉, with an S-level two-level Wagner tree (Z₃ compression
  middle, then the order-9 equations) and a compiled D-check;
- about 2^44 per core, 75 core-hours, calibrated on H(940) at v = 235.

Its cost floor is the D-deficit D_def = C(D|S) − h_D ≈ 26 bits (pinned by the S-level energy equation) times the
S-step per S-solution. Read §3.3–3.6 and §9 for why it cannot go lower.

**Exact screens, all in the repo** (`report/slice_screens_20261013.md`; read it in full):
1. lag 0 + rows + Z₂ compression for GP slices (`theory/gp_z2screen.py`);
2. the α-congruence Σ_k PAF(β_k) ≡ −4q (mod 4f), with β = α + f·[s odd] (`theory/gp_alpha_screen.py`). It is the
   general form of the 523 "mod-8 XOR" condition and explains why naive stage-1 joins fail;
3. the exact compressed stage 1 of **any** slice along any Hall factor Z_K: fibre-orbit letters, M-invariant on
   Z_{m/K}, solved by CP-SAT (`theory/slice_compress_cpsat.py`).
   - **Composite slices are intersections [P]:** compressing along each prime factor gives a GP-type stage 1 with
     extra invariance.
4. earlier screens:
   - the w-screen, Theorem 5 and Lemma 6 for Paley slices (`theory/paley_wscreen.py`);
   - the θ_q criterion (`scripts/theta_screen.py`);
   - h0 and qnorm for c = 4 (`scripts/paley_bge3.py`).

**Exactly empty** [P/M]:
- 955 ×289 (rows);
- 971 ×81 (α reduction, plus an independent CP-SAT INFEASIBLE);
- 917 ×17 (rows + Z₂);
- 907 ×85 and ×91, 823 ×73, 789 ×29, 971 ×221 and ×101 (lag 0);
- 739 ×127 (lag 0);
- 537 ×9 sextic (stage 1 exhaustive), ×59 and ×15 (compression);
- every Paley slice of 523, 599, 659 and 739 (w-screen).

**Calibration warning:**
- slices at Gaussian +46 to +67 were exactly empty seven times;
- the Gaussian undercounted 573's stage 2 by 2^5–2^10 (there it worked).

So screen exactly first. An estimate is evidence only after the exact screens pass, and only when calibrated.

**Live candidates [E unless marked]:**

| v | slice | structure | log₂E | exact status |
|---|---|---|---|---|
| 739 (H(2956)) | ×215 | fixes Z₂; −1 on Z₉; index 8 on Z₄₁ (f = 5); 90 orbits per block | +61.9 | survives all screens; compressed stage 1 on Z₁₈ (10 orbit letters per block; letters α + 10j, j = −4..4) UNKNOWN after 900 s CP-SAT |
| 647 (H(2588)) | ×353 | fixes Z₂; index 4 on Z₁₇ (f = 4) × index 6 on Z₁₉ (f = 3), diagonal order 12; 70 orbits | +36.2 | compressions along Z₁₇ (to Z₃₈) and Z₁₉ (to Z₃₄) both UNKNOWN after 1200 s |
| 947 (H(3788)) | ×49 | **the H(2292) shape**: Paley(11) on Z₁₁ × order-3 multiplier on Z₈₆ = Z₂ × Z₄₃ (index 14, f = 3); 30 orbits on Z₈₆ | +37.3 (msweep); Paley split [E] (`theory/paley_split_slice.py 947 11 86 49`, `logs/nd/split_947_mu49.log`): stage 1 +33 … +48, **stage 2 per y −19 … −23 with 107–109 u-bits**, core +14 … +24.4 (K = 29–34) (573 had stage 2 −9.6 with 73 u-bits and core +23, and stage 2 was undercounted there by 2^5–2^10) | w-screen alive (`theory/paley_wscreen.py one 947 1 11 86`); compression along Z₄₃ **FEASIBLE** (`logs/gp_screen/c947_49_K43_y.json`); along Z₁₁ (the Paley stage 1 on Z₈₆) UNKNOWN at 1200 s |
| 971 | ×469 | −1 on Z₅; index 8 on Z₉₇ | +36.2 | compressed to Z₁₀ (6 orbits per block) UNKNOWN at 600 s |
| 947 | ×21 | composite (Z₁₁ index 5, Z₄₃ index 6) | +28.3 | compression along Z₄₃ UNKNOWN at 600 s |
| 859 | ×289 | composite (Z₁₁, Z₁₃) | +27.4 | both compressions UNKNOWN at 600 s |
| 647 | ×183 | composite | +24.6 | compression along Z₁₇ UNKNOWN at 600 s |

Also exactly empty: 947 ×45 (s = 3 odd forces |y| ≥ 13, so lag 0 fails [P]) and 971 ×79 (compressed
stage 1 INFEASIBLE [M]).

- 739 is below 3000 and so the most valuable target. The CP-SAT UNKNOWNs are not evidence of emptiness: CP-SAT also
  returned UNKNOWN on 523's stage 1 at 120 s / 2 workers and solved it at 323 s / 4 workers.
- Unexplored: noncyclic abelian G where m has a square factor (750 = 2·3·5³, 882 = 2·3²·7², 918 = 2·3³·17,
  990 = 2·3²·5·11, 788 = 4·197). The sweeps covered cyclic Z_m only. Also b ≥ 3 borders, which
  `report/unattempted_2000_3000.md` closed for v < 750 below 300 bits.

## 3. Questions (reorder freely; first record exact status for every live slice)

1. **Decide stage 1 exactly for 739 ×215.**
   - The slice is GP(Z₁₈, q = 41, s = 8) restricted to labels symmetric under Z₉ negation [P].
   - Use the α-congruence (here mod 4f = 20, with α symmetric on Z₁₈), the per-block row, Z₂ and Z₃ conditions, and
     a residue-aware symmetric 4-list join in the manner of `theory/q523_ymitm.py`. Or prove it empty.
2. **Stage 2 for 739 ×215.**
   - Port the 523 level budget (`theory/q523_levels.py`, `q523_gauss.py`) to s = 8 with the symmetry restriction.
     K is the octic field of conductor 41; −1 ∉ H since (q − 1)/s = 5 is odd.
   - Compute C₂, the S/D budget on Z₉ (the Z₃ middle level still exists) and D_def. Then cost the 523 route here.
   - Does the symmetry, which halves the variables and pairs the equations, raise or lower D_def? If the route is
     ≤ 2^45 per core, give the run plan and a positive control (a smaller order of the same shape, e.g. Z₁₈ or Z₆ ×
     index 8 on Z₁₇, v = 307 or 103, with the same Z₉ or Z₃ symmetry).
3. **647 ×353 as a doubly-staged slice.**
   - Each compression is exact. Both further compress to the same Z₂ data.
   - Design a three-level split: a joint stage 1 of the Z₃₈ and Z₃₄ compressions, consistent through Z₂; then the
     cyclotomic stage 2 on each factor; then the cross term on Z₁₇ × Z₁₉.
   - Cost it, or find an exact obstruction. The 17-side has quartic classes with f = 4 (letters ±1, ±7, ±9, ±15, ±17)
     and the 19-side sextic classes with f = 3.
4. **947 ×49, the H(2292) shape at a larger cofactor.**
   - The tools that built H(2292) apply as they stand: `theory/paley_slicepipe.py guided1 947 11 86 49 K …` for
     stage 1, `scripts/paley_s2join.py` for stage 2 (read `report/RESULT_2292.md` and `report/paley_design_answer.md`
     §7.5).
   - The obstacle is stage 2. There are about 108 u-bits per y against 73 at 573, and the class splits are weaker:
     Z₈₆ = Z₂ × Z₄₃ gives only Z₂, where Z₅₂ gave Z₄. So the plain join is about 2^54 per y.
   - Design a hierarchical stage 2 (the S/D split on Z₂, as in `theory/p571_hsolve.py` and the 523 D-check; or a
     split by the cubic classes of Z₄₃ under the order-3 multiplier), cost it per core, and calibrate on 573 and 419.
   - Recall that the Gaussian undercounted 573's stage 2 by 2^5–2^10.
5. **New slice families.**
   - Run the exact screens (§2) over every cyclic slice with log₂E > 0 that is not yet screened. The msweep files
     list them.
   - Then consider:
     - noncyclic G for the square-factor m above;
     - pairs (G, M) where M fixes a larger factor;
     - b > 1 borders for 750–1000.
   - Rank by exact status first and by a calibrated cost second.
6. **Is there a better route than slices?** The user is open to other constructions (Williamson-type or propus
   arrays, Ito's or Yamada's families, Miyamoto- or Goethals–Seidel variants). Mark recollections [L] and check
   applicability exactly for the open list before recommending any.

## 4. Tools (read before writing new code)

| file | what it does |
|---|---|
| `theory/gp_z2screen.py`, `theory/gp_alpha_screen.py`, `theory/gp_screen_candidates.py` | the exact GP screens |
| `theory/slice_compress_cpsat.py m mu K` | exact compressed stage 1 for any slice; controls: `40 9 5` (OPTIMAL), `536 59 67` (INFEASIBLE) |
| `theory/q971_stage1.py` | the α-reduction search with a planted self-test (a template for small F) |
| `scripts/cyclo_stage1_cpsat.py`, `scripts/cyclo_stage1_exact.py`, `scripts/cyclo_stage2_cpsat.py`, `scripts/cyclo_split.py` | GP stage 1 and 2 by CP-SAT, the exact reduction, the stage split [E] |
| `theory/q523_*.py` (read only; copy to `theory/nd_*`) | 523 pipeline: `q523_lib` (`Cyclo`), `q523_levels`, `q523_gauss`, `q523_ymitm`, `q523_sdsolve`, `q523_wagner*`, `q523_dcheck*`, `q523_prod`, `q523_finish`, `q523_planted` |
| `scripts/paley_s2join.py` | the stage-2 join that built H(2292) (Z_d class splits, prime-residue passes) |
| `theory/p571_hsolve.py`, `theory/p571_levels.py` | hierarchical S/D solver (571 session) |
| `report/paley_523_answer.md`, `report/tractability_after_2292.md`, `report/unattempted_2000_3000.md`, `report/slice_screens_20261013.md` | the record |

## 5. Rules

- A construction counts only if it is verified exactly. Separate theorems, proofs and estimates.
- For each route give the exact objects, the bit costs, the smallest checkable instance and the first computation to
  run.
- Calibrate every Gaussian you rely on against exact counts at small orders.
- **Compute:** 4 cores and about 15 GB. Install `ortools` and `numba` with pip if they are missing. Containers
  restart: write restart-safe runners and commit results as they arrive.
- Do not use `pkill -f` with a pattern that matches your own shell; kill by PID.
- Commit trailers as the repository's existing commits use them. No pull requests.
