# Prompt: construct H(2396) (v = 599), led by the b = 5, c = 6 Paley(11) slice (self-contained, maximum reasoning)

*For a maximum-reasoning Claude Code session on `middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`.*

- **Tags:** [M] verified by exact computation in the repo (name the script), [P] proved (give the proof),
  [E] Gaussian estimate (name the estimator and its calibration), [L] literature recollection (check it; arXiv and
  OEIS are blocked from these containers).
- **Branch:** other sessions push to it. Fetch and merge before every push. Never kill processes you did not start,
  and never `pkill -f` with a pattern that matches your own shell; kill by PID.

## 0. Goal and deliverables

Construct a Hadamard matrix of order **2396 = 4·599**. 599 is prime and the order is open.

This repo has built six previously open orders:
- H(2092), H(2292), H(2588): b = 1;
- H(2884), H(3668): Cooper–Wallis;
- H(3732): b = 0.

It has also built a new H(1388) in the published b = 5 slice of 347 as a positive control.

Deliverables, in order of value:
- **(A)** A verified H(2396): H Hᵀ = 2396·I, checked by an independent script that shares no code with the search
  (patterns: `scripts/verify_packed_h.py`, `scripts/verify_h2092.py`). Also a bundle `pub/H2396/`: core on G,
  `border.json` with layout, K, Kc and the corner written out, the normalized matrix, and a standalone `verify.py`.
  The b ≥ 3 bundle writer of the H(1388) control does this (commit d7c7c6a; see `report/alpoge_design_answer.md` §5).
- **(B)** For the lead slice: a stage-1 generator that delivers bordered y's at a measured rate, plus a calibrated
  stage-2 cost per core (≤ about 2^45 operations, ≤ about 100 core-hours on 4 cores), plus the production run.
- **(C)** Exact obstructions, slice by slice.
- **(D)** If the lead fails, a ranked list of the other 599 routes in §4 with exact status and next computation.

Write to `report/b599_answer.md`, with scripts in `theory/b599_*.py` and logs in `logs/b599/`. Commit and push as
results arrive; containers restart. Commit trailers follow the repository's existing commits; no pull requests.

## 1. The b ≥ 3 construction (exact; `report/alpoge_design_answer.md` §1, `report/bge3_construction.md`)

- **Parameters:** c = b + 1 even, c | v + 1, h = (v + 1)/c − 1, G abelian of order m = c·h, H ≤ G of index c,
  Q = G/H.
- **Core:** four ±1 functions x_k on G with
  - Σ_k PAF_k = T, where T(0) = 4m, T = −4b on H \ 0 and T = +4 off H;
  - in characters, Σ_k |X̂_k(χ)|² = 4v, except 4b when χ ≠ 1 is trivial on H.
- **Border:** a frame K (4b × 4c, rows balanced in each of the four blocks, K Kᵀ = 4c·I) inflated by h, and the
  same for Kc. The corner is forced, C = −K·Cb(z)ᵀ·Kcᵀ/(4c), where z is the core's Q-compression, and C must be ±1.
- **Border classes for c = 6:**
  - the 347 template class, which transfers to 599: its row equation has 528–1,056 solutions, every one checked gives
    a ±1 corner [M], §3;
  - the monomial / Singer class (`theory/ab_mono_border.py`, §3.2, §11.3).

  `theory/ab_borderfit.py` decides whether a given y's compression is bordered by a template frame.

## 2. The lead: v = 599, b = 5, c = 6, q = 11, s = 2, H′ = Z₉, U = −1 on Z₉

**Group.**
- m = 594 = 2·3³·11 and h = 99.
- G = Q × H with Q of order 6 and H = H′ × Z₁₁, H′ = Z₉. G′ = Q × H′ has order 54.

**Slice.**
- Paley(11) on Z₁₁: index 2, f = 5, and 11 ≡ 3 (mod 4).
- U = −1 on Z₉ on G′.
- Per fibre over b ∈ G′ the labels are α(b) at z = 0, p(b) on the QRs and n(b) on the non-residues. That gives
  30 G′-orbits and 90 orbits per block.

**Stage 1 [P]** (§1.2): y = α + 5(p + n) ∈ {±1, ±9, ±11} on the 30 G′-orbits, with

    Σ_k PAF_{G′}(y_k) = 2176 (b = 0),   −220 on H′ \ 0,   +44 off H′.

- **Lag 0:** Σ_{k,b} y² = 2176 over 216 entries.
- **Q-level:** |ẑ_k(ψ)|² ≤ 4b = 20 for ψ ≠ 1, and Σ_k a_k² = 4v = 2396.
- **α-congruence:** mod 4f₀.
- **S/D split** on the Z₂ of Q.

**Stage 2 [P].** For an index-2 slice with q ≡ 3 (mod 4) it is a PAF problem on G′ plus a linear system (§11.3).
Write u = p − n and a = α − σ/2. Given y, the ±2-valued u_k on the free points must satisfy:
- Σ_k PAF_{G′}(u_k) = prescribed;
- Σ_k CC_odd(a_k, u_k) = 0.

The 523 machinery (S/D split on a Z₂ factor, middle levels, compiled D-check) applies unchanged (§1.3). The
generic production script is `theory/ab_prod_generic.py`.

**Budget [E]** (`theory/ab_levels_batch.py`, typical-y; §11.2):
- whole slice +28.8; stage 1 +48.1; stage 2 −2.2 per y (C₂ 108);
- **D_def 21.5**: about 2^19.3 S-solutions per y and **0.22 expected cores per y**;
- one middle level: F = 43.3, sides of 2^23 per middle value, about 2^21 pair visits per S-solution;
- **about 2^43 per core**, the size of the 523 run.
- Calibration of this cost model: 347 had D_def 16.2 and took 28 minutes; 523 had D_def 26.1 and cost about 2^44.

**What is decided [M].**
- The same G′ and q with U = ⟨2⟩ (order 6) on Z₉: 231,363 α-quadruples, 0 S-solutions. Empty in both border
  classes.
- Five other multiplier groups (§7): stage 1 empty.
- The ⟨2⟩ sub-slice's stage-1 solutions would have been stage-1 solutions of the lead. **The lead's stage 1 is
  undecided.**

**Why it is stuck.** With 30 G′-orbits per block (2^30 α-vectors per block), the α-quadruple-first enumerator
(`theory/ab_s1sd.py`) does not scale.

**But enumeration is not required.** At about 0.22 cores per y, a core needs only about 5–50 bordered y's.
- Stage 1 is estimated at +48 bits [E]: abundant if the estimate holds. Still, stage-1 estimates have been wrong in
  both directions, so this must be checked exactly.
- 523 was solved the same way: a residue-aware 4-list join generated y's at a few per minute
  (`theory/q523_ymitm.py`) after naive joins and annealing had failed. They failed because of exact mod-8 and PSD
  structure, which is the α-congruence and the Q-level.

## 3. Plan (reorder if you find better)

1. **Stage-1 generator.**
   - Build a sampler or join that outputs exact stage-1 solutions of the lead, bordered from the start: Q-compressions
     restricted to rows admitted by screen A and fitted by a 347-class or monomial border.
   - Ingredients to combine:
     - the α-congruence residues;
     - the Q-level per-character bound;
     - the S/D split on Z₂ ⊂ Q;
     - per-block spectral bounds applied during list recursion (the report's suggestion, §11.4);
     - energy-split hashed joins (`theory/nd_ejoin.py`);
     - and the `q523_ymitm.py` pattern: random block lists with the composition and rows built in, a middle value on a
       few lags, a 2+2 hash on the rest.
   - Validate on the 347 slice (stage 1 must contain 347's own y) and on small same-shape slices with known counts
     (`theory/ab_exact_small.py`).
   - Alternatively, prove the lead's stage 1 empty.
2. **Calibrate stage 2** on the first real y's.
   - Exact S-solution counts per middle value against the Gaussian.
   - D-success per S-solution against 2^{h_D − C(D|S)}, measured as the 523 session did
     (`theory/q523_dfine_calib.py`).
   - Use the index-2 linear structure (§2) if it shrinks the D-check or the S-join.
   - Bias the S-join toward S-solutions with large h_D: the cost per core is about 2^{C(D|S) − h_D/2}.
3. **Cost, then production** on 4 workers, restart-safe, with a tally of expected cores. H(2092) arrived at about
   3 expected cores, a 5% wait; do not stop a calibrated run early. Any core: assemble with the border and forced
   corner, verify independently, commit at once, then bundle.
4. **If the lead dies or costs more than 2^47, work through §4.**

## 4. Other routes to 599 (record exact status for each one you touch)

| route | status | next |
|---|---|---|
| b = 5, c = 6, H′ = Z₃² or other U | D_def 30–32 [E], stage 1 +88 | only if the lead's stage 1 is empty |
| b = 5, c = 6, Z₁₁ index 5 or 10 (f = 2, 1) | not costed | screens A and B, level budget |
| b = 5, c = 6, uniform cyclotomy over F₉ on Z₃² ≤ H (the F_{p²} family, `theory/nd_fq*.py`, `report/next_design_answer.md` §4) | not tried | derive targets, screen |
| b = 1, doubly staged slices on m = 598 = 2·13·23 | **not screened**. Same shape as 647's 646 = 2·17·19, which built H(2588). Whole-slice Gaussians mean nothing for these (§9.4 there). | the sweep's slices acting on both Z₁₃ and Z₂₃ (`logs/msweep_record/m598.jsonl`; whole-slice log₂E in brackets): ×3 (order 33, 30 orbits, −3.6), ×49 (18 orbits, −27.9), ×31 (24, −33.6), ×25 (42, −49.8), ×45 (70, −108.6), ×51 (40, none), ×{45, 47} (48, none), ×{183, 231} (72, none). The ×353 slice of 647 had a whole-slice +36 but an exact joint result unrelated to it, so these numbers decide nothing. Run the exact joint test `theory/nd_joint_exhaust.py 598 μ 13 23`; then the cross level as at 647 (`report/paley_647_answer.md`, `theory/p647_*.py`) |
| b = 1 GP slices of 598 | Paley(23) ×27 dead (w-screen); others screened in `report/slice_screens_20261013.md` and `report/next_design_answer.md` §7 | re-check coverage |
| b = 3, c = 4 | every configuration dead (screen A, stage 1; §7, §11.3) | – |
| b = 7, c = 8 (h = 74 = 2·37) | q = 37, s = 4: screen A fails | other s |
| c = 10 (h = 59), 12 (h = 49 = 7²; F₄₉ cyclotomy), 20 (h = 29), 24 (h = 24), 30 (h = 19) | frames exist for every even c (§3.1); no corner known for c ≥ 10 (the Singer ansatz is exhausted at c = 10) | the corner-flatness condition of §3.1 for a non-Singer frame; only if a much cheaper core level appears |
| b = 0 over Z₅₉₉, T-matrices | b = 0: row-sum lemma; T(599): no route (`report/tmat_routes.md` §10) | – |

## 5. Tools (read before writing new code)

| file | what it does |
|---|---|
| `theory/ab_lib.py`, `ab_stage1.py`, `ab_stage2.py`, `ab_s1sd.py`, `ab_s1_direct.py`, `ab_s1_sample.py`, `ab_s1_bordered.py`, `ab_s1_anneal.py` | staged-system library (any abelian G, b = 0/1/≥ 3 targets); exact and sampled stage 1 |
| `theory/ab_sweep.py`, `ab_cbatch.py`, `ab_c6_batch.py`, `ab_exact_small.py` | screens A and B, exact batches, brute-force validation (0 false negatives on 100 slices) |
| `theory/ab_levels.py`, `ab_levels_batch.py`, `ab_gauss.py`, `ab_typical_levels.py` | level budgets, D_def, lattice Gaussians |
| `theory/ab_prod347.py`, `ab_prod_generic.py`, `ab_small_exact.py` | hierarchical stage-2 production (347 control: new H(1388) in 28 minutes); exact small stage 2 |
| `theory/ab_assemble.py`, `ab_template.py`, `ab_border_transfer2.py`, `ab_borderfit.py`, `ab_mono_border.py`, `ab_singer.py` | assembly, templates, border transfer and fit, monomial / Singer borders |
| `theory/q523_*.py`, `theory/p647_*.py` | 523 and 647 production kernels (Wagner tree, `mjoin_rad2`, compiled D-check, checkpoints). Read only; copy into `theory/b599_*` to adapt. |
| `theory/nd_joint_exhaust.py`, `nd_cross.py`, `nd_ejoin.py`, `nd_fq*.py` | doubly staged b = 1 tools; F_{p²} tools |
| `logs/ab/p1/`, `logs/ab/mono6/`, `logs/ab/c6_*` | existing 599 c = 6 logs |

## 6. Rules

- A construction counts only if it is verified exactly. Separate theorems, proofs and estimates.
- Calibrate every Gaussian against exact counts. Lattice and 2-adic obstructions are real.
- **Compute:** 4 cores and about 15 GB of RAM. Install `ortools`, `numba` and `sympy` with pip if missing.
- Write restart-safe runners. Arm a check-in so the session never idles while workers run (an idle session was
  reclaimed once and killed production workers). Commit results as they arrive.
