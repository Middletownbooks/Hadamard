# Prompt: construct H(2092) (v = 523) in the quartic-Paley slice (self-contained, for a fresh Claude Code session)

*For a maximum-reasoning Claude Code session on `middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`.
Tags: [M] verified by exact computation in the repo, [P] proved, [E] Gaussian estimate, [L] literature recollection
(check it; arXiv and OEIS are blocked from these containers).*

## 0. Goal and deliverables

Construct a Hadamard matrix of order **2092 = 4·523** (523 is prime). It is open.
- 2092 is the smallest order on our list, and competitors are active: every open order below 2000 was closed in
  August 2026.
- We built H(2292) and H(2884) in this repo (`report/RESULT_2292.md`, `cw/RESULT_2884.md`).

The target is a b = 1 bordered Goethals–Seidel core in a **quartic generalisation of the Paley-coset slice**:
- its first stage has already been solved exactly;
- its second stage is the open problem (§3, §4).

Deliverables, in order of value:
- **(A)** a verified H(2092): the core, the border and H Hᵀ = 2092·I exactly (pattern: `scripts/verify_h2292.py`);
- **(B)** a stage-2 method costing ≤ about 2^45 per core [E], validated on a smaller order of the same shape (§5);
- **(C)** an exact obstruction proving the slice (or every y of it) empty.

Write to `report/paley_523_answer.md`, with scripts in `theory/q523_*.py`. Commit and push.
- Mark claims [M] (naming the script), [P] (with a proof) or [E].
- Do not kill processes you did not start: slice sweeps for other orders may be running (`theory/p571_msweep.py`).

## 1. Architecture (exact)

**Bordered GS, b = 1.**
- Four ±1 functions x_k on an abelian group G of order m = v − 1 = 522, with Σ_k PAF_k(g) = −4 for g ≠ 0 and row sums
  (2, 0, 0, 0) after reordering and negating blocks.
- Completion: the GS array of the group-developed blocks, plus border rows and columns S ⊗ 1ᵀ and Tᵀ ⊗ 1 (S, T ∈ H(4))
  and the forced corner −S Σᵀ T/4.
- `theory/paley_b1pipe.py` (`finish`, `bordered_gs`) and `scripts/verify_h2292.py` (cyclic version) do this.

**The slice (523 ×181).**
- G = Z₅₂₂ ≅ Z₁₈ × Z₂₉ by CRT. The multiplier 181 satisfies 181 ≡ 1 (mod 18) and 181 ≡ 7 (mod 29).
- So the slice group M = ⟨181⟩ (order 7) **fixes Z₁₈ pointwise** and acts on Z₂₉ by H = ⟨7⟩, the 7th powers: the
  index-4 subgroup of (Z/29)*.
- Cyclotomic classes: C_j = 2^j·H for j = 0..3 (2 is a primitive root mod 29), each of size f = 7. −1 = 2¹⁴ lies in
  C₂.
- On each Z₂₉-coset g ∈ Z₁₈ of each block k, x is: α at z = 0, and a constant c_j = ±1 on each class C_j. That is
  5 signs per coset, 90 orbits per block and 360 free bits.

**Exact two-stage split.**
- **Stage 1** = the Z₂₉-coset sums y_k(g) = α + 7·Σ_j c_j ∈ {±1, ±13, ±15, ±27, ±29}. They satisfy exactly:

      Σ_k PAF_{Z18}(y_k)(h) = 4m − 4(q − 1) = 1976 at h = 0,   = −4q = −116 for h ≠ 0.

- **Stage 2:** given y, choose (α, c_j) per coset with α + 7Σc_j = y(g), such that for every h ∈ Z₁₈ and class t:

      R(h, t) = Σ_k Σ_g [ α_g c_{g+h,t} + c_{g,t′} α_{g+h} + Σ_{i,j} N_t(i, j) c_{g,i} c_{g+h,j} ] = −4,

  where:
  - t′ is the class of −s for s ∈ C_t, so t′ = t + 2 mod 4;
  - N_t(i, j) = #{z ∈ C_i : z + s ∈ C_t-rep shifted into C_j} are the order-4 cyclotomic numbers of F₂₉.
  The lags (h, 0) follow from stage 1. There are 36 independent equations after R(−h, t) = R(h, t′).
- Verified [M]: `scripts/cyclo_stage2_cpsat.py` uses exactly this reduction. Given the y of a verified H(164) core
  (v = 41, Paley(5) on Z₈), it rebuilds a full core, FFT-checked on Z₈ × Z₅.

**Algebra (for design).**
- span{δ, 1_{C₀}, …, 1_{C₃}} is the Schur ring of order-4 cyclotomy on Z₂₉.
- Its characters are the trivial one plus 4 orbits of nontrivial characters of Z₂₉ under H. Conjugation pairs orbit
  t with t + 2.
- The non-trivial part lives in the cyclic quartic field of conductor 29 (the Gauss periods of order 4).
- So stage 2 is a family of norm-sum conditions, one per character of Z₁₈ and per conjugate pair of period orbits, with
  values in that quartic field tensored with Q(ζ₁₈).
- This is the analogue of the Q ⊕ Q(√−p) split that made the Paley case work (`report/paley_design_answer.md` §2).

## 2. What is established

**Exact screens** [M]:
- The plain Paley(29) slice of 523 (×91) is dead by the w-screen (`theory/paley_wscreen.py`).
- Every slice mapping onto (Z/29)* is dead by the θ₂₉ criterion (`scripts/theta_screen.py`).
- The quartic slice is covered by neither.

**Gaussian estimates** [E]:
- `theory/p571_msweep.py 522` gives the quartic slice ×181 **log₂E = +92.1** (constraint bits 247).
- `scripts/cyclo_split.py 18 29 4 92.1` gives **stage 1 = +44.8, stage 2 = +47.3 per y**.
- So stage 2 is rich: each y should have many completions. That is the opposite of 571, where stage 2 was starved.
- Calibration warnings:
  - the same estimator gave +7.2 for 537's sextic slice, which is exactly empty;
  - it gave positive values for Paley slices that are exactly dead;
  - it undercounted 573's stage 2 by about 2^5–2^10.

**Stage 1 is non-empty [M].**
- `scripts/cyclo_stage1_cpsat.py 18 29 4` (OR-tools CP-SAT; install it with `pip install ortools` if the container
  lacks it) found a verified y in 323 s: `logs/cyclo/stage1_cpsat_523.json`.
  - Rows (2, 0, 0, 0).
  - Letters used: 62 × |1|, 6 × |13|, 4 × |15|.
- Controls:
  - v = 41 (Paley(5) on Z₈): solved in 3.8 s;
  - 537's sextic slice: proved INFEASIBLE in 3.5 s, matching the exhaustive enumerator `scripts/cyclo_stage1_exact.py`.
- More y's: rerun with other `--seed` values. One CP-SAT run gives one y, and about 2^45 exist [E].

**Stage 2 is the open problem [M].**
- `scripts/cyclo_stage2_cpsat.py` returns UNKNOWN after 1800 s on the y above.
- A **planted** instance of the same size (random α, c; targets computed from it) is also UNKNOWN after 600 s, while a
  planted v = 41 instance is solved instantly. So CP-SAT is not a stage-2 solver at this size.
- Annealing and tabu on these reduced systems also fail their positive controls (`scripts/slice_sa.py` on 419).

**Stage-2 size for the found y [E]:**
- per coset, |y| = 1 gives 6 (α, c) choices and |y| = 13 or 15 gives 4, so H₂ ≈ 180 bits;
- yield ≈ +47, so C₂ ≈ 133;
- a plain meet in the middle therefore costs about **2^66 per y**.

## 3. Tools in the repo

| file | what it does |
|---|---|
| `scripts/cyclo_stage1_cpsat.py`, `scripts/cyclo_stage1_exact.py` | stage 1 (CP-SAT; exhaustive for small F) |
| `scripts/cyclo_stage2_cpsat.py` | the exact stage-2 reduction (`solve`, `planted`); cyclotomic classes and numbers (`classes`, `cyclo_numbers`) |
| `scripts/cyclo_split.py` | the stage split [E] |
| `scripts/paley_s2join.py` | the exact stage-2 join that built H(2292): numba hashed 2+2 join with **Z_d class splits** (additive compressed equations select compatible class pairs; at 573 this cut 2^34 → 2^23–2^25 pairs per side), signature dedup, prime-residue passes |
| `theory/p571_hsolve.py`, `theory/p571_levels.py` | hierarchical exact solver of the 571 session: when a factor Z₂ is fixed pointwise, the condition splits into **S = x₀ + x₁ and D = x₀ − x₁ levels** with exact parity links, solved level by level; verified H(172), H(316); its cost floor is (C₂ − h_D)/2 |
| `theory/p571_z10profiles.py` | exact rational-character profile enumeration (Đoković–Kotsireas compression) |
| `report/paley_571_answer.md`, `report/paley_design_answer.md`, `report/tractability_after_2292.md` | the theory and the record |

## 4. Questions (reorder freely)

1. **A hierarchical stage 2.**
   - Z₁₈ = Z₂ × Z₉ is fixed pointwise. So, exactly, the condition splits into an S level and a D level on Z₉ × Z₂₉
     (entries in {0, ±2}, with exactly one of S and D nonzero at each point).
   - Z₉ adds the Z₃-compression (rational characters) as further class splits.
   - Design the level order, the class keys and the sublist sizes. The stage-2 slack (+47 per y [E]) lets you take
     random sublists freely; 571 lacked that. Cost it per core.
   - Is ≤ 2^45 reachable? If not, what is the floor?
2. **Use the quartic-field structure.**
   - Can stage 2 be split further by the conjugate pair of period orbits (the cyclic quartic field of conductor 29), the
     way the Paley case split into a rational and a √−p part?
   - For example: a "stage 2a" on the quadratic subfield Q(√29) (the characters of order 2 restricted to the classes,
     i.e. the Paley(29) part c₀ + c₂ versus c₁ + c₃), then a "stage 2b" for the rest.
   - The Paley(29) slice itself is dead. But here the quadratic part is only a projection, so the w-screen does not
     apply. Check what it constrains.
3. **Choose y for stage 2.** Stage 1 is cheap (minutes per y by CP-SAT). Are some y's provably better for stage 2? For
   example, more |y| = 1 cosets give more freedom, and certain placements of the big letters give exact stage-2
   obstructions or simplifications. Give screens.
4. **Alternative split.**
   - Another option is a joint stage that fixes y together with α (or with the Z₂-compression of c), then joins the rest.
   - Any exact reformulation that balances the two stages better is welcome.
5. **739 (H(2956))** has the same shape: Z₁₈ × quartic on Z₄₁ (×127). Gaussian: whole slice +57, stage 1 +29, stage 2
   +29 per y. Note where the method transfers.

## 5. Positive controls (required before claiming a method)

Use small orders with the same shape: a fixed cofactor Z_F and the order-4 cyclotomic classes on a prime q ≡ 1 mod 4,
with v = Fq + 1.
- **Candidates:** q = 13 (f = 3), q = 17 (f = 4), q = 29 (f = 7), with F ∈ {2, 4, 6, 8, …} and gcd(F, q) = 1. For
  example v = 53 (4·13 + 1), 69 (4·17 + 1), 79 (6·13 + 1), 103 (6·17 + 1), 117 (4·29 + 1), 175 (6·29 + 1).
- Show the slice is non-empty at some such v by your method, and verify the full Hadamard matrix.
- Planted tests at full 523 size are also expected, as `cyclo_stage2_cpsat.py planted` does: a method must recover a
  planted stage-2 solution.

## 6. Rules

- A construction counts only if verified exactly (H Hᵀ = N·I).
- Separate theorems, proofs and estimates.
- For each route give the exact objects, their bit cost, the smallest checkable instance and the first computation to
  run.
- Calibrate any Gaussian you rely on against exact counts at small orders. It has been wrong in both directions here.
