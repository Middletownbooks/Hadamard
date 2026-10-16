# Twisted (block-permuting) symmetries in bordered b = 1 cores at the open orders

Date: 2026-10-16. Branch `claude/focused-planck-j1645f`. Answers `report/prompt_b1_twist_design.md` (item 2).
Scripts: `scripts/b1tw_*.py`. Logs: `logs/b1tw/`.

**Tags.** [M] exact computation (script named); [P] proved here; [E] Gaussian estimate (estimator and calibration
named).

---

## 0. Bottom line

1. **(A) No H(4v) was built.** No twisted b = 1 slice at an open order is within reach. Every slice that survives the
   exact screens costs n ≥ 155 constraint bits; the b = 0 route that just worked at 883 had n = 112 (§4, §5).

2. **(E) The parity-character lemma [P]** (§2) closes the 3-cycle twist at eight open orders for every group and
   every automorphism.
   - Let |G| = v − 1 = 2h with h odd. A b = 1 core whose blocks B, C, D are cycled by an automorphism needs
     v = 4a² + 3e² with e odd.
   - That rules out **599, 659, 719, 839, 863, 947, 955 and 971**.
   - At the other eight orders the sum of B over the odd-order part H is pinned: |e| = 13, 9, 15, 5, 7, 5, 13, 13
     at 571, 643, 739, 751, 823, 859, 907, 991.
   - This is the b = 1 counterpart of 919's b = 0 coincidence 4·919 = 1 + 3·35². The prompt's premise that "row sums
     cost nothing" at b = 1 was wrong: the coincidence moves from the row sums to the order-2 character.

3. **(B) Estimator and validation** (§3).
   - `scripts/b1tw_lib.py` handles twisted b = 1 slices on Z_m and `scripts/b1tw_ab.py` on any abelian G: orbits with
     stabilisers, twists 'BC' and 'BCD', exact row-sum DP and sampling, and a Gaussian on the classes.
   - `scripts/b1tw_exact.py` gives exact counts, cross-checked by an independent direct-PAF counter.
   - **The Gaussian alone is unreliable on twisted slices.** It predicted 2^8.7–2^10.9 solutions in slices with exactly
     none (30/×7 and 30/×17, 'BC'), and ran 3 bits low (42/×11 'BC': 2^9.9 against exactly 2^12.9).
   - **The exact quotient screens are decisive.** They are sound (proof) and, on 62 small slices, exactly as strong
     as the exact count: 60 killed, all with 0 solutions; 2 kept, both with solutions.

4. **(C) Screen of the open orders** (§4).
   - **Cyclic G = Z_{v−1}:** 551 slices (205 'BC', 117 'BCD', 229 untwisted for reference); all 19 orders;
     8–150 orbits.
     - The exact quotient screens kill 410 slices, 229 of them at the parity quotient Z₂.
     - Surviving twisted slices with log₂E > 0: the cheapest is 571, u = 103 (d = 12), 'BC': 66 orbits, 37 classes,
       log₂E = 42.5, **n = 198**.
   - **Noncyclic G:** Z₂² × Z₁₆₇ (669), Z₂² × Z₁₉₇ (789), Z₂³ × Z₆₇ (537) and Z₂ × Z₃² × Z_q (739, 955, 991),
     with GL(2,2), GL(3,2) and GL(2,3) automorphisms.
     - 470 slices; the exact screens kill 357.
     - The cheapest twisted survivors are 789's GL(2,2) 3-cycle ('BCD', 58 orbits, 29 classes, log₂E = 8.7,
       **n = 155**) and 739's order-4 GL(2,3) element with μ ('BC', 46 orbits, 29 classes, log₂E = 7.5, n = 156).

5. **(D)** No search was run: nothing is within budget (§5).

6. **Side check: b = 0 at the composite orders** (§6).
   - 537, 669, 789 and 955 with multipliers acting on both factors; nothing new.
   - The one cheap slice, 955 'BC' with n = 89.7, is exactly empty: its Z₅ quotient has no solution.

## 1. Twist types and containment [P]

- **Twist types.** The ansatz is X_{π(i)}(u x) = ε_i X_i(x) with u ∈ Aut(G).
  - 0 is a fixed point of u, so its label must be g-fixed. That forces the signs around every cycle of π to multiply
    to +1, and forbids ε = −1 on any π-fixed block.
  - Negating blocks therefore removes all signs.
  - With the ±2 block fixed (the prompt's restriction), the twists are:
    - 'id', the multiplier slice;
    - 'BC', with C(ux) = B(x) and B(ux) = C(x), needing even d;
    - 'BCD', with C(ux) = B(x), D(ux) = C(x) and B(ux) = D(x), needing 3 | d.
  - Slices depend only on ⟨u⟩, up to relabelling C ↔ D.
- **Containment.**
  - In a 'BC' slice every block is u²-invariant, so the slice lies inside the untwisted ⟨u²⟩ slice. Likewise 'BCD'
    lies inside ⟨u³⟩.
  - Bits and classes: a 'BC' slice has the same number of free bits (4 per U-orbit, counting B's two U²-halves) and the
    same classes (⟨U, −1⟩-orbits) as the untwisted ⟨u⟩ slice.
  - **So twisting does not lower the Gaussian cost.** It selects a different family at the same cost, with different
    exact (quotient) behaviour.

## 2. The parity-character lemma [P]

**Setting.** Let |G| = m = 2h with h odd, so G = ⟨t⟩ × H with |H| = h. This covers every prime v ≡ 3 (mod 4), 955,
and every G of order v − 1 for those v. Let χ be the order-2 character: χ(t) = −1 and χ|_H = 1. For a block X write
e_X = Σ_{x∈H} X(x) and o_X = Σ_{x∈tH} X(x). Both are odd, since each sums h signs. Then X̂(χ) = e_X − o_X and
s_X = e_X + o_X.

**Lemma.** For a b = 1 core (Σ_X |X̂(ψ)|² = 4v for ψ ≠ 1, block sums (±2, 0, 0, 0)):

    v = 4a² + e₁² + e₂² + e₃²,   with e₁, e₂, e₃ the H-sums of the three zero-sum blocks (all odd).

*Proof.*
- A zero-sum block has o_X = −e_X, so X̂(χ) = 2e_X ≡ 2 (mod 4).
- The ±2 block has X̂(χ) = 2e_X ∓ 2 = 4a, since e_X is odd.
- Then 4v = 16a² + 4(e₁² + e₂² + e₃²). ∎

**Corollary 1 (3-cycle twist).**
- Every automorphism maps H, the elements of odd order, onto itself.
- If u cycles B → C → D, then e_C = Σ_{x∈H} C(ux) = Σ_{x∈H} B(x) = e_B, and likewise e_D = e_B.
- Hence **v = 4a² + 3e², with e odd.**
- For prime v ≡ 3 (mod 4) this holds if and only if v ≡ 1 (mod 3):
  - such a prime is x² + 3y², and v ≡ 3 (mod 4) forces x even and y odd;
  - a prime v ≡ 2 (mod 3) is not of the form x² + 3y².
- For 955 = 5·191 there is no representation either: 5 ≡ 2 (mod 3) divides it to an odd power.
- So **599, 659, 719, 839, 863, 947, 955 and 971 have no 3-cycle-twisted b = 1 core** for any G of order v − 1 and
  any automorphism u.
- At 571, 643, 739, 751, 823, 859, 907 and 991 the representation is unique (`logs/b1tw/parity_lemma.txt`). It pins
  |e_B| to 13, 9, 15, 5, 7, 5, 13 and 13.

**Corollary 2 (transposition).**
- 'BC' gives e_B = e_C, so v = 4a² + 2e_B² + e_D².
- This form represents every open v, with 3 to 15 representations each, so there is no global obstruction.
- The kills at the parity quotient come instead from the slice's letters. With orbits of size d on H ∖ 0, e_X is
  confined to X(0) + (d or d/2)·Z, and for large d no representation fits.

**Check against the screen [M]** (`logs/b1tw/screen_open_w*.jsonl`).
- All 31 cyclic 'BCD' slices in range at the lemma-closed orders with such slices (599, 659, 955, 971) die at Z₂, as
  the lemma requires.
- At 537 (15 slices) and 907 (8) they die at Z₂ by their letters.
- 571, 739, 859 and 991 keep survivors.

**For 4 | m** (537, 669, 789) with G ⊇ Z₂², the three order-2 characters together give
v = a² + 4(b₁² + b₂² + b₃²) with a odd, under a GL(2,2) 3-cycle. That is always solvable, so the twist is compatible
there. This is the 479 mechanism.

## 3. Stage A: estimator and validation

**Tools.**
- `scripts/b1tw_lib.py` (G = Z_m) and `scripts/b1tw_ab.py` (any G = product of coprime cyclic and elementary
  components, `theory/abgroup.py`).
- Slice construction covers orbits of any size, with stabilisers: odd-size orbits force B = C under 'BC', and sizes
  prime to 3 force B = C = D under 'BCD'.
- Exact row-sum counts (DP) and exact uniform sampling conditioned on the row sums.
- **Gaussian.** It runs on the class sums. The class {m/2}, of size 1, is dropped because its value is implied by the
  row sums. Lattice steps per class are measured from the samples. The output is log₂E for the whole slice (all
  placements of the ±2 block) and the conditional constraint bits n.
- **Exact counter** (`scripts/b1tw_exact.py`): complete enumeration of each row-sum group, joined by 64-bit keys, every
  join re-checked. An independent direct-PAF counter (Python dicts, all lags) agrees: 7,680 at 42/×11 'BC' and 0 at
  30/×17 'BC'.
- **Exact quotient screens** (`scripts/b1tw_lag0.py`, `scripts/b1tw_ab.py` [P]).
  - Compress each block onto G/H for every u-stable subgroup H. Under 'BC' and 'BCD', C and D compress to B's
    compression composed with u⁻¹ and u⁻².
  - The compressed blocks must satisfy Σ_k PAF_{y_k}(j) = −4|H| + [j = 0](4|G| + 4) on G/H, with the row sums.
  - The all-lags test enumerates every attainable compressed block (cells × letters) and joins them. When that is too
    large it falls back to the lag-0 energy DP.

**Validation on small m** (`logs/b1tw/validate_small.jsonl`, `logs/b1tw/calib_small.jsonl`, `logs/b1tw/calib_alive.jsonl`).
- **The quotient screens are sound and, at small m, complete [M].** Over 62 slices with m ≤ 186:
  - 60 are killed, every one with exact count 0;
  - the 2 kept (42/×11, 'id' and 'BC') have 24,576 and 7,680 solutions.
- **Twisted slices are almost always empty at small m [M].** 291 small slices with m ≤ 326 and every
  row-sum group of at most 20 bits were decided exactly. Solutions occur in 3:
  - 2 of 167 untwisted slices;
  - 1 of 100 'BC' slices;
  - 0 of 24 'BCD' slices.
- **The Gaussian is not usable on its own [M]:**

| slice (m, u, twist) | orbits | classes | Gaussian log₂E | exact | exact quotient screen |
|---|---|---|---|---|---|
| 30, ×7, BC | 12 | 7 | +8.7 | 0 | killed (Z₆ energy) |
| 30, ×17, BC | 10 | 7 | +10.9 | 0 | killed (Z₁₅ all lags; energy passes) |
| 66, ×5, BC | 12 | 7 | +3.9 | 0 | killed (Z₃₃ all lags) |
| 70, ×3, BC / BCD | 10 | 7 | +2.8 / +2.2 | 0 / 0 | killed |
| 38, ×7, id | 14 | 7 | +29.5 | 2^28.95 | alive |
| 42, ×11, id | 12 | 7 | +9.6 | 2^14.58 | alive |
| 42, ×11, BC | 12 | 7 | +9.9 | 2^12.91 | alive |

  - Where solutions exist the Gaussian runs 0.6 bits high to 5 bits low, the same spread as the 3-cycle b = 0 family
    (`report/b0_919_answer.md` §3).
  - Where the exact screens kill a slice, the Gaussian still predicts up to 2^11 solutions.
  - **Rule used in Stage B:** exact quotient screens first; the Gaussian only on survivors, read as ±5 bits.

## 4. Stage B: the screens

### 4.1 Cyclic G = Z_{v−1}, all 19 open orders [M]

Source: `scripts/b1tw_screen.py`; logs `logs/b1tw/screen_open_w{0,1,2}.jsonl`, summary in
`logs/b1tw/screen_open_summary.txt`.
- **Scope.** Every cyclic subgroup ⟨u⟩ of (Z/m)* with 8–150 orbits, in 'BC' (d even) and 'BCD' (3 | d), plus the
  untwisted slices for reference.
  - 551 slices: 205 'BC', 117 'BCD', 229 untwisted.
  - The exact quotient screens kill 410, 229 of them at the parity quotient Z₂.
- **Per order** (twisted slices only):

| v | twisted slices | killed at Z₂ | alive, log₂E > 0 | cheapest alive twisted (u, d, twist, orbits, classes, log₂E, n) |
|---|---|---|---|---|
| 537 | 36 | 29 | 6 | ×97, 6, BC, 96, 59, +75.3, 285.9 |
| 571 | 45 | 0 | 13 | ×103, 12, BC, 66, 37, +42.5, **198.2** |
| 599 | 26 | 18 | 2 | ×321, 6, BC, 120, 71, +69.3, 384.7 |
| 643, 659, 669, 751, 789, 839, 863, 955 | 3, 14, 2, 9, 9, 2, 2, 20 | most | 0 | none |
| 739 | 32 | 10 | 4 | ×215, 10, BC, 90, 49, +56.9, 267.5 |
| 823 | 9 | 7 | 2 | ×41, 8, BC, 106, 71, +44.7, 361.8 |
| 859 | 40 | 0 | 2 | ×289, 15, BCD, 90, 45, +23.5, 288.2 |
| 907 | 20 | 17 | 2 | ×283, 10, BC, 96, 63, +35.5, 327.8 |
| 947 | 4 | 0 | 1 | ×21, 14, BC, 84, 47, +27.2, 276.2 |
| 971 | 17 | 7 | 6 | ×469, 12, BC, 86, 53, +43.1, 282.8 |
| 991 | 32 | 0 | 1 | ×71, 10, BC, 150, 79, +24.3, 527.6 |

- **Why the cheap twisted slices die [P + M].**
  - n is about 5.3 bits per class, and a slice with few classes has large orbits.
  - Large orbits make the letters at the parity quotient coarse: e_X ∈ X(0) + d·Z, or + (d/2)·Z for B under 'BC'. The
    lemma's representation (§2) then has no admissible solution.
  - Composite m adds small-layer orbits, so every slice carries extra classes: about 37 classes for 66 orbits at 571,
    against 21 for 43 orbits at the prime 883.
- **Reference.** The untwisted survivors of my screens include 599 ×3 (n = 94) and 537 ×9 (n = 144). Both were
  already proved empty by the repo's stronger exact compressions (599 ×27 by the w-screen, `report/b599_answer.md`;
  537 ×9 stage 1, `report/tractability_after_2292.md`).
  - So my quotient screens are necessary conditions, weaker than CP-SAT compressions.
  - The twisted survivors above may be dead too. At n ≥ 198 that no longer matters.

### 4.2 Noncyclic G [M]

Source: `scripts/b1tw_ab_screen.py`; logs `logs/b1tw/screen_ab*.jsonl`, summary in `logs/b1tw/screen_ab_summary.txt`.
- **Scope.** u runs over every cyclic subgroup of Aut(G) = ∏ Aut(component), with GL parts up to conjugacy. The screens
  use every u-stable subgroup H, taken as a product of u-stable component subgroups. All-lags tests run when
  |G/H| ≤ 64, lag-0 energy above that.

| v | G | slices (id / BC / BCD) | killed | at a 2-group quotient | alive with log₂E > 0 | cheapest alive twisted (d, twist, orbits, classes, log₂E, n) |
|---|---|---|---|---|---|---|
| 537 | Z₂³ × Z₆₇ | 40 / 27 / 26 | 88 | 57 | 2 twisted, 1 untwisted | 6, BC, 94, 71, +14.8, 341.8 |
| 669 | Z₂² × Z₁₆₇ | 5 / 3 / 1 | 9 | 9 | none | – |
| 739 | Z₂ × Z₃² × Z₄₁ | 64 / 63 / 12 | 99 | 57 | 15 twisted, 9 untwisted | 20, BC, 46, 29, +7.5, 156.4 |
| 789 | Z₂² × Z₁₉₇ | 18 / 14 / 6 | 35 | 35 | 1 twisted, 1 untwisted | 21, BCD (GL(2,2) 3-cycle × μ of order 7), 58, 29, +8.7, **154.7** |
| 955 | Z₂ × Z₃² × Z₅₃ | 33 / 31 / 8 | 60 | 45 | 4 twisted, 5 untwisted | 8, BC, 134, 81, +72.0, 445.3 |
| 991 | Z₂ × Z₃² × Z₅ × Z₁₁ | 54 / 53 / 12 | 66 | 0 | 6 twisted, 6 untwisted | 8, BC, 144, 79, +106.4, 449.8 |

- **669.** 166 = 2·83 leaves only slices of 4–12 orbits, which die at Z₂ or Z₂², or slices of 168 or more orbits.
- **The GL(2,2) 3-cycle** (Alpöge's 479 mechanism, §2) is compatible at 789 and survives the screens. It still needs
  n ≈ 155 in a two-list problem.
- **GL(2,3) on Z₃²** (739, 955, 991) adds slices but no cheap survivor.

## 5. Stage C/D: costs of the best survivors [E]

**Cost scale, calibrated on the 883 run** (`report/b0_883_answer.md`).
- A four-list flat-pool k-tree at n = 112 found its core after about 6 hours on 4 cores, roughly 2^37 work.
- List methods cost about 2^{n/3} with four independent lists and 2^{n/2} with two.

**'BC' (571 ×103, n = 198, log₂E = 42.5).**
- The lists are A and D (U-invariant, 66 bits each) and B (U²-invariant, about 132 bits): a three-list problem.
- Lemma H's P/R split would give four lists, but only per support.
  - The support costs F ≈ 60 bits (the even-size U-orbits), against log₂E = 42.5.
  - Lemma G (F ≤ log₂E) fails: each support holds about 2^−17.5 solutions.
- Even with four lists the cost is ≥ 2^66, about 2^29 times the 883 run.

**GL(2,2) 3-cycle at 789 (n = 155, log₂E = 8.7).**
- Two lists: A with 58 bits and B with 116 bits.
- About 2^77, and the A-table (≈ 2^54 entries) cannot be stored.

**Verdict.** No twisted b = 1 slice at an open order is searchable. Stage D (control and production) was not run.

## 6. Side check: b = 0 slices at the composite open orders [M, P]

Prime-order b = 0 screens exist (`report/b0_919_answer.md` §4). The composite orders were examined only through the 933
dual/cross pipeline (`report/alpoge_design_answer.md` §6), which costs them as not viable. So
`scripts/b0tw_composite.py` screens G = Z_p × Z_q with u = (a, μ) acting on both factors: twists 'id', 'BC', 'BCD';
b = 0 row-sum patterns (Σ s² = 4v, with s_B = s_C [= s_D]); quotient tests onto Z_p and Z_q; the Gaussian.
Logs: `logs/b0tw/composite_{537,669,789,955}.jsonl`, summary in `logs/b0tw/composite_summary.txt`, up to 70 orbits.

- **537 and 789.** q − 1 = 2·89 and 2·131 leave only slices of 6–9 orbits, which die on row sums, or slices above
  the range.
- **669.** Slices up to 70 orbits die on row sums. The only survivor (μ of order 6, 114 orbits) has n ≈ 390.
- **955.**
  - 'id' with μ of order 19 on Z₁₉₁ (55 orbits, 27 classes): n = 143. That is beyond the flat-pool k-tree, which was
    at its limit at 883 (n′ ≈ 79 with 2^24 pools).
  - **'BC' with u = (−1 on Z₅, μ of order 19)**: 33 orbits, 17 classes, Gaussian log₂E = +18.3, **n = 89.7**. This
    was the cheapest slice found at any open order.
    - **It is exactly empty [P + M]** (`scripts/b0tw_955_z5.py`, `logs/b0tw/v955_bc_z5_quotient.txt`).
    - Compress along Z₁₉₁. Then y_A and y_D are symmetric on Z₅, y_C is y_B reversed, and every entry lies in
      ±1 + 38ℤ.
    - The Z₅ condition PAF_{y_A} + PAF_{y_D} + 2PAF_{y_B} = 3820·δ, with the slice's row sums (37, 35, 1) or
      (1, 35, 37), has no solution among the 216 × 216 × 352 admissible vectors.
- **So the side check adds no candidate.** Its one cheap slice is another case of a positive Gaussian killed by a
  small quotient with coarse letters, the pattern of §3.

## Reproduce

```
python3 scripts/b1tw_exact.py 42 11 BC 30 17 BC 30 7 BC            # exact counts + Gaussian (7680, 0, 0)
python3 scripts/b1tw_lag0.py 30 17 BC 42 11 BC                      # quotient screens (lag-0 energy) per slice
python3 scripts/b1tw_validate.py OUT.jsonl 30 42 54 66 70 78 90     # small-m exact vs Gaussian
python3 scripts/b1tw_calib.py OUT.jsonl CAL.jsonl                   # + all-lags quotient verdicts
python3 scripts/b1tw_calib_alive.py OUT2.jsonl $(seq 30 4 330)      # exact counts of every small slice
python3 scripts/b1tw_screen.py S.jsonl 571 643 739 --rmax 150 --id  # cyclic screen of open orders
python3 scripts/b1tw_ab_screen.py N.jsonl 789 E2.2+C197 --rmax 150 --id   # noncyclic (Z_2^2 x Z_197)
python3 scripts/b1tw_ab_screen.py N.jsonl 739 C2+E3.2+C41 --rmax 150 --id # noncyclic (Z_2 x Z_3^2 x Z_41)
```
