# The b ≥ 3 examples (v = 479, 347, 359): construction analysis

Scripts: `scripts/bge3_border.py`, `scripts/bge3_skeleton.py`, `scripts/bge3_skeleton2.py`, `scripts/b359_aut.py`,
`scripts/bge3_components.py`. Logs: `logs/bge3_*.log`, `logs/b359_aut.log`. All statements are measured [M] on the
full matrices, after undoing the hoa64 normalisation on the core rows and columns. After that correction each full matrix
satisfies H Hᵀ = 4v I exactly. My extracted cores match the original published sequences: block A is identical, and
blocks B, C and D each map onto the original under one shared automorphism of G plus a translation and a sign. No
character twist is involved. Placed in my group coordinates, the originals give the two-valued sums {−4b, +4}.

## 1. One architecture, used identically in all three

Notation: c = b+1, H ≤ G of index c, h = |H|, m = v − b = c·h, and v + 1 = c(h + 1).

| v | b | c | h | G | H | G/H | extension |
|---|---|---|---|---|---|---|---|
| 479 | 3 | 4 | 119 = 7·17 | Z2²×Z119 | Z119 (cyclic) | Z2² | split |
| 347 | 5 | 6 | 57 = 3·19 | Z6×Z57 | Z57 (cyclic) | Z6 | split |
| 359 | 7 | 8 | 44 = 4·11 | Z2²×Z8×Z11 | Z44 = 2Z8×Z11 (cyclic) | Z2³ | **non-split** |

Every one of the three matrices has the same four layers:

1. **Corner.** The 4b × 4b corner is itself a Hadamard matrix, of order 12, 20 and 28 respectively.
2. **Border × core.** Each border row, restricted to a core block, is constant on every H-coset. So the border is a
   4b × 4c sign matrix K, inflated by repeating each entry h times. The same holds for the border columns (Kc).
3. **K is a Hadamard matrix of order 4c with 4 rows removed.** The rows of K are mutually orthogonal. The 4 missing rows
   are exactly **H4 ⊗ 1_c**: a 4 × 4 Hadamard matrix acting on the four core blocks, constant across each block's
   c cosets. This shows up after removing one real-character twist, which is the sign ambiguity in fixing a
   group-developed block (f and χ·f give equivalent matrices).
   - Equivalently, every border row is **balanced inside each core block**: equal numbers of + and − cosets.
   - This forces c to be even, so b is odd, which is the other session's theorem, seen here directly.
   - Kc has the same structure and spans the same space as K.
   - **K is not a Kronecker product H4 ⊗ L in any of the three.** No border row has rank 1 as a 4 × c array.
     For c = 6 it cannot be one: two ±1 vectors of length 6 are never orthogonal.
4. **Core.** A Goethals–Seidel quadruple over G with Σ_k PAF_k = −4b on H\{0} and +4 off H. The core's coset sums,
   i.e. its compression onto G/H, are "constant per block plus a small deviation". The deviation is forced: the total
   power at the nontrivial characters of G/H is exactly 4b.

So each matrix is literally an **inflated Hadamard matrix of order 4c**. Take W of order 4(b+1) that contains the rows
H4 ⊗ 1_c. Keep its other 4b rows as the border, repeating each entry h times. Replace the H4 ⊗ 1_c rows by a GS core
over G, whose compression mod H plays the role of those rows. Put an H(4b) in the corner. The b "extra" units of
PAF deficit on H are exactly what the border rows absorb.

**Similarity verdict.** At the level of architecture, the three examples are the same construction: identical layering,
identical border law, identical core condition. That is strong evidence that one method or program produced all three.

## 2. Why those borders

Given v, the border is almost forced once c is chosen.
- It needs an H(4b) for the corner.
- It needs 4b rows of an H(4c) that are orthogonal to H4 ⊗ 1_c.
- Both exist for every odd b in range. So the border itself is not what limits the choice.

The real choice is c (equivalently b) together with G ⊃ H. The admissible options are the even c ≥ 4 with c | v+1
(so h = (v+1)/c − 1); see `logs/bge3_border_options.log`.

| v | options c: h | chosen |
|---|---|---|
| 479 | **4: 119 = 7·17**, 6: 79 (prime), 8: 59 (prime), 10: 47 (prime), 12: 39 = 3·13, 16: 29 (prime) | c = 4 |
| 347 | 4: 86 = 2·43, **6: 57 = 3·19**, 12: 28 = 4·7 | c = 6 |
| 359 | 4: 89 (prime), 6: 59 (prime), **8: 44 = 4·11**, 10: 35 = 5·7, 12: 29 (prime) | c = 8 |

**All three choices fit one rule: take the smallest c ≥ 4 for which h is neither prime nor twice a prime.** In each
case the next composite option (c = 12, 12 and 10) was passed over, so the choice looks like "smallest b", not
"best-looking H".

Plausible reason [E]: H is cyclic in all three, and a composite h lets a symmetry group act with orbit sizes that
match the lag classes. The observed stabilisers are:
- 479: order 12, a twisted GL(2,2) × multiplier group, 46 affine orbits;
- 347: order 6 (scalar), 70 orbits;
- 359: order 5, acting on Z11 only, 96 orbits.

The rule is fitted after the fact to 3 data points, so treat it as a hypothesis, not a finding.

Applied to the open orders (details in the log):

| v | rule picks c (b), h | remark |
|---|---|---|
| 515 | 4 (3), 128 = 2⁷ | |
| 523 | 4 (3), 130 = 2·5·13 | |
| 537 | **none**: v+1 = 2·269 has no even divisor ≥ 4 | **no b ≥ 3 border exists at all** |
| 571 | 22 (21), 25 | c = 4 gives h = 142 = 2·71 |
| 573 | 14 (13), 40 | |
| 599 | 6 (5), 99 = 9·11 | |
| 643 | 4 (3), 160 | |
| 647 | 4 (3), 161 = 7·23 | same shape as 479 (Z2² × Z_{7·p}) |
| 659 | 4 (3), 164 = 4·41 | |
| 669 | 10 (9), 66 | |
| 719 | 6 (5), **119** | **same H = Z119 as the solved v=479** |
| 721 | 38 (37), 18 | |
| 739 | 4 (3), 184 = 8·23 | |

## 3. Where the three differ (method variants within one architecture)

- **Symmetry used.** 479 has a twisted automorphism of order 12 that permutes the three nontrivial Z2²-cosets together
  with a multiplier. 347 has only the scalar multipliers (order 6). For 359 I searched the full Aut(G), all
  384 × 10 = 3840 automorphisms, and it has only the scalar order-5 multipliers on Z11; there is no hidden twist.
  So the amount of symmetry imposed differs a lot: 46, 70 and 96 orbits.
- **Shape of the compression.** In 359, all four blocks deviate from constant on **the same single coset** of G/H:
  A = (−2⁷, 0), B = (−2⁷, −4), C = (−4⁷, −2), D = (0⁷, 4). In 479 the deviating coset differs from block to block,
  and one block (D) is exactly constant. In 347 one block (A) is constant and the others deviate on 2 or 3 cosets.
- **An extra exact relation in 359.** D(1,0,e,f) = (−1)^e · D(0,1,e,f): two Z2²-slices of block D are copies of each
  other twisted by the order-2 character of Z8. No random block with D's own symmetry shows this. Nothing comparable
  appears in 479 (0 slice relations beyond its symmetry) or in 347 (one exact B–C slice match; 0 in 10 null samples).
- **The extension type differs.** 359's G/H is a non-split extension, because H contains the index-2 subgroup of Z8.
  479 and 347 are split, with G = (G/H) × H.
- **479 in Z2²-components.** For each character χ of Z2², the four blocks' χ-components on Z119 are complementary:
  their power spectra sum to 4v at every ω ≠ 0. So a 479 core is four interlocked complementary quadruples over Z119,
  with entries in {0, ±2, ±4}. No component is shared between blocks. Blocks A and B nevertheless agree on 83% of
  positions under the best affine map, against a null maximum of 76%. That is a near-copy, not an exact shared piece.

**Reading.** One architecture throughout (layers 1–4), with 479 and 347 closest to each other (split G, cyclic H,
per-block compression deviations). 359 is the variant: a non-split extension, weak symmetry, one special coset shared
by all blocks, and an exact twisted slice relation. That pattern fits a later or differently seeded version of the
same program, perhaps one that fixes a "special coset" first. That last part is interpretation [E].
