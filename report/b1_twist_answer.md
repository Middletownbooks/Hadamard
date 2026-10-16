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
   - **Cyclic G = Z_{v−1}:** 521 slices, 195 'BC' and 107 'BCD' plus the untwisted ones for reference; all
     19 orders; 8–150 orbits.
     - The exact quotient screens kill 381 slices, 229 of them at the parity quotient Z₂.
     - Surviving twisted slices with log₂E > 0: the cheapest is 571, u = 103 (d = 12), 'BC': 66 orbits, 37 classes,
       log₂E = 42.5, **n = 198**.
   - **Noncyclic G** (Z₂² × Z₁₆₇, Z₂² × Z₁₉₇, Z₂³ × Z₆₇, Z₂ × Z₃² × Z_q, with GL(2,2), GL(3,2) and GL(2,3)
     automorphisms): FILL.

5. **(D)** No search was run: nothing is within budget (§5).

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

**Check against the screen [M].** Every 'BCD' slice at 599, 659 and 907 dies at Z₂ in `scripts/b1tw_screen.py`.
- 599 and 659 are ruled out by the lemma.
- 907's slices die by their letters.
- 571, 739 and 859 keep survivors.

**For 4 | m** (537, 669, 789) with G ⊇ Z₂², the three order-2 characters together give
v = a² + 4(b₁² + b₂² + b₃²) with a odd, under a GL(2,2) 3-cycle. That is always solvable, so the twist is compatible
there. This is the 479 mechanism.

## 3. Stage A: estimator and validation

FILL

## 4. Stage B: the screens

FILL

## 5. Stage C/D: costs of the best survivors

FILL

## Reproduce

FILL
