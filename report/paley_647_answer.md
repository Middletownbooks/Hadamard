# H(2588) (v = 647) in the 647 ×353 doubly-staged slice: exact cross-level solver, calibration, cost and run

Date: 2026-10-14. Branch `claude/focused-planck-j1645f`. Answers the 647 ×353 prompt (validate and calibrate the cross
level, cost it, run it, find a same-shape control).
Scripts: `theory/p647_*.py`, `scripts/verify_h2588.py`. Logs: `logs/p647/`.

**Tags.** [M] exact computation (script named); [P] proved here (proof in the text); [E] Gaussian estimate (estimator
and calibration named); [L] literature recollection, unchecked.

---

## 0. Bottom line and run status

*Status: in progress (this section is updated as results arrive).*

- **(A)** No H(2588) yet.
- **(B)** Calibration runs on 4 joint solutions of y_B #0 are under way (§2).
- **(C)** No obstruction found; every check so far is consistent with a non-empty cross level.

---

## 1. Tools and their validation [P, M]

### 1.1 Exact identities used by the solver [P]

Per block k and Z₂-layer e the labels of the 35 M-orbits of the layer form a 5 × 7 ±1 matrix with prescribed weighted
row sums (y_B) and column sums (y_A); `nd_cross.Cross.layer_labellings` lists them exactly (lists L0, L1).
With a₀, a₁ the two layers at the 35 orbit representatives of Z₃₂₃, σ = (a₀ + a₁)/2 and δ = (a₀ − a₁)/2.

- **P1 (additive σ code).** For a ∈ {±1}, (a₀ + a₁)/2 + 1 = (a₀ + 1)/2 + (a₁ + 1)/2. Hence
  code(σ) = Σ_b (σ_b + 1)·3^b = c₀(a₀) + c₁(a₁) with c_e(a) = Σ_b ((a_b + 1)/2)·3^b, and code is a bijection from
  {−1, 0, 1}³⁵ onto [0, 3³⁵) ⊂ [0, 2⁶³). So the distinct σ of a block are the distinct values of c₀[i] + c₁[j]: one
  int64 `unique` over n₀n₁ ≤ 2.1·10⁶ values, no row sorting.
- **P2 (D-level by even lags).** a₀ = σ + δ and a₁ = σ − δ give PAF(a₀) + PAF(a₁) = 2(PAF(σ) + PAF(δ)). So given
  an S-solution (Σ_k PAF(σ_k) = T_S), the D-equations Σ_k PAF(δ_k) = T_D are equivalent to
  Σ_k (PAF(a₀,k) + PAF(a₁,k)) = T_E = 2(T_S + T_D) = (2584, −4, …, −4), the even-lag PAF of x on Z₆₄₆. The D-check
  needs only the per-layer PAFs (n₀ + n₁ vectors per block), no per-pair D-signature.
- **P3 (block energy is fixed mod 24).** Write σ as a 5 × 7 matrix (rows {0, A₁..A₄}, columns {0, B₁..B₆}); its
  weighted margins are s_B(i) = σ_{i0} + 3Σ_{j≥1} σ_{ij} and s_A(j) = σ_{0j} + 4Σ_{i≥1} σ_{ij}, fixed by the joint
  solution. Since σ_{ij} ∈ {−1, 0, 1} and these are distinct mod 3 and mod 4, σ_{i0} ≡ s_B(i) (mod 3) and
  σ_{0j} ≡ s_A(j) (mod 4) fix the whole first row and column. Since σ ≡ σ² (mod 2), the number n₁₂ of nonzero cross
  entries satisfies n₁₂ ≡ Σ_{i≥1} (s_B(i) − σ_{i0})/3 (mod 2). With orbit sizes 1, 3, 4, 12,
  ‖σ‖² = σ₀₀² + 3n₃ + 4n₄ + 12n₁₂ is therefore fixed mod 24 per block. [M: every list of every calibration joint
  solution lies in one class mod 24; the S-solutions sit in left-energy classes 24 apart.] The join's energy split
  already respects this, so it costs nothing extra and saves nothing extra.

### 1.2 The exact S-join (`theory/p647_join.py`) [M]

- Four lists of distinct σ with 20-coordinate signatures (PAF at the 20 ⟨M, −1⟩-lag classes of Z₃₂₃; coordinate 0
  is the energy ‖σ‖²). Pairing (a, b | c, d).
- Every row gets a 64-bit additive hash and a residue mod a prime P = 4093 (both linear, so additive). For every left
  energy EL and every group of residues, the left and right pairs of that class are emitted (inner lists sorted by
  (energy, residue) with a CSR index), scattered into 2^rbits partitions by the top key bits, and each partition is
  joined in cache. The exact pair count of every (EL, residue) class is precomputed by a 2-D FFT, so passes have a
  known size. Packed 8-byte entries (24-bit fingerprint, 2 × 20-bit indices); a fingerprint match is confirmed on the
  64-bit key, a key match on the exact integer signature.
- **Self-test** (`python3 theory/p647_join.py --selftest`): 40/40 (16-byte kernel) and 40/40 (packed kernel) against a
  brute-force 2+2 dictionary join on random instances, and `nd_ejoin.ejoin4` agrees 40/40.
- **Speed** [M]: 22–26 ns per emitted pair on one energy class of a real joint solution (EL = 319 of j = 9902,
  1.05·10¹⁰ pairs, 151,950 S-solutions, identical hit sets for five kernel settings); about 30 ns under a 4-worker load.
  The per-x overhead and the scatter dominate; the 523 kernel ran at a similar 17–26 ns.
- The old `nd_ejoin.ejoin4` re-enumerates every pair in each of its up to 64 partition passes, which is why it could
  not run the cross level at this size.

### 1.3 The cross pipeline (`theory/p647_cross.py`) [M]

Lists → distinct σ (P1) → S-join → D-check of every S-solution (P2; a compiled 2+2 on 64-bit hashes over the product of
its σ-preimages, every match re-checked exactly) → x_k on Z₆₄₆ → full cyclic PAF check.

**Validation.**
- Plants on real joint solutions (one random pair per block; targets = its own S and even-lag sums; the join run on
  the plant's left-energy class): 2/2 so far, S-level and full level (j = 9902: 8,252 S-solutions in the class,
  2.3·10⁷ D-candidates, 1 D-hit = the plant; j = 6526: 1,292, 575,903, 1).
- Reduced instances against a brute force that never uses the pipeline's signatures (every kept pair → its sequence on
  Z₆₄₆ → its full 646-lag PAF by FFT; a 2+2 dictionary on the exact 646-vector): 3/3 agree (`--brute`).
