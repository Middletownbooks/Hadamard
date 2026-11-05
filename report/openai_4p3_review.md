# Review of the OpenAI model's n = 4p + 3 answer: normal form, universal k4 border, elliptic seeds

*Claude Code review, 2026-11-05, of `report/openai_4p3_answer.md` (branch `codex/order-3452-compression`, commit
27cfc29), as asked by `report/prompt_4p3_review.md`.*

**Credit.**
- The theorems, the proofs, the elliptic-character ansatz, the certificates `logs/openai_4p3/border_certificates.json`
  and the scripts `theory/openai_4p3_*.py` are by the OpenAI model.
- The bordered construction is Kumar's form of Theorem 7 of Egan–Kotsireas–London–Ó Catháin (EKLO).
- The frame-search border reduction that produced the certificates is the repository's (`theory/bt_border.py`).
- The real cores are:
  - the 479 core of Alpöge, Voinov and Reynolds-Haertle's published H(1916), as extracted in this repository
    (`data/core_v479_b3.txt`; `forensics/originals/matrix_H1916_v479_w_b3.txt.gz`);
  - Kumar's 3964 core (`hadamard/3964/certificate.json` in https://github.com/vatskishan03/hadamard, Apache-2.0).
- The review scripts `theory/4p3_review_*.py` are new.

**Independence.** No review script imports anything from `theory/openai_4p3_*.py`. The matrix checks
(`4p3_review_border.py`, `4p3_review_hverify.py`) import no repository module at all. The GS array, group,
normalisation, assembly, seed and join code are written from their definitions. The verifier runs as a separate process and
shares no code with the assembly. Two scripts use repository code, both outside any check of OpenAI code:
`4p3_review_qxcheck.py` reads the census Q-level lists (those lists are the reference data being compared), and the
estimates call the repository's lattice Gaussian `theory/q523_gauss.gauss_bits`, as the prompt asks.

**Tags.** [P] proved, [M] exact computation (script named), [E] estimate. A cap or a timeout is never a zero.

## Bottom line

| claim | verdict | evidence |
|---|---|---|
| **1. Compression normal form** [P] | **holds** | The proof is correct for both Q. It has no gap; three terse steps are written out in §1.1. An exhaustive brute force for every odd m ≤ 199 on both Q finds every solution in normal form: the ordered counts equal the normal-form prediction 4096·R(m) exactly [M]. All 2,742 census Q-level solutions in 11 structures are in normal form [M]. |
| **1′. Small-repair obstruction** [P, M] | **holds** | The minimum-change table is reproduced exactly at all 16 targets, both columns [M]. The per-index gate (e ∈ {2, 4, 6, 8, 12, 16}) gives the same caps [M]. The fixed-budget inequalities are correct [P] and agree with the exact gate for e ∈ {1, 2, 3, 4, 6, 8}, c ≤ 4 and every prime p ≤ 1200 [M]. |
| **2. Universal border** [P, M] | **holds; recorded as a repository result (Theorem B, §2.1)** | Both certificates pass every identity under independent integer code [M]. Matrices of orders **1916** and **3964** were assembled from real cores with only the universal certificate and verified [M]. So were 6 randomly scrambled copies of those cores and 9 elliptic cores (orders 60–124). Every corrupted core or matrix was rejected. |
| **3. Elliptic seeds** [M] | **reproduced; ansatz closed [M + E]** | The five controls rebuild identically and border to Hadamard matrices [M]. The seed gate is exactly empty at p = 11, 17, 19, 23, 29, 31 (two implementations, no shared code) and at every prime 37–97 (one implementation) [M]. The calibrated expectation of seed-gate solutions is at most **2^−363** at every target prime, against the closing threshold 2^−10 [E]. The phase gate alone would also need ≈ 2^−1000 luck at p = 149 [E]. At p = 13, where the seed gate holds, an exhaustive phase search finds no core on either Q [M], so among primes up to 97 the ansatz has cores only at 3, 5 and 7. |

**What changes in the repository.**
1. **Theorem B (k4 border theorem, OpenAI model, verified here).** Every valid core of the k4 family with odd |K|
   completes to a Hadamard matrix, using one fixed certificate per Q. The border problem of that family is closed.
   - It covers the k4 hosts of 38 of the 85 open orders below 10000: those with n ≡ 7 (mod 8). These include all 16
     n = 4p + 3 targets and 3452.
   - A README line and a note in `report/b3_push_answer.md` §7 record it.
2. **Q-level in closed form.** The compression normal form solves the Q-level of every such host in closed form.
3. **Ansätze closed.** The elliptic ansatz is closed at the targets [M + E]. The small-repair cyclotomic classes
   (fixed e, c ≤ 3) are closed at all 16 targets, and fixed (e, c) cannot give an infinite family [P, M].
4. **A correction to the answer's table.** The elliptic ansatz also has cores at **p = 7 on ℤ₄**; the answer found
   none there. One was found with the first seed choice tried and bordered to H(124).
5. **No Hadamard matrix at an open order** is claimed, by the answer or by this review.

---

## 1. Claim 1: the compression normal form (Task 1)

### 1.1 The proof, line by line

**The statement checked.**
- Let m = |K| be odd, Q = ℤ₂² or ℤ₄, and z₀, …, z₃ odd vectors on Q with Σᵢ PAF_{zᵢ} = 4m + 12 at 0 and 4m off 0.
- Then, after permuting and negating whole rows, one row is constant (a, a, a, a) and three rows are (tⱼ + 2, tⱼ, tⱼ,
  tⱼ) up to the position of the exceptional entry, with m = a² + Σⱼ tⱼ(tⱼ + 1).
- Conversely, every such quadruple satisfies the Q-equations.

The answer writes p for m. The proof never uses primality, so it holds for every odd m.

**ℤ₂²: correct.**
1. *Integrality, left implicit.* Every Walsh coefficient ẑ(η) is a signed sum of four odd integers, hence even, so
   w = ẑ/2 is an integer.
2. *Column norms.* The Fourier form of the Q-equations gives Σᵢ ẑᵢ(η)² = 12 for each nontrivial η. So
   Σᵢ wᵢ,η² = 3, i.e. three entries ±1 and one 0 in each nontrivial column.
3. *Parity.* For η ≠ θ, wη − wθ = ±z(q₁) ± z(q₂), where q₁, q₂ are the two elements on which η and θ differ. The
   answer says "difference"; it is a signed sum, and even either way. So a row's three nontrivial w's have one parity.
4. *Constant row.* A zero in some column forces that row's other two nontrivial entries to be even, and they lie in
   {0, ±1}, so they are 0. Hence all three zeros lie in one row, and that row is constant.
5. *Singleton rows, left implicit.* In the other rows the nontrivial entries are ±1, so ẑ(η) = ±2. Inverse Walsh
   gives z = (ẑ₀ ± 2(4δ_s − 1))/4, where ẑ₀ is the principal coefficient, the sign is the product of the three
   nontrivial signs, and s is determined by them. That is three equal entries and one differing by ±2. Negate if
   needed so that it is +2.

**ℤ₄: correct.**
1. With A = (z₀ − z₂)/2, B = (z₁ − z₃)/2 and C = (z₀ − z₁ + z₂ − z₃)/2: |ẑ(χ₁)|² = 4(A² + B²) and ẑ(χ₂) = 2C.
   So Σ(A² + B²) = 3 and ΣC² = 3.
2. C − (A − B) = z₂ − z₃ is even.
3. Three rows have C odd. Each of them needs A² + B² ≥ 1, and the budget is 3, so each has exactly one of A, B equal
   to ±1 and the fourth row has A = B = C = 0.
4. *Left implicit: "solving these relations".* Write u for the repeated value.
   - (A, B) = (±1, 0) gives z₁ = z₃ = u and (z₀, z₂) = (u + C + A, u + C − A).
   - (0, ±1) gives z₀ = z₂ = u and (z₁, z₃) = (u − C + B, u − C − B).
   - In all four sign cases this is three entries u and one u ± 2.
   - The C = 0 row is constant.

**Energy and the converse.**
- A singleton row has energy 4t² + 4t + 4 and off-zero PAF 4t² + 4t on either Q. A constant row has 4a² at every
  lag.
- The zero-lag sum then gives m = a² + Σ t(t + 1), and the off-zero sum gives 4m. This proves both directions.

**Verdict [P].** The theorem holds as stated, for both Q and every odd m. Each pronic value k(k + 1) has exactly one
odd t, so each singleton class is ±(t·1 + 2δ_s): 8 rows. Each constant class is ±a·1: 2 rows. Hence the number of
ordered solutions is 4·2·8³·R(m) = 4096·R(m), where R(m) = #{(a, k₁, k₂, k₃) : a > 0 odd, kⱼ ≥ 0,
a² + Σ kⱼ(kⱼ + 1) = m}. The brute force tests exactly this count.

### 1.2 Brute force [M]

**Method** (`theory/4p3_review_nf.py`, logs `logs/4p3_review/nf_{v4,z4}.{json,log}`; no repository imports).
1. List every odd 4-vector of energy ≤ 4m: 196,912 rows at m = 199.
2. Group the rows by PAF vector: 20,047 classes for ℤ₂² and 9,221 for ℤ₄. No class mixes shapes.
3. Join the classes 2 + 2 on exact linear integer keys, for each odd m. This gives:
   - N(m), the number of ordered solutions, by the weighted convolution Σ W(k) W(T − k);
   - every multiset of PAF classes that solves the target, each checked for normal form.
4. Second method, for m ≤ 9: a 3 + 1 count with dictionary lookups, no keys and no classes. It agrees.

**Result.**
- For every odd m from 3 to 199, on both Q: N(m) = 4096·R(m) exactly.
  - Zero multisets outside normal form.
  - The multiset count equals the normal-form count (up to 33 multisets per m).
- This covers m = 119 (the 1916 core's |K|) and the targets 149 and 179:
  - R(149) = 75: 307,200 ordered solutions, 14 multisets;
  - R(179) = 90: 368,640 ordered solutions, 17 multisets.
- Run times: 551 s for ℤ₂² and 136 s for ℤ₄.

The counts per m, as m:R(m); the ordered count is N(m) = 4096·R(m) on both Q.

```
3:3 5:3 7:4 9:7 11:6 13:9 15:13 17:9 19:10 21:15 23:15 25:13 27:19 29:18 31:16 33:30 35:21 37:19 39:27 41:21
43:31 45:31 47:24 49:25 51:39 53:33 55:28 57:48 59:30 61:35 63:54 65:33 67:34 69:52 71:42 73:45 75:51 77:39
79:45 81:55 83:51 85:50 87:70 89:45 91:46 93:78 95:48 97:54 99:80 101:57 103:63 105:78 107:54 109:55 111:75
113:84 115:58 117:79 119:60 121:61 123:117 125:63 127:74 129:87 131:72 133:81 135:91 137:75 139:77 141:121
143:93 145:81 147:99 149:75 151:76 153:126 155:90 157:79 159:117 161:81 163:99 165:126 167:93 169:98 171:115
173:105 175:95 177:130 179:90 181:91 183:171 185:93 187:94 189:144 191:105 193:124 195:150 197:114 199:108
```

### 1.3 Cross-check against the repository's census Q-level lists [M]

`theory/4p3_review_qxcheck.py`, logs `logs/4p3_review/qxcheck.{json,log}`.
- The lists come from `theory/bt_gates.Gates.q_level` on the census structures: an exact join that knows nothing
  about the normal form.
- Every listed solution was expanded to full Q-vectors and checked by new code against the structure's own Q
  addition table. Every one passes the quotient correlations and is in normal form.

| structure | G | \|K\| | Q-level solutions | in normal form | distinct (a; t) |
|---|---|---:|---:|---:|---:|
| 3964 control (Kumar's structure) | ℤ₂²×ℤ₂₄₇ | 247 | 2 | 2 | 2 |
| 1916 split structure | ℤ₂²×ℤ₁₁₉ | 119 | 2,112 | 2,112 | 12 |
| 863 S7 host 2e93be9718324148 | ℤ₂²×ℤ₅×ℤ₄₃ | 215 | 576 | 576 | 9 |
| 863 ℤ₂² S168 / S84 (0f46…, 4b53…) | ℤ₂²×ℤ₅×ℤ₄₃ | 215 | 1 + 1 | 2 | 1 |
| 863 ℤ₄ S56 / S28 (459f…, 9d13…) | ℤ₄×ℤ₅×ℤ₄₃ | 215 | 8 + 8 | 16 | 1 |
| 839 ℤ₂² S270 / S135 (3e80…, 8f2d…) | ℤ₂²×ℤ₁₁×ℤ₁₉ | 209 | 1 + 1 | 2 | 1 |
| 839 ℤ₄ S90 / S60 (1a81…, 5131…) | ℤ₄×ℤ₁₁×ℤ₁₉ | 209 | 8 + 24 | 32 | 1 + 3 |
| **total (11 structures)** | | | **2,742** | **2,742** | |

### 1.4 Claim 1′: the minimum-change table and the fixed-budget theorem [M, P]

`theory/4p3_review_repair.py`, logs `logs/4p3_review/repair_table.{json,log}`.

**The table.**
- I enumerated every normal form with new code and recomputed both columns of the answer's §4 table at all 16
  targets. They agree exactly (29/4, 29/4, 37/5, …, 61/7), and also with its log `compression_host.json` at the 23
  primes 5–97.
- Per index e ∈ {2, 4, 6, 8, 12, 16}, the least common budget c for which a normal form has all 16 cell sums
  attainable equals the balanced cap at every target.
- The attainable-sum formula S(p, e, c) equals the set of sums of explicit cells at p ∈ {5, 7, 13, 17, 19, 29, 31,
  37}, e ∈ {1, 2, 3, 4, 6} and c ≤ 2. The explicit cells use a primitive root, every class-sign vector, both zero
  signs and every flip set of size ≤ c.

**The fixed-budget theorem [P].**
Write s = (p − 1)/e.
- *e even.* The bases are η + 2sk, |k| ≤ e/2. When 2s − B₀ > √(4p + 12), only k = 0 survives the per-cell energy
  bound z² ≤ 4p + 12. So |z| ≤ B₀ = 2c + 1, hence 4p + 12 ≤ 16B₀², and there is no core when p > 4B₀² − 3.
- *e odd.* Every |b| ≥ s − 1. With L = s − 1 − 2c > 0 every |z| ≥ L, so 16L² > 4p + 12 excludes everything.
- Both conditions hold for all large p at fixed (e, c).

**Numeric check.** For e ∈ {1, 2, 3, 4, 6, 8}, c ≤ 4 and every prime p ≤ 1200, the exact gate never passes where an
inequality says "no core" (0 contradictions) [M]. The passing primes stop early. For example, at (e, c) = (2, 3) the
largest is p = 163, against the bound 4B₀² − 3 = 193; at (2, 4) it is 281, against 321. The largest over all these
(e, c) is 281.

**Scope (as the answer says).** These are compression-level obstructions: zeros for the stated classes. A pass is
not evidence for cores.

---

## 2. Claim 2: the universal border (Task 2)

### 2.1 Theorem B (k4 border theorem; OpenAI model, verified here) [P, M]

**Statement.**
- Let G be abelian of order 4m with m odd, K ≤ G of order m, and Q = G/K ∈ {ℤ₂², ℤ₄}.
- Let x₀, …, x₃ : G → {±1} with Σᵢ PAF_{xᵢ} = 4|G| at 0, −12 on K ∖ 0 and +4 off K.
- Then permuting the rows, negating rows and translating rows by elements of G yields a core whose Goethals–Seidel
  array (Kumar's layout) is completed to a Hadamard matrix of order 4(|G| + 3). The border is the fixed certificate
  (X, Rh, Ch) for Q in `logs/openai_4p3/border_certificates.json`, lifted constantly along K.

**Proof.**
1. *Splitting.* m is odd and |Q| = 4, so G = Q′ × K with Q′ ≅ Q. Translating xᵢ by an element of Q′ translates its
   compression zᵢ(q) = Σ_{k∈K} xᵢ(q + k) on Q. Permutation, negation and translation preserve Σ PAF.
2. *Compressions.* The zᵢ are odd and satisfy the Q-equations: sum the defining equations over each coset.
3. *Normal form.* By Claim 1 (§1), after these operations z₀ = c₀·1 and zⱼ = cⱼ·1 + 2δ₀ for j = 1, 2, 3.
4. *Transfer.* The GS array is linear in its four sequences. So Mc(z) − Mc(z*), with z* = (0, 2δ₀, 2δ₀, 2δ₀), is
   the array of four constant sequences, and every 4 × 4 block of it is constant. Every 4-entry block of every row of
   Rh sums to 0 [M]. Hence Rh Mc(z)ᵀ = Rh Mc(z*)ᵀ, and the certificate identity X Chᵀ + Rh Mc(z*)ᵀ = 0 [M] holds for
   z.
5. *Assembly* (Kumar's reduction, `theory/bt_border.py`). Lift R[x, (i, g)] = Rh[x, 4i + q(g)] and
   C[(i, g), x] = Ch[4i + q(g), x].
   - Summing a GS block over a K-coset of columns gives the GS block of the compressions. So R Mᵀ and X Cᵀ are the
     lifts of Rh Mcᵀ and X Chᵀ, and the off-diagonal block of H Hᵀ vanishes.
   - The top-left block is X Xᵀ + |K| Rh Rhᵀ = (12 + 16m) I = 4n I.
   - The bottom-right block is C Cᵀ + M Mᵀ. Here C Cᵀ is the lift of I₄ ⊗ Γ_Q (12 within a coset, −4 across). Over
     an abelian group the GS blocks commute as required, so M Mᵀ = I₄ ⊗ [Σᵢ PAFᵢ(h − g)]: 4|G| on the diagonal,
     −12 on K ∖ 0 and +4 off K. The sum is (4|G| + 12) I. ∎

Primality and cyclicity of K are never used.

**Scope.**
- Theorem B covers every n ≡ 7 (mod 8) (m = (n − 3)/4 odd). Below 10000 that is 38 of the 85 open orders:
  - the 16 n = 4p + 3 targets: 2396, 2876, 4124, 4348, 5308, 5756, 5884, 6236, 7324, 7388, 7484, 8156, 8348, 8380,
    9148, 9596;
  - 3004, 3452, 4156, 4412, 4892, 5276, 5468, 5692, 5948, 6172, 6268, 6332, 6652, 6876, 7132, 7516, 8252, 8828,
    8860, 9308, 9532, 9692.
- It covers both ℤ₂² and ℤ₄ hosts.

### 2.2 The certificates [M]

`theory/4p3_review_border.py`, log `logs/4p3_review/border.json`. Mc was built from the GS-array definition over Q.

| check | ℤ₂² (frame 75) | ℤ₄ (frame 75) |
|---|---|---|
| all entries ±1 (X 12×12, Rh 12×16, Ch 16×12) | yes | yes |
| X Xᵀ = 12 I | yes | yes |
| Rh Rhᵀ = 16 I | yes | yes |
| Ch Chᵀ = I₄ ⊗ Γ_Q (12 diagonal, −4 off) | yes | yes |
| X Chᵀ + Rh Mc(z*)ᵀ = 0 | yes | yes |
| each 4-entry block of each Rh row sums to 0 | yes | yes |
| Rh Mc(1ᵢ)ᵀ = 0 for each row i (transfer, direct) | yes ×4 | yes ×4 |
| own Mc(z*) equals the Mc stored in the file | yes | yes |

### 2.3 Real cores, assembled with only the universal certificate [M]

**Pipeline**, per core:
1. Check the core by integer PAF on every lag.
2. Normalise: row permutation, signs, and Q-translations by the element of the coset with zero odd part.
3. Check the normalised core again.
4. Assemble H.
5. Verify H Hᵀ = 4n I in a separate process (`theory/4p3_review_hverify.py`) by two exact methods: a bit-packed
   popcount Gram over all pairs, and a float64 Gram, which is exact here since every partial sum is an integer below
   2⁵³.

**Negative controls**, per core: one sign of the core flipped (the core check must reject it), one entry of H flipped
(the verifier must exit non-zero), and a seeded random flip inside the verifier.

| core | G | \|K\| | Q | order | compressions as given | normalisation (row, sign, q-shift) | H sha256 (int8, first 16) | result |
|---|---|---:|---|---:|---|---|---|---|
| **1916** (Alpöge et al., `data/core_v479_b3.txt`) | ℤ₂²×ℤ₁₁₉ | 119 | ℤ₂² | 1916 | (3,1,1,1) (1,1,1,3) (9,11,9,9) (5,5,5,5) | (3,+,0) (0,+,0) (1,+,3) (2,+,1) | a9c9abcdc57ab5d4 | **Hadamard** |
| **3964** (Kumar) | ℤ₂²×ℤ₂₄₇ | 247 | ℤ₂² | 3964 | (−11,−13,−13,−13) (−5,−7,−7,−7) (1,−1,−1,−1) (7,7,7,7) | (3,+,0) (0,+,0) (1,+,0) (2,+,0) | 8fedf214123dac1f | **Hadamard** |
| 1916, 3 scrambled copies | ℤ₂²×ℤ₁₁₉ | 119 | ℤ₂² | 1916 | random (e.g. (−5,−5,−5,−5) (−1,−3,−1,−1) …) | negations and shifts used | 8ff5…, e403…, 9f2e… | Hadamard ×3 |
| 3964, 3 scrambled copies | ℤ₂²×ℤ₂₄₇ | 247 | ℤ₂² | 3964 | random (e.g. (13,11,13,13) …) | negations and shifts used | 06e7…, 0736…, 09a1… | Hadamard ×3 |
| answer's examples, rebuilt by new code | ℤ₂²/ℤ₄ × ℤ₃, ℤ₅, ℤ₇ | 3, 5, 7 | both | 60, 60, 92, 92, 124 | already normal | identity | de44…, ec93…, a3f5…, e637…, b6ea… | Hadamard ×5 |
| random elliptic cores (§3.1) | ℤ₂²×ℤ₅, ℤ₄×ℤ₅, ℤ₂²×ℤ₇, ℤ₄×ℤ₇ | 5, 7 | both | 92, 92, 124, 124 | already normal | identity | db1b…, 1737…, 52ce…, 991d… | Hadamard ×4 |

**About the table.**
- The scrambled copies apply a random row permutation, random row signs, random translations by arbitrary elements of
  G and a random automorphism of G (GL(2, 2) on ℤ₂², a unit on the odd part). This exercises every branch of the
  normalisation; the unscrambled cores need no negation.
- All 17 corrupted cores and 17 corrupted matrices were rejected.
- Total run time: 24 s.

**Pass criterion met.** 1916 and 3964 both verify, so the border problem of the k4 family with odd |K| is closed, and
the README line and the b3 §7 note were added. The matrices of orders 1916 and 3964 are known orders (Alpöge et al.
and Kumar). They are regenerable from the committed script and inputs, so only their hashes are logged.

**Not tested end to end: a non-cyclic K.** The proof does not depend on the structure of K. No core with K = ℤ₃² was
found to test it on, by either of two searches (`theory/4p3_review_smallcores.py`; logs
`logs/4p3_review/smallcores_*`):
- simulated annealing over ℤ₂² × ℤ₃² and ℤ₄ × ℤ₃², 10 minutes and 29 restarts each;
- exhaustive searches of seven slices invariant under an order-4 automorphism, 11–12 orbits per block, all empty
  [M for those slices].

The same annealer finds a core over ℤ₂² × ℤ₅ in 0.3 s. So this case is untested, not failed.

---

## 3. Claim 3: the elliptic seeds (Task 3)

`theory/4p3_review_elliptic.py` re-implements the ansatz from the text of the answer's Proposal 2. It imports no
repository code; the estimates load the repository's Gaussian.

### 3.1 Reproduction [M]

- **The five controls.** From the parameter objects stored in `logs/openai_4p3/core_p*_*.json` (data only), the new
  seed definition rebuilds sequences identical to the stored ones at (p, Q) = (3, ℤ₂²), (3, ℤ₄), (5, ℤ₂²), (5, ℤ₄) and
  (7, ℤ₂²).
  - All pass the integer core check, and all one-sign corruptions are rejected.
  - Bordered by Theorem B, they give Hadamard matrices of orders 60, 60, 92, 92 and 124 (§2.3).
  - Log: `logs/4p3_review/elliptic_examples.json`.
- **Random cores** (`cores --mode random`, pick seed 2026). The pick is a seeded uniform choice of a gate solution and
  a seed realisation, then a uniform choice among all its shift solutions; the answer's examples are excluded. Cores
  were found at p = 5 (both Q) and p = 7 (both Q), each on the first realisation tried.
  - The realisation at p = 7 on **ℤ₄** has cores. The answer reported only that its chosen seeds fail there.
  - Logs: `logs/4p3_review/elliptic_random_p{5,7}_{v4,z4}.json`.
- **The 16 compression witnesses (C).** The answer's witness tuples hold under new trace code, and their parameters
  rebuild to seeds with exactly those sums (`logs/4p3_review/elliptic_coverage_check.{json,log}`).

### 3.2 The exact seed gate (K) [M]

**Method.**
- Every parameter tuple (λ, ε, mask, u) is enumerated, together with every negative root for the switched rows.
- The distinct vectors P₀ and D = (3P + P′)/4 are kept, over t = 1, …, (p − 1)/2.
- The join is an exact 2 + 2 join: pairs (Dᵢ, Dⱼ), i ≤ j, on a linear 64-bit hash, which is exact under wraparound,
  so nothing is missed. Every hit is verified on the full integer vector.
- All solutions are counted, not just the first.

**The answer's run.** `python3 theory/openai_4p3_search.py --primes 23,29,31 --phase-max-p 0 --seconds 540 --out
logs/openai_4p3/next_seeds.json` (log `logs/4p3_review/openai_next_seeds.log`) returned `no_seed_quadruple` at 23,
29 and 31, uncapped. The new implementation confirms each zero.

| p | ordinary / switched vectors | gate solutions: distinct multisets; ordered parameter tuples | answer's script |
|---:|---|---|---|
| 3 | 2 / 3 | 1; 497,664 | pass (2/3) |
| 5 | 4 / 10 | 4; 22,394,880 (2^24.42) | pass (4/10) |
| 7 | 9 / 35 | 15; 27,869,184 (2^24.73) | pass (9/35) |
| 11 | 45 / 175 | **0** | no quadruple (45/175) |
| 13 | 87 / 264 | 12; 2,985,984 (2^21.51) | pass (87/264) |
| 17 | 156 / 496 | **0** | no quadruple (156/496) |
| 19 | 204 / 612 | **0** | no quadruple (204/612) |
| 23 | 308 / 924 | **0** | no quadruple (308/924) |
| 29 | 532 / 1,540 | **0** | no quadruple (532/1540) |
| 31 | 580 / 1,740 | **0** | no quadruple (580/1740) |
| 37 … 97 (14 primes) | 876 / 2,556 … 6,176 / 18,336 | **0 at every one** | not run |

**Reading the table.**
- The two implementations share no code. They agree on every vector count and every zero for p ≤ 31.
- Below 100, p = 13 is the only prime above 7 where the gate is nonempty. No positive appeared at 11, 17 or 19, so
  the answer's exclusion stands.
- Logs: `logs/4p3_review/elliptic_gate_small.json` (3–7), `elliptic_gate.{json,log}` (11–31) and
  `elliptic_gate_ext.{json,log}` (37–97, 222 s at p = 97).

### 3.3 Scaling of the seed gate [E]

**Model.**
- N_par is the number of ordered parameter quadruples: (16(p−1)(p−2)) · (24(p−1)(p−2))³. That is λ, ε, u and mask
  per row, and the negative-root choice in rows 1–3.
- "bits" is −log₂ of the lattice-Gaussian probability that the (K) vector hits −3·1. It is
  `theory/q523_gauss.gauss_bits` read at the true target, on 4,000 sampled parameter quadruples, uniform over
  (seed, root) pairs. The rank is always full, and the target lies at Mahalanobis distance 1.2–1.6.
- E_par = N_par · 2^−bits.
- E_dist = E_par / 2^17.92 estimates distinct-vector multisets. 2^17.92 = 6 · 24 · 12³: 6 orders of rows 1–3, 24
  parameter tuples per ordinary vector and 12 per switched vector.
  - It is the exact ratio at p = 13.
  - The per-vector ratios are exactly 24 and 12 at every p ≡ 3 (mod 4) from 19 to 83, and 22.7–23.6 and 11.6–11.9 at
    p ≡ 1 (mod 4).

| p | equations | log₂ N_par | bits | log₂ E_par | exact ordered count | log₂ E_dist | exact distinct |
|---:|---:|---:|---:|---:|---|---:|---:|
| 5 | 2 | 32.1 | 8.3 | 23.8 | 2^24.42 | 5.8 | 4 |
| 7 | 3 | 37.4 | 12.2 | 25.2 | 2^24.73 | 7.3 | 15 |
| 11 | 5 | 43.7 | 20.5 | 23.2 | **0** | 5.3 | 0 |
| 13 | 6 | 45.9 | 24.1 | 21.8 | 2^21.51 | 3.9 | 12 |
| 17 | 8 | 49.4 | 37.0 | 12.4 | 0 | −5.5 | 0 |
| 19 | 9 | 50.8 | 41.0 | 9.8 | 0 | −8.2 | 0 |
| 23 | 11 | 53.2 | 50.7 | 2.4 | 0 | −15.5 | 0 |
| 29 | 14 | 56.0 | 67.7 | −11.7 | 0 | −29.6 | 0 |
| 31 | 15 | 56.8 | 72.9 | −16.1 | 0 | −34.0 | 0 |
| 37–97 | 18–48 | 59.0–70.4 | 90.5–270.5 | −31.5 to −200.1 | 0 | −49.5 to −218.0 | 0 |

**Calibration** (exact counts at p ≤ 31, extended to 97).
- Where solutions exist (p = 5, 7, 13), the Gaussian matches the exact parameter count within 0.7 bits, in both
  directions: 0.6 bits low at p = 5, and 0.5 and 0.3 bits high at p = 7 and 13.
- At p = 11 it overestimates: about 2^5 distinct solutions predicted, none exist.
- Wherever it predicts fewer than one distinct solution (p ≥ 17), the exact gate is empty.
- No underestimate larger than 0.6 bits was seen.

**At the 16 targets** (`logs/4p3_review/elliptic_estimate_targets_{a,b}.json`):

| 4n | p | equations | log₂ N_par | bits | log₂ E_par | log₂ E_dist |
|---:|---:|---:|---:|---:|---:|---:|
| 2396 | 149 | 74 | 75.4 | 438.3 | **−362.9** | −380.8 |
| 2876 | 179 | 89 | 77.5 | 540.7 | −463.2 | −481.1 |
| 4124 | 257 | 128 | 81.7 | 810.6 | −728.9 | −746.8 |
| 4348 | 271 | 135 | 82.3 | 858.9 | −776.6 | −794.5 |
| 5308 | 331 | 165 | 84.7 | 1073.1 | −988.4 | −1006.3 |
| 5756 | 359 | 179 | 85.6 | 1174.4 | −1088.8 | −1106.7 |
| 5884 | 367 | 183 | 85.9 | 1202.8 | −1117.0 | −1134.9 |
| 6236 | 389 | 194 | 86.5 | 1282.3 | −1195.7 | −1213.7 |
| 7324 | 457 | 228 | 88.4 | 1534.8 | −1446.4 | −1464.3 |
| 7388 | 461 | 230 | 88.5 | 1547.7 | −1459.2 | −1477.1 |
| 7484 | 467 | 233 | 88.7 | 1571.1 | −1482.4 | −1500.4 |
| 8156 | 509 | 254 | 89.7 | 1728.6 | −1639.0 | −1656.9 |
| 8348 | 521 | 260 | 89.9 | 1772.1 | −1682.2 | −1700.1 |
| 8380 | 523 | 261 | 90.0 | 1780.0 | −1690.1 | −1708.0 |
| 9148 | 571 | 285 | 91.0 | 1961.7 | −1870.7 | −1888.7 |
| 9596 | 599 | 299 | 91.5 | 2067.5 | −1975.9 | −1993.9 |

**Why.**
- Each K-lag equation costs about 6 bits (5.9 at p = 149, 6.9 at p = 599), and there are (p − 1)/2 of them.
- The seed family has only about 8 log₂ p + 18 bits of freedom.
- So the expectation falls by about 3 to 3.5 bits per unit of p. It crosses one distinct solution at p ≈ 15 and is at
  most 2^−363 at every target.

### 3.4 The phase gate (S) [M, E]

The answer leaves this open. The estimate is a second, independent reason the ansatz cannot reach the targets.

**Method** (`phase-estimate`, log `logs/4p3_review/elliptic_phase_estimate.{json,log}`).
- Draw a random realisation that satisfies (C); a random seed choice is outside the lattice, because Σ_t of (S) is
  fixed by the compressions.
- Sample 3,000 uniform shift choices per row and apply the Gaussian at the target 4·1.
- This estimates the number of shift solutions per realisation, out of p¹² shift choices.
- Calibration: exact counts of the same realisations by the full 2 + 2 shift join, at p ≤ 13.

| p | log₂ p¹² | rank | bits | log₂ E[shift solutions] | exact |
|---:|---:|---:|---:|---:|---|
| 5 | 27.9 | 6 | 13.7, 13.5 | 14.2, 14.4 | 7,008; 16,128 (2^12.8, 2^14.0) |
| 7 | 33.7 | 9 | 21.2, 20.3 | 12.5, 13.4 | 3,120; 4,488 (2^11.6, 2^12.1) |
| 13 | 44.4 | 18 | 56.8, 53.3 | −12.4, −8.9 | 0; 0 |
| 17, 19, 23 | 49–54 | 24–33 | 78–123 | −29 to −69 | — |
| 29, 31, 37, 43, 53, 61 | 58–71 | 42–90 | 157–390 | −99 to −319 | — |
| 149 | 86.6 | 222 | 1108.2, 1110.7 | **−1021.5, −1024.1** | — |
| 179 | 89.8 | 267 | 1375.7, 1367.6 | −1285.9, −1277.8 | — |

**Reading the table.**
- At p = 5 and 7 the estimate is within 1.5 bits of the exact counts, and above them in all four samples.
- At the targets, a seed quadruple that passed (K) would still need about 2^−1000 luck in its shifts.

**p = 13.**
- The 12 gate solutions have 28,512 seed realisations (2^14.8), each with about 2^−9 to 2^−12 expected shift
  solutions. The answer's single failed quadruple is therefore expected, and the estimate predicts cores at p = 13
  only of order 2^4 in the whole ansatz.
- Global dilation k ↦ uk maps cores to cores, and rows 1–3 may be reordered. This leaves 2,376 orbits of
  realisations, exactly 28,512/12. The exhaustive search over them (`cores --mode exist-orbit`) is reported in §5.

### 3.5 Decision: **close the ansatz** [M + E]

The prompt's closing rule is met.
- The seed gate is exactly empty at p = 23, 29 and 31 [M], by two implementations. It is also empty at every prime
  from 37 to 97, at 11, 17 and 19, and nonempty only at 3, 5, 7 and 13. At 13 the phase gate is exactly empty on
  both Q (§5), so cores exist only at 3, 5 and 7.
- The calibrated expectation is at most 2^−363 at every target prime [E], about 350 bits below the threshold 2^−10.
- The model's largest underestimate at p ≤ 97 was 0.6 bits (p = 5). Its large error (p = 11) was an overestimate.
- The phase gate adds a further factor of about 2^−1000 per realisation at p = 149 [E].

No target prime is worth an exact seed gate for this ansatz.

**For reference: cost of an exact gate.**
- The pair table has about 2(p − 1)²(p − 2)² entries: 9.5·10⁸ at p = 149 (about 20 minutes with a chunked join) and
  2.6·10¹¹ at p = 599 (days) [E, from the measured 222 s at p = 97].
- At 2^−363 there is no reason to run it.

---

## 4. Next step for the n = 4p + 3 orders

**What these results settle.**
1. *Border.* Theorem B: a core alone now gives a Hadamard matrix at each of the 16 orders. Assembly and verification
   take seconds (`theory/4p3_review_border.py`, `theory/4p3_review_hverify.py`; or `theory/openai_4p3_border.py`).
2. *Q-level.* With K = ℤ_p prime, Q is the only exact compression level, and the normal form solves it in closed
   form. Every core has cell sums (a, a, a, a) and (tⱼ + 2, tⱼ, tⱼ, tⱼ) with p = a² + Σ tⱼ(tⱼ + 1). There are
   R(p) = 75 ordered forms at 149 and 90 at 179 (14 and 17 multisets). A search at these orders starts from these
   forms; no Q-level join is needed.

**What remains.**
1. *Seed-plus-shift ansätze fail.* Both proposed families are closed: the small-repair cyclotomic seeds (Claim 1′)
   and the elliptic seeds (§3). The scaling in §3.3–3.4 gives the reason.
   - An ansatz with polynomially many seeds per row and generic correlations pays about 6 bits per K-lag equation
     ((p − 1)/2 of them) and about 5 bits per off-K equation (about 3(p − 1)/2 of them).
   - Its freedom is O(log p) bits: about 8 log₂ p for the seeds and 12 log₂ p for the shifts.
   - Under this model, any such family falls below one expected solution by p ≈ 20, unless its correlations are far
     from generic.
2. *A viable proposal must be algebraic.* It needs an identity that makes the K-lag sums identically −12 and the off-K
   sums identically +4 (as Paley's quadratic character makes PAFs constant), leaving at most O(log p) equations for
   search.
3. *Screen new proposals cheaply.* The next request to a proposer should ask for that, and should screen any proposal
   at p ≤ 97 before any target run. The review's tooling does this in minutes:
   - the exact seed gate and its calibrated Gaussian (`4p3_review_elliptic.py gate / estimate`);
   - the phase-gate estimate (`phase-estimate`).
4. *Search.* The only search-based candidates with positive estimates at these orders remain the multiplier slices
   with r = 225–360 at 599 and 719 (`report/b3_push_answer.md` §0; estimates up to 2^76.7).
   - They have no exact level between Q and G. The normal form now fixes their Q-level, but adds no finer level.
   - By the calibration warnings (every exactly decided positive estimate at the open orders has been zero, up to
     2^85), they do not justify a long run either.
5. *Recommendation.* No compute on the n = 4p + 3 orders until an algebraic proposal passes the p ≤ 97 screens.
   Theorem B applies to every k4 host with odd |K|, so a core found by any method at any of the 38 orders with
   n ≡ 7 (mod 8) completes at once.

---

## 5. Supplementary runs

<!-- SUPPLEMENT -->
**p = 13, exhaustive phase search over the 2,376 dilation orbits** (`cores --mode exist-orbit`, two workers per Q;
logs `logs/4p3_review/elliptic_exist_p13_*`).
- **No core in the elliptic ansatz at p = 13, on either Q [M].** All four workers finished exhaustively (1,188 orbit
  representatives each; 3,078 s and 3,089 s for ℤ₂², 3,150 s and 3,155 s for ℤ₄). That is every realisation of every
  gate solution and every shift.
- The estimate (§3.4) predicted about 2^4 cores in the whole ansatz from 2^14.8 realisations at about 2^−10 each.
  The exact zero on ℤ₂² is another case of the phase-gate estimate erring high.
- The answer's own statement (one selected quadruple failed) is therefore strengthened to an exact zero on both Q.
- So the ansatz yields cores at p = 3, 5 and 7 only, among all primes up to 97 [M].

## 6. Reproduce

All runs are deterministic. Seeds are fixed in the commands and logged. Run from the repository root.

```sh
# Task 1
python3 theory/4p3_review_nf.py --q v4 --mmax 199 --naive-mmax 9 --seconds 1200 --out logs/4p3_review/nf_v4.json
python3 theory/4p3_review_nf.py --q z4 --mmax 199 --naive-mmax 9 --seconds 1200 --out logs/4p3_review/nf_z4.json
python3 theory/4p3_review_qxcheck.py --out logs/4p3_review/qxcheck.json --seconds 1200
python3 theory/4p3_review_repair.py --out logs/4p3_review/repair_table.json
# Task 2 (needs Kumar's repository at /home/user/vatskishan03/hadamard, or --kumar DIR)
python3 theory/4p3_review_elliptic.py examples --out logs/4p3_review/elliptic_examples.json
python3 theory/4p3_review_elliptic.py cores --mode random --p 5 --q v4 --seconds 300 --pick-seed 2026 \
    --exclude logs/openai_4p3/core_p5_v4.json --out logs/4p3_review/elliptic_random_p5_v4.json   # same for p5_z4, p7_v4
python3 theory/4p3_review_elliptic.py cores --mode random --p 7 --q z4 --seconds 240 --pick-seed 2026 \
    --out logs/4p3_review/elliptic_random_p7_z4.json
python3 theory/4p3_review_elliptic.py collect --examples logs/4p3_review/elliptic_examples.json \
    --random logs/4p3_review/elliptic_random_p5_v4.json,logs/4p3_review/elliptic_random_p5_z4.json,logs/4p3_review/elliptic_random_p7_v4.json,logs/4p3_review/elliptic_random_p7_z4.json \
    --out logs/4p3_review/elliptic_cores_for_border.json
python3 theory/4p3_review_border.py --out logs/4p3_review/border.json --scramble 3 \
    --elliptic logs/4p3_review/elliptic_cores_for_border.json --scratch /tmp
python3 theory/4p3_review_smallcores.py --mods 2,2,3,3 --q v4 --seconds 600 --tmin 6 --seed 11 --out logs/4p3_review/smallcores_v4_3x3.json
python3 theory/4p3_review_smallcores.py --mods 4,3,3 --q z4 --seconds 600 --tmin 6 --seed 12 --out logs/4p3_review/smallcores_z4_3x3.json
python3 theory/4p3_review_smallcores.py --mods 2,2,3,3 --q v4 --slice "1,0,0,0;0,1,0,0;0,0,0,2;0,0,1,0" --out /tmp/slice.json   # 7 slices: logs/4p3_review/smallcores_slices.log
# Task 3
python3 theory/openai_4p3_search.py --primes 23,29,31 --phase-max-p 0 --seconds 540 --out logs/openai_4p3/next_seeds.json
python3 theory/4p3_review_elliptic.py coverage --out logs/4p3_review/elliptic_coverage_check.json
python3 theory/4p3_review_elliptic.py gate --primes 3,5,7 --out logs/4p3_review/elliptic_gate_small.json
python3 theory/4p3_review_elliptic.py gate --primes 11,13,17,19,23,29,31 --out logs/4p3_review/elliptic_gate.json
python3 theory/4p3_review_elliptic.py gate --primes 37,41,43,47,53,59,61,67,71,73,79,83,89,97 --out logs/4p3_review/elliptic_gate_ext.json
python3 theory/4p3_review_elliptic.py estimate --primes 5,7,11,13,17,19,23,29,31 --N 4000 --out logs/4p3_review/elliptic_estimate_small.json
python3 theory/4p3_review_elliptic.py estimate --primes 37,41,43,47,53,59,61,67,71,73,79,83,89,97 --N 4000 --out logs/4p3_review/elliptic_estimate_mid.json
python3 theory/4p3_review_elliptic.py estimate --primes 149,257,331,367,457,467,521,571 --N 4000 --out logs/4p3_review/elliptic_estimate_targets_a.json
python3 theory/4p3_review_elliptic.py estimate --primes 179,271,359,389,461,509,523,599 --N 4000 --out logs/4p3_review/elliptic_estimate_targets_b.json
python3 theory/4p3_review_elliptic.py phase-estimate --primes 5,7,13,17,19,23,29,31,37,43,53,61,149,179 --q v4 \
    --reps 2 --N 3000 --exact-max-p 13 --out logs/4p3_review/elliptic_phase_estimate.json
for q in v4 z4; do for w in 0 1; do
  python3 theory/4p3_review_elliptic.py cores --mode exist-orbit --p 13 --q $q --part $w/2 --seconds 3300 \
      --out logs/4p3_review/elliptic_exist_p13_${q}_part$w.json
done; done
```
