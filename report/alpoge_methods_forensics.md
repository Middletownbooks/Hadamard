# What the original matrices say about how Alpöge et al. found them

Date: 2026-10-09.

**Data.** The analysis covers all 12 orders:
- the authors' own (unnormalised) matrices, supplied by the user, for 8 orders: v = 167, 179, 311, 347, 359, 443, 479, 491;
- their record sequences for v = 223, 283 and 487;
- for v = 419, records rederived from the normalised hoa64 copy (§1.1). These are exact up to, at most, a common translation of all four sequences, which none of the statistics below depend on.

**Files.** Scripts and data are in `forensics/` (see `forensics/README.md`).

**Tags.**
- **[M]** measured on the files.
- **[P]** proved here.
- **[E]** Gaussian constraint-bit estimate.
- **[S]** secondary public source, not verified here.
- **[I]** interpretation.

## Bottom line

1. **The search method has not been published [S].** The release is the matrices plus a decoder.
   - Secondary accounts describe it as a sed-obfuscated shell script with twelve records and five builder routines ([VibeMathed](https://vibemathed.com/problem/hadamard-matrix-of-order-668)).
   - The five `func` letters in the files (g, h, u, v, w) are presumably those five builders.
   - Credit is given to Levent Alpöge, Philippe Voinov, Saul Reynolds-Haertle and Claude ([Epoch AI](https://epoch.ai/frontiermath/open-problems/hadamard), [arXiv 2609.18380](https://arxiv.org/abs/2609.18380)).
   - A news report says the specifics await a technical report ([pedaily](https://news.pedaily.cn/202608/567669.shtml)).
   - I could not open any primary page from this container; these are search-result summaries.
2. **The builders are textbook arrays [M].**
   - `g` is the Goethals–Seidel array, with R the back-identity.
   - `h` is the bordered Goethals–Seidel array (b = 1).
   - `u`, `v`, `w` are group-developed Goethals–Seidel arrays with an inflated border of 4b rows and columns, over:
     - Z₂×Z₂×Z₈×Z₁₁ (u, b = 7);
     - Z₂×Z₃×Z₅₇ (v, b = 5);
     - Z₂×Z₂×Z₁₁₉ (w, b = 3).

   The novelty is in finding the four sequences, not in the arrays.
3. **Every record sits in a sharply restricted search space, and the restriction is chosen order by order [M].**
   - **Multiplier orbits, for 10 of 12 orders.** The groups have d = 3–15 and give 42–128 orbits per block. In 6 of the 9 orbit-type records in the authors' own coordinates, the orbits are centred exactly at the origin; 443, 487 and 491 are shifted by 2.
   - **167 and 179 use no multiplier group. They use a half-shift template.**
     - In natural Goethals–Seidel coordinates every block is x = [p | q | −p | q] with |p| ≈ h/2 (h = (v−1)/2), apart from 1–2 entries per block (5 per core).
     - The chance probability is 10^−71.4 for 167 and 10^−78.0 for 179.
     - This is the "half-shift structure" that `REPORT.md` could not reproduce. Its tests summed x(i)·x(i+h), which a two-step sign pattern cancels to about 0.
   - **The two cores whose group contains Z₁₁ (359 and 419) are both nearly Legendre-twisted on that factor.** They are almost anti-invariant under the non-squares of Z₁₁:
     - 359: 116 of 128 orbit pairs, p ≈ 10^−21. Equivalently, x(g, z) = ε(g)·(z/11) for z ≠ 0 in 116 of 128 Z₁₁-cosets.
     - 419: 48 of 56 orbit pairs, p ≈ 10^−7.

     No such twist appears for the other primes tested (7, 13, 17, 19) or in any other core.
   - **The b = 0 cores have individually flat blocks** (`REPORT.md` §0; not re-measured here).
4. **What the template buys [P/E].**
   - Fixing σ(y) = x(y)·x(y+h) splits the b = 1 equations exactly into two independent problems on Z_h: one periodic, one negacyclic.
   - Their supports are complementary arcs, so both problems are nearly aperiodic.
   - At v = 167 each half has about 165 sign bits and 136–138 constraint bits; at 179, 177–179 sign bits and 148–150 constraint bits. The unrestricted core at 167 has 664 sign bits and 333 constraint bits. That puts each half in the same size class as the orbit slices used at 419 and 443.
   - Parity rules out four exact arcs [P]. The published exceptions sit where parity allows them.
   - At 167 and 179 every multiplier slice of useful size is predicted empty (log₂E between −56 and −17, from `logs/b1_mult_validation.log`). That is presumably why a non-orbit reduction was used there.
   - The template's expected solutions per half fall:

     | v | expected solutions per half |
     |---|---|
     | 167, 179 | 2^31 |
     | 419 | 2^6 |
     | 443 | 2^2 |
     | 523 | 2^−13 |

     This fits the switch to multiplier slices at 419 and 443. It also means **the template does not reach the open orders** [E].
5. **Not inferable from the outputs [I]:**
   - which solver was used (local search, k-tree/MITM, SAT or other);
   - how much compute;
   - who did what.

## 1. Provenance and verification [M]

| file (user-supplied) | header | HHᵀ = NI | normalising gives the hoa64 copy | cross-checks |
|---|---|---|---|---|
| `matrix_H668_v167_h_b1` | v=167 N=668 func=h b=1 m=166 | yes | yes | identical to the matrix in google-deepmind/formal-conjectures commit 91ebb29 (`Hadamard668.lean`) |
| `matrix_H716_v179_h_b1` | v=179 … func=h b=1 | yes | yes | |
| `matrix_H1244_v311_g_b0` | v=311 … func=g b=0 | yes | yes | |
| `matrix_H1388_v347_v_b5` | v=347 … func=v b=5 m=342 | yes | yes | core = `data/orig_bge3` core, up to the inversion convention |
| `matrix_H1436_v359_u_b7` | v=359 … func=u b=7 m=352 | yes | no | core = `data/orig_bge3` core, up to the inversion convention |
| `matrix_H1772_v443_h_b1` | v=443 … func=h b=1 | yes | no | core = `data/core_v443_b1.txt` up to translation and the R convention |
| `matrix_H1916_v479_w_b3` | v=479 … func=w b=3 m=476 | yes | no | core = `data/orig_bge3` core, up to the inversion convention |
| `matrix_H1964_v491_g_b0` | v=491 … func=g b=0 | yes | no | core = `data/core_v491_b0.txt` up to translation and the R convention |
| `seq_v223_g.txt` | v=223 N=892 func=g b=0 | rebuilt with the `g` layout: yes | – | records, not a matrix |
| `seq_v283_g.txt` | v=283 N=1132 func=g b=0 | rebuilt with the `g` layout: yes | – | byte-identical to `data/seq_v283_g.txt` |
| `seq_v487_g.txt` | v=487 N=1948 func=g b=0 | rebuilt with the `g` layout: yes | no | core = `data/core_v487_b0.txt` after a common translation by 2 |
| `seq_v419_h_rederived.txt` | v=419 N=1676 func=h b=1 | rebuilt with the `h` layout: yes | by construction | rederived from `data/hadamard_1676.csv.gz` (§1.1) |

So the repo's existing cores (hoa64-derived, plus `data/orig_bge3`) describe the same matrices, re-labelled. For 668, 716, 1244 and 1388 this was checked through the normalised hoa64 copies; for the others, through the cores.

What hoa64 changed:
- For 443, 487 and 491, the copies are the authors' matrices with all four sequences translated by 2. This moves the multiplier fixed point from 2 to 0.
- For 359 and 479, the cores are identical to the authors'; the difference lies outside the core.
- 347, 668, 716 and 1244 are pure normalisations.

### 1.1 Rederiving the 419 records from the normalised copy [M]

The `h` builder's border is the same in all three b = 1 originals, so it does not depend on the records. Normalisation is hoa = diag(H[:,0])·H·diag(H[0,0]·H[0,:]), and H's first row and first column belong to that fixed border. So the original is recovered exactly as

  H = diag(col₀)·hoa·diag(H₀₀·row₀),

using the builder's own row₀ and col₀ (`forensics/hbuilder.py`).

Validation:
- On the hoa64 copies of 668 and 716 this returns the authors' originals exactly.
- On 1772 it returns an exact `h`-builder matrix, but with the sequences translated by 2, which is hoa64's relabelling.
- On 1676 the result is an exact `h`-builder matrix and Hadamard. Its records are in `forensics/originals/seq_v419_h_rederived.txt`.

The remaining uncertainty is a possible common translation. hoa64 re-centred exactly the records whose fixed point is 2, and the rederived 419 has its fixed point at 0. So the authors' 419 record is centred at either 0 or 2. Every statistic in this report is translation-invariant.

## 2. The builders [M]

**`g` (b = 0).** The layout was inferred block by block from H(1244). It is reproduced exactly by H(1964), and it rebuilds Hadamard matrices from the 223 and 283 records:

```
[  A     BR     CR     DR  ]
[ -BR    A     DᵀR   -CᵀR  ]
[ -CR  -DᵀR    A      BᵀR  ]
[ -DR   CᵀR  -BᵀR     A    ]
```

- A, B, C, D are circulant.
- R is the back-identity (R_{ij} = 1 iff i + j = m − 1).
- The record stores A's first row and, for B, C and D, the first row of BR etc. That row is g(s) = f(−1−s), where f is the circulant's own first row.

"Natural coordinates" below always undo this reversal.

**`h` (b = 1).** Inferred from H(668); 716 and 1772 are rebuilt exactly from their records.
- Four border rows and columns come first. They are record-independent.
- The 4×4 corner has rows (− + + −), (+ − + −), (+ + − −), (− − − −).
- Border rows 0–3 are constant on the core column-blocks, with signs (−,−,−,+), (−,−,+,−), (−,+,−,−), (+,−,−,−).
- Border columns 0–3 are constant on the core row-blocks, with signs (+,−,−,+), (+,−,+,−), (+,+,−,−), (−,−,−,−).
- The core is the same Goethals–Seidel array as `g`.
- Row sums of the records are (±2, 0, 0, 0).

**`u`, `v`, `w` (b = 7, 5, 3).**
- The core starts at offset 4b.
- Blocks are group-developed over a mixed-radix group, last index fastest.
- Block (0,0) is f(j−i); blocks (0,1..3) are f(i+j).
- The groups are exactly those `scripts/extract_group.py` found: (2,2,8,11), (2,3,57) and (2,2,119).
- The authors' records agree with `data/orig_bge3/cores_original.txt`; B, C and D there are stored as f(−g − (1,…,1)).

## 3. Symmetry census of the records [M]

`forensics/census_originals.py` searches for the largest multiplier group U fixing all four blocks about one common point. It works in natural coordinates and searches every candidate fixed point.

| v | N | b | builder | group | d = \|U\| | orbits/block | U's fixed point | other structure |
|---|---|---|---|---|---|---|---|---|
| 167 | 668 | 1 | h | Z₁₆₆ | 1 | 166 | – | **half-shift template** (§4) |
| 179 | 716 | 1 | h | Z₁₇₈ | 1 | 178 | – | **half-shift template** (§4) |
| 223 | 892 | 0 | g | Z₂₂₃ | 3 | 75 | 0 | flat blocks |
| 283 | 1132 | 0 | g | Z₂₈₃ | 3 | 95 | 0 | flat blocks |
| 311 | 1244 | 0 | g | Z₃₁₁ | 5 | 63 | 0 | flat blocks |
| 347 | 1388 | 5 | v | Z₂×Z₃×Z₅₇ | 6 | 70 | 0 | – |
| 359 | 1436 | 7 | u | Z₂²×Z₈×Z₁₁ | 5 | 96 | 0 | **Legendre(11) twist** (§5) |
| 419 | 1676 | 1 | h | Z₄₁₈ | 15 | 42 | 0 (rederived; 0 or 2) | **Legendre(11) twist** (§5) |
| 443 | 1772 | 1 | h | Z₄₄₂ | 12 | 50 | 2 | – |
| 479 | 1916 | 3 | w | Z₂²×Z₁₁₉ | 4 scalar (12 with the GL(2,2) twist, `bge3_construction.md`) | 128 (46) | 0 | – |
| 487 | 1948 | 0 | g | Z₄₈₇ | 9 | 55 | 2 | flat blocks |
| 491 | 1964 | 0 | g | Z₄₉₁ | 7 | 71 | 2 | flat blocks |

Notes:
- The d and r values agree with the census in `REPORT.md`.
- The new information is the anchoring: the authors' records put U's fixed point at the origin in every case except 443, 487 and 491, which are shifted by 2.
- A common translation of all four blocks is a symmetry of the problem, so the shift by 2 is a storage detail, not a method signal.

## 4. The half-shift template (167, 179)

### 4.1 Statement [M]

Write m = 2h (h = 83 or 89, odd) and work in natural coordinates. Define σ_k(y) = x_k(y)·x_k(y+h), a ±1 function on Z_h. In every block σ_k is a two-step function: −1 on [0, M_k) and +1 on [M_k, h−1]. Equivalently,

  x_k = [p_k | q_k | −p_k | q_k],  |p_k| = M_k,

apart from the exceptions listed below. Here "sign changes" counts sign changes of σ around Z_h, and "log₁₀ P" is the probability of that few changes for a uniform random σ: P(2j changes) = 2·C(h, 2j)/2^h.

| v | block | sign changes (null ≈ h/2) | log₁₀ P | σ = +1 arc | exceptions (position, kind) |
|---|---|---|---|---|---|
| 167 | A | 4 | −18.4 | [43, 82] | 61 hole |
| 167 | B | 6 | −16.1 | [39, 82] | 19 isolated +1; 41 hole |
| 167 | C | 4 | −18.4 | [42, 82] | 39 isolated +1 |
| 167 | D | 4 | −18.4 | [42, 82] | 60 hole |
| 179 | A | 6 | −17.7 | [46, 88] | 23 isolated +1; 66 hole |
| 179 | B | 4 | −20.1 | [44, 88] | 23 isolated +1 |
| 179 | C | 4 | −20.1 | [46, 88] | 43 isolated +1 |
| 179 | D | 4 | −20.1 | [44, 88] | 46 hole |

- **Combined.** log₁₀ P = −71.4 (167) and −78.0 (179).
- **The b = 1 cores from orbit slices show nothing.** 443 (original) gives −1.4 and 419 (rederived) gives −3.9.
- **Anchoring.** All arcs end at h−1 in natural coordinates. In the stored records, B, C and D are reversed (§2), so there the arcs start at 0. The template is anchored to the authors' origin.

**Why `REPORT.md` missed it.** Its statistics were sums, such as agreement of x(i+h) with x(i) and PAF(h). For σ ≈ ±(two equal halves), those sums are about 0. The structure is in the run pattern, not the mean. `scripts/b1_cyclotomy.py` did see "±1 purity" ≈ 0.96–0.98, because an arc is reflection-symmetric about its own centre, but it did not identify the arcs.

### 4.2 Exact decoupling [P]

Let u = x[0, h) and w = x[h, 2h). In Z[z]/(z^{2h} − 1), reduction mod z^h − 1 gives u + w = 2c, and reduction mod z^h + 1 gives u − w = 2d. Here c = u·1_S and d = u·1_{S^c}, with S = {σ = +1}.

The b = 1 condition is

  Σ_k x_k(z)x_k(z⁻¹) = (4m + 4) − 4·Σ_{g∈Z_m} z^g.

Reducing it both ways gives:

- Σ_k PAF_{c_k}(t) = −2 for t ≠ 0 (periodic on Z_h), with Σ_k |S_k| = 2h − 1;
- Σ_k NPAF_{d_k}(t) = 0 for t ≠ 0 (negacyclic on Z_h), with Σ_k |S_k^c| = 2h + 1.

For a fixed σ the two systems share no variables. When the supports are arcs of length ≤ (h+1)/2, the periodic and negacyclic correlations equal the aperiodic ones. In the records, both systems hold term by term in aperiodic form for lags t ≤ 18 (167) and t ≤ 20 (179). Beyond that, the exceptions bring in wrap-around terms.

### 4.3 Parity: exceptions are forced, and where they can go [P, M]

**Lemma.** For a fixed support, PAF_{c}(t) ≡ |S ∩ (S + t)| (mod 2).

If all four S_k are proper non-empty arcs, then

  Σ_k PAF_{c_k}(1) ≡ Σ_k (|S_k| − 1) = 2h − 5,

which is odd. But the sum must be −2. So no b = 1 core is an exact four-arc template.

**Exhaustive check** (`forensics/template_parity.py`). Arc lengths range over h/2 ± 7, and every lag t = 1..n is checked, along with the row-sum parities (|S_A| odd, the others even). Results:
- 0 flips: infeasible.
- 1 flip: feasible, but **only** at the exact centre of an arc or of its complement.

The reason: a centre hole removes 0 or 2 pairs at every lag, so it changes |S| without changing any lag parity. Any other flip changes the parity on a window of lags. That must be cancelled elsewhere, by other flips or by the arcs' own parity pattern.

The records have 5 exceptions per core, not the minimum 1, so the template was not minimal:
- 6 of the 10 exceptions lie within 2 positions of a region centre (offsets −1.5, 0, −2, +0.5, −1, +1.5);
- the other 4 lie exactly 2 positions from an arc end.

Each exception flips the parity of a window of lags:
- a near-centre exception flips only 0–4 lags around h/4;
- an exception 2 positions from an arc end flips almost every lag, from 3 to about n.

Each core has exactly one near-end pair (167: B's hole at 41 and C's extra point at 39; 179: C's 43 and D's 46), whose windows nearly coincide. The near-centre flips and the arcs' own parity pattern take care of the rest, and together every lag ends up even, as it must in a solution.

### 4.4 Cost [E]

`forensics/template_cost.py` uses the published supports, random signs and lattice spacing 2 (the parities are fixed by the support).

| v | periodic half: sign bits − constraint bits | negacyclic half | unrestricted core |
|---|---|---|---|
| 167 | 165 − (136 + ~11 row-sum bits) ⇒ log₂E ≈ +17 | 167 − 138 ⇒ ≈ +29 (its one row-sum condition not subtracted) | 664 − 333 ⇒ +331 |
| 179 | 177 − (148 + ~11) ⇒ ≈ +17 | 179 − 150 ⇒ ≈ +29 (same) | 712 − 362 ⇒ +350 |

**Scaling.** Using plain arcs of lengths (n, n+1, n+1, n), log₂E per half is:

| v | log₂E per half |
|---|---|
| 167 | +31 |
| 179 | +31 |
| 251 | +28 |
| 307 | +22 |
| 419 | +6 |
| 443 | +2 |
| 523 | −13 |
| 739 | −66 |

**Reading [I].**
- The template turns one 333-bit problem into two independent problems of about 140 bits each.
- That is the size class of the orbit slices the same team solved at 419 and 443 (`logs/b1_mult_validation.log`):
  - 419: d = 15, 124 constraint bits, log₂E = +18.3;
  - 443: d = 12, 159 constraint bits, log₂E = +15.4.
- **Why not orbits at 167 and 179: every multiplier slice of useful size is predicted empty there** (same log).
  - 167: {±1}, 84 orbits per block, log₂E = −48.5. The only other group, of order 41, leaves 6 orbits.
  - 179: d = 2, 4 and 8 give −55.8, −30.7 and −17.4.
  - 179's coarse d = 11 slice (18 orbits) is predicted at +24. But a partial join over 3000-element pools found no matches at depth 8, where 2^11 were predicted (`logs/b1_partialjoin_179_d11.txt`).

  So at the two smallest b = 1 orders the usual orbit reduction is predicted empty, while the half-shift template has room to spare. At 419 and 443 the reverse holds: orbit slices work, and the template is down to 2^6 and 2^2 solutions per half.
- **For the open orders this is a dead end.** At v = 523 each half already has log₂E ≈ −13 and about 535 constraint bits. That agrees with the earlier, cruder "hopeless" verdict in `unattempted_2000_3000.md`.

## 5. The Legendre(11) twist in 359 and 419 [M]

`forensics/signed_orbits.py` tests every core for signed relations beyond its stabiliser U. For each multiplier class a with a ∉ U but a² ∈ U, the U-orbits pair up as (O, aO). The test counts pairs with x(aO) = −x(O) ("anti") against a fair-coin null.

**Result.** Exactly one class per core stands out, and only in the two cores whose group has a Z₁₁ factor. In both, it is the class that acts as a non-square on Z₁₁ and trivially modulo U elsewhere.

| v | group | U | class a (action) | anti pairs, per block | total | log₁₀ p (corrected) |
|---|---|---|---|---|---|---|
| 359 | Z₂²×Z₈×Z₁₁ | order 5 (squares on Z₁₁) | 17 (identity on Z₂²×Z₈, non-square on Z₁₁) | A 24/32, B 30/32, C 32/32, D 30/32 | 116/128 | −21.8 (−21.0 over 7 classes) |
| 419 | Z₂×Z₁₁×Z₁₉ | order 15 (squares on Z₁₁, order 3 on Z₁₉) | 7 (non-square on Z₁₁, in U's image on Z₁₉) | A 13/14, B 9/14, C 12/14, D 14/14 | 48/56 | −7.3 (−6.9 over 3 classes) |

The other classes are null, including the non-square classes of the other primes. Examples:
- 419, a = 27 (non-square on Z₁₉): log₁₀ p = −0.4;
- 347 (Z₁₉ non-square): −0.5;
- 479 (Z₇ non-square): 0.0;
- 443 (Z₁₃/Z₁₇ classes): ≥ −1.0;
- every b = 0 core: ≥ −0.3.

**What it means for 359.** U acts on the Z₁₁ coordinate alone, so the relation reads coset by coset: in 116 of 128 Z₁₁-cosets, x(g, z) = ε(g)·(z/11) for z ≠ 0. The other 12 cosets are constant on Z₁₁*.

**Reading [I].** Both searches seem to have been biased toward a Kronecker-like "sequence on G/Z₁₁ ⊗ Paley(11)" form. Under the exact form, x(g, z) = δ(g) at z = 0 and ε(g)·(z/11) otherwise. That reduces each core to sequences on G/Z₁₁, of order 32 for 359 and 38 for 419: for 359, 4·64 = 256 free bits instead of 4·96 = 384. Neither solution is exact: 12 of 128 cosets (359) and 8 of 56 orbit pairs (419) deviate. That is the same "template plus a few deviations" picture as the half-shift cores.

**Relevance to the open orders [I].** Z₁₁ appears in the theory-arm's one surviving b = 1 slice, 573 over Z₂²×Z₁₁×Z₁₃ (n ≈ 167; see `report/theory_arm_2000_3000.md`). Z₁₁ is also available for 599 (b = 5, m = 594 = 2·3³·11) and 669 (b = 9, m = 660 = 2²·3·5·11), and with larger borders for 571, 643, 659 and 719. Whether an exact or near Legendre(11) restriction keeps solutions there is untested.

## 6. What the evidence supports [I]

The records look like one strong generic solver applied to a menu of classical reduced search spaces, chosen order by order:

- orbit slices (multiplier groups, in 479 with a GL(2,2) twist);
- the half-shift template at the two smallest b = 1 orders;
- a Legendre(11)-biased slice at 359 and 419, the two orders with a Z₁₁ factor;
- for b = 0, an objective or pool that yields individually flat blocks.

`REPORT.md` §0.5 found "flat pools + k-tree join" consistent at v = 487, but not proven.

In the Gaussian model, the reduced problems land at roughly 120–210 constraint bits:
- template halves: 136–150;
- 419 and 443 slices: 124 and 159;
- b = 0 slices: 130–207, per `REPORT.md`.

Nothing in the outputs identifies the solver itself. The `func` letters name builders (array types), not search methods.

## Reproduce

```
cd forensics
python3 load_originals.py      # verify every original; compare with hoa64 and the formal-conjectures 668
python3 gbuilder.py            # g layout from H(1244); same in H(1964); rebuild 892/1132/1948 from records
python3 hbuilder.py            # h layout and fixed border; rederive the 419 records from the hoa64 copy (validated on 668, 716, 1772)
python3 census_originals.py    # §3 table
python3 halfshift_b1.py        # §4.1 table (167, 179, 443 originals; 419 rederived)
python3 template_parity.py     # §4.3 parity search
python3 template_cost.py       # §4.4 costs and scaling
python3 signed_orbits.py       # §5: the Legendre(11) twist in 359 and 419, null elsewhere
```

## Sources (secondary; primary pages were not reachable from this container)

- VibeMathed, "Hadamard Matrix of Order 668": https://vibemathed.com/problem/hadamard-matrix-of-order-668
- Epoch AI, FrontierMath open problems, "Hadamard": https://epoch.ai/frontiermath/open-problems/hadamard
- arXiv 2609.18380 (cites "L. Alpöge, P. Voinov, S. Reynolds-Haertle, and Claude-Anthropic, Explicit construction of Hadamard matrices for all unknown orders under 2000"): https://arxiv.org/abs/2609.18380
- Wikipedia, "List of mathematical discoveries by artificial intelligence": https://en.wikipedia.org/wiki/List_of_mathematical_discoveries_by_artificial_intelligence
- pedaily news report (Aug 2026): https://news.pedaily.cn/202608/567669.shtml
- Lean formalisation of H(668): google-deepmind/formal-conjectures, commit 91ebb29 (`FormalConjectures/Wikipedia/Hadamard668.lean`)
