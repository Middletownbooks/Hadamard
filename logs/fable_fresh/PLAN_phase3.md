# Plan: what is left to try at the open Hadamard orders (third idea phase, for an implementer)

## Context

Repo `middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`. Goal unchanged: a Hadamard matrix at an open
order, or a proved closure of a natural family. Deliverable so far: `report/fable_fresh_answer.md`; state file
`logs/fable_fresh/NOTES.md`.

Where we stand (2026-11-14, commit fd52dc6b, nothing running):
- Closed with proofs or exhaustive runs: formal products (§3, §4), Theorem 7 with tiny S (§6), every border of a
  Goethals–Seidel core is a character border (§6.4), the 2396 slice has no border with Q₈, ⟨i⟩/⟨j⟩/⟨k⟩ or any
  order-3 symmetry (§7.1, §7.2), tier-0 slices of the old screen (§8a).
- The earlier screens (b0, b1, k4, and Kumar's censuses) are exhausted on this hardware: cheap slices are empty,
  rich ones cost ≥ 2⁵⁴. Three slices with a non-empty stage 1 remain (5028 b0, 6852 b0, 8108 b1), but their G′ is
  ℤ₃ or ℤ₂, so stage 1 is nearly vacuous and they are 2⁷²–2⁹¹ problems. **Not worth more compute.**
- The one family that is new, cheap to decide exactly, and only screened in a corner: **character borders with a
  general character set T over Q = ℤ₂^a** (§7). Screened so far: a = 4 only, cyclic σ, cyclic S, four primes
  (`theory/fable_tcensus2.py`). The corner produced the 2396 slice (78 orbits per block: too big, and no symmetric
  border).

**The new observation that drives this plan.** For a = 5, 6 there are σ-groups on ℤ₂^a with very few orbits and
non-classical invariant T, and S is then small, so whole slices have 12–40 orbits per block (48–160 sign bits) and
can be decided *exactly*, core included. Examples (hand-checked sizes; the census must recompute them):
- 719 = 64·11 + 15: Q = F₈ ⊕ F₈ with F₈\* × F₈\* (orbits 1, 7, 7, 49), T = the two 3-flats through 1 (15
  characters), S = ℤ₁₁ with U of order 5: **12 orbits per block**. Also Q = F₁₆ ⊕ F₄ with σ of order 15 (orbits
  1, 15, 3, 15, 15, 15), T = a mixed 15-orbit: 18 orbits.
- 1487 = 64·23 + 15 (S = ℤ₂₃, U of order 11): the same two shapes, 12 and 18 orbits.
- 863 = 64·13 + 31 (T = 15 + 15 + 1 under σ₁₅; S = ℤ₁₃, U of order 3): 30 orbits.
- 1579 = 32·49 + 11: Q = F₈ ⊕ F₄ with σ of order 21 (orbits 1, 7, 3, 21), T = F₈ ∪ F₄ (11 characters),
  S = ℤ₄₉ with U of order 21: 20 orbits.
Most of these will die at the Q-level by the big-cell Parseval mechanism of §7. That is fine: the point is a
complete, exact map of the only unscreened cheap family, and a real core search for anything that survives.

Three work items, in priority order. Item 1 is the main one.

---

## Item 1. Exact census of small general-T slices over ℤ₂^a, a = 4, 5, 6

New script `theory/fable_tcensus3.py` (generalises `theory/fable_tcensus2.py`; reuse its `qlevel` meet-in-the-middle
and `fable_tclasses.kind` for "classical or not").

**Enumerate**, for every open prime p (the 64 entries 4p of `data/open_orders_lt10000.json`) and a ∈ {4, 5, 6}:
1. *Group on Q.* All groups Σ = ∏ C_i acting on Q = ⊕ F_{2^{d_i}} (Σ d_i = a, summands with d_i = 1 fixed), C_i a
   subgroup of odd order of F_{2^{d_i}}\* acting by multiplication on its summand, plus their cyclic "diagonal"
   subgroups (odd-order elements of GL(a, 2) are semisimple, so this covers every cyclic σ up to conjugacy; build
   them from field multiplications, do not brute-force GL(a, 2)). Keep those with at most 12 orbits on Q.
   Characters carry the transpose-inverse action (as in `gl_odd_classes` of `fable_tcensus2.py`).
2. *T.* Every Σ-invariant set of characters with |T| ≡ p (mod 2^a), |T| ≤ 2^a, s = (p − |T|)/2^a ≥ 2 an integer;
   record its affine kind (classical = point, flat, flat minus point, or a character-twist of one; keep classical
   ones in the output but flag them, they are Kumar's families).
3. *S and U.* Every abelian group S of order s (cyclic and elementary abelian factors), every subgroup U of
   Aut(S) of the form: units of each cyclic factor, or a subgroup of F_{ℓ^m}\* on ℤ_ℓ^m; any order (odd preferred,
   record parity). Orbits per block r = (#Σ-orbits on Q) · (#U-orbits on S). **Keep r ≤ 64.**
   Skip at once (gate of §8a) any U that is the full unit group of a prime factor ℓ ≥ 5.

**Exact gates per slice, cheapest first** (each is necessary; record which one kills):
- G1 *cell energy*: Σ over the 4·2^a cells of y² = 4(2^a − |T|)s + 4|T| must be reachable with the cell alphabet
  {Σ_O ε_O |O|} (O the U-orbits of S).
- G2 *Q-level* (as in `fable_tcensus2.qlevel`): cells constant on Σ-orbits, Σ_k ŷ_k(ρ)² = 4|T| on T and 4p off T.
  Exact count; save the solutions (format of `logs/fable_fresh/q599_solutions.json`).
- G3 *S-level* (new, the dual compression): w_k(s) = Σ_q x_k(q, s), U-invariant, alphabet {Σ over Σ-orbits of Q of
  ε·|orbit|}; Σ_k |ŵ_k(ψ)|² = 4p for every nontrivial ψ of S, and the trivial character as in G2. Exact
  enumeration in the PAF domain on one lag per U-orbit (integers only).
- G4 *consistency*: the row sums of G2 and G3 solutions must agree.
Output: `logs/fable_fresh/tcensus3.jsonl`, one record per slice (p, a, Σ, T, kind, S, U, r, gate reached, counts).
Restart-safe: skip keys already present; commit every hour (the container restarts when the session idles).

**For every slice that survives G1–G4 with r ≤ 40:**
- *Border.* Admissible block functions = ±1 functions on ℤ₂^a with Walsh support in T: enumerate by fixing the
  values on |T| independent points and solving the linear system (2^|T| candidates; DFS with pruning if |T| > 24).
  Then search borders with the existing engines, adapted to take (a, T, M) as input:
  `theory/fable_border14.py` (Q₈ and the three cyclic groups; clique engine `theory/fable_clique.py`), and a
  Σ-invariant search modelled on `theory/fable_border15.py`/`fable_border16.py` (number of fixed rows = trace of
  the element on the border space; fixed rows first, or orbit rows first when they are few).
- *Core, exact.* New `theory/fable_tcore.py`: for each Q-level solution enumerate per block all lifts (U-orbit
  signs in every cell with the given cell sum), compute the PAF vector on one lag per orbit of Σ × U × {±1} on
  G ∖ 0 (exact integers), and do a 2 + 2 join on the full vector with a 64-bit hash and exact recheck
  (pattern of `theory/h4892_stage1_any.py`). Report a cap as UNDECIDED, never as zero.
- A core plus a border is a matrix: assemble with `fable_border.assemble`, verify H Hᵀ = nI with
  `scripts/verify_h2092.py` (independent code), commit at once with a bundle in `pub/` in the format of
  `pub/H2092/`. **Do not post outside the repository.**

**For survivors with 40 < r ≤ 64:** border as above; cost the core with `theory/open_screen_deficit.py` if its
slice object can express the group, otherwise report r, the Q-level count and the fibre size only.

**Also finish the a = 4 loose ends from the first census:** the remaining primes ≡ 7, 11 (mod 16) (only 599, 907,
971, 1019 were run), and the Q-levels that were alive but never followed up: 907 (s = 56, σ of order 5, |T| = 11,
52–120 orbits, 42 Q-level solutions) and 971 (s = 60, 76–96 orbits, 52 solutions). For those two: dump the
Q-level solutions and run the border engines (b = 44).

**Controls (must pass before any open-order verdict is trusted):**
- The pipeline on ℤ₂³ × ℤ₃ with T = {0, 1, 2, 4} and trivial multiplier group must find Q-level solutions, a
  border and a core, and assemble a verified H(112) (the instance of §7).
- `fable_tcensus3.py` restricted to a = 4, cyclic σ, cyclic S must reproduce `logs/fable_fresh/tcensus2.jsonl`
  for 599 (10 T with 240 solutions at r = 78; zero at r = 30, 42, 60).
- G3 must pass on the real cores in `data/` that have a ℤ₂^a × S structure (e.g. Kumar's 3964 over
  ℤ₂² × ℤ₂₄₇ if present; otherwise on the annealed H(112)/H(116) cores).
- `fable_tcore.py` with a planted quadruple (target = its own PAF sum) must recover it.

## Item 2. Decide the 2396 border without any symmetry, for single-character rows; shape atlas for |T| = 7, 11

Cheap, and it removes the last CP-SAT "UNKNOWN" of §7.1.
- Add to `theory/fable_border14.py` an option `sub=id` (no block symmetry: one image per row) and `mode=single`
  (block functions = the 14 signed characters of T only). Rows: 14⁴ = 38,416, i.e. 19,208 up to sign; 28 rows per
  side. The coupled search (`fable_clique.coupled`, needz = 28, neede = 27) should finish in minutes: after two
  Z-rows the common compatible set is already below 27 [E]. Run for solution 0 (constant profile) and solution
  2160 (non-constant) of `logs/fable_fresh/q599_solutions.json`. Keep the first-row symmetry reduction.
- If that is exhaustive-negative, widen stepwise and stop at the first class that does not finish in 6 hours:
  rows with at most one non-character block (about 2.1·10⁵ rows), then at most two, then all 810,000.
- *Atlas.* For each affine class of 7-subsets and 11-subsets of F₂⁴ (3 and 2 classes, `theory/fable_tclasses.py`)
  and the constant profiles GS(4, 2, 2, 2) resp. GS(6, 2, 2, 0) (use the synthetic M = Γ ⊗ I as in
  `theory/fable_border16.py`): number of admissible block functions; does an orthogonal admissible basis exist at
  all (a Γ-free necessary condition); Q₈, cyclic and (where it applies) single-character borders. One table in the
  report. This says which shapes at a = 4 can ever carry a border, before more census time goes into them.
- Control: the H(108) case (Q₈ none, each cyclic group finds one) and H(116) must come out as in §7.2; `sub=id
  mode=single` must find the known single-character border of H(112).

## Item 3 (only if Items 1 and 2 leave time). Involutions that move Q, for the 2396 border

Completes "a border of this shape has no symmetry at all". Group (Q₈ ⋊ ℤ₃) × Aut_Q, order 24 · 1152 (§7.2).
- Enumerate conjugacy classes of elements g with g² = ±1, g ≠ ±1, whose Q-part is not trivial (block part 1 or i;
  Q-part an involution of the 1152 affine maps). Classes via brute-force conjugation on the 27,648 elements.
- For each class generalise `theory/fable_border16.py` to order 2: rows with xg = ±x count as fixed rows,
  orbits have two rows; with d₊, d₋ the dimensions of the ±1 eigenspaces on the border space, the admissible
  splits are n₂ orbits, d₊ − n₂ and d₋ − n₂ fixed rows of each sign, 0 ≤ n₂ ≤ min(d₊, d₋). Fixed-first or
  orbit-first as in that script; node caps reported as UNDECIDED.
- Control: the class of i (Q-part trivial) must reproduce "none" of the `sub=z1` run.

---

## Not to do (decided, with the reason)

- More time on 5028 b0, 6852 b0, 8108 b1 or on the 28 capped tier-1/2 slices: tiny G′, so stage 1 is not a gate;
  12 GB already failed.
- General (no-symmetry, all 810,000 rows) borders for 2396 beyond the stepwise attempt of Item 2: about 10² core-days
  by the tree count [E], and the slice is out of reach even with a border.
- Classical-shape borders at index k ≥ 10 with Singer × Q₈: reduces to a quaternion sequence of length k − 1 with
  ideal autocorrelation; the only algebraic candidates are single characters, which Kumar proved impossible for
  k ≥ 12; by count about 2⁻¹⁶ expected at k = 16 [E].
- Products, templates, schemes, non-abelian class functions, time–frequency symmetries: closed in the report (ideas
  1–3, 16–18).

## Rules that still apply (from `report/prompt_fable_fresh.md` §4)

- Branch `claude/focused-planck-j1645f` only; before each push `git fetch -q origin claude/focused-planck-j1645f &&
  git merge -q --no-edit origin/claude/focused-planck-j1645f`, then `git push -u origin HEAD:claude/focused-planck-j1645f`.
  No pull requests. Commit at least hourly: the container is reclaimed when the session idles.
- Never `pkill -f`; stop only your own jobs by PID. Long runs restart-safe (checkpoint files in `logs/fable_fresh/`).
- Tags [P], [M], [E]; a cap or timeout is not a zero; every new search calibrated on a solved case first.
- Credit: Kumar for the character-border family and `kt4.c`; EKLO; Alpöge et al.; the OpenAI model for Theorem B.
- Report: add a section "8b. Census of small general-T slices" and update §7.1/§7.2, the bottom line, §9 and §10
  of `report/fable_fresh_answer.md`; state file `logs/fable_fresh/NOTES.md`.

## Verification (end to end)

1. `python3 theory/fable_clique.py` (self-test) and the controls listed under each item pass.
2. `logs/fable_fresh/tcensus3.jsonl` has one record per enumerated slice; a summary script prints, per prime, the
   number of slices by gate reached and the survivors with r, |T|, kind.
3. For every survivor with r ≤ 40: a border verdict (found / none in the named classes / undecided) and a core
   verdict (found / exactly none / undecided at a named cap).
4. Any matrix: `scripts/verify_h2092.py`-style independent check of H Hᵀ = nI, novelty against
   `data/open_orders_lt10000.json`, committed with a `pub/` bundle.
5. Report sections updated; everything pushed.
