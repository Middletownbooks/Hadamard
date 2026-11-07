# Prompt (high effort, targeted): decide the 4892 slice exactly (v = 1223, b = 1), then clear the remaining cheap slices

*Self-contained, for a fresh Claude Code session on `middletownbooks/hadamard`, branch `claude/focused-planck-j1645f`.*

**Why this session exists.** It is the first open-order slice in a long while with a **non-empty stage 1** at an
affordable cost. Deciding the whole slice exactly is cheap: estimated at about 2^31–2^32 elementary operations.

- **If a core exists:** assemble H(4892), verify it independently, and publish it.
- **If not:** the slice is closed exactly. Then run the second part (§4), which clears the other cheap slices with the
  same tools.

**Time box:** about 12 hours of wall-clock time. Commit after each step.

**Read first:**
1. this file;
2. `logs/open_screen/followup_20261106/README.md` (what was run on 2026-11-06, and its results);
3. `report/open_screen_answer.md` §0–§4 (the screen, its cost model, its stage-1 tool);
4. `report/paley_523_answer.md` §0 (the b = 1 route that built H(2092): staged slice, border, verification);
5. `pub/H2092/` (bundle format, `verify.py`, `border.json`).

## Rules

- **Pushing:** before every push run
  `git fetch -q origin claude/focused-planck-j1645f && git merge -q --no-edit origin/claude/focused-planck-j1645f`,
  then `git push -u origin HEAD:claude/focused-planck-j1645f`.
- **Branches:** never push to `main` or to `codex/order-3452-compression`. Open no pull requests.
- **Commits:** trailers follow the existing commits. Scripts go in `theory/h4892_*.py`, logs in `logs/h4892/`.
- **Processes:** never kill processes you did not start. Never `pkill -f` with a pattern that matches your own
  shell; kill by PID. Never `rm` a path built from an unchecked shell variable.
- **Restart safety:** arm one `send_later` check-in while workers run. Runners resume from their output files.
- **Tags:** mark claims [P] proved, [M] exact computation (script named), [E] estimate. A cap or a timeout is never a
  zero.
- **Verification:** a core counts only if verified by integer PAF on all lags by a script sharing no code with the
  search. A Hadamard matrix counts only if H Hᵀ = 4892·I is verified exactly by a separate script that shares no code
  with the construction. Use the pattern of `scripts/verify_h2092.py`.
- **Compute:** 4 cores and about 15 GB of RAM. Install `numpy`, `numba` and `ortools` with pip if they are missing.

## 1. The slice (exact)

**The family and the group.**
- **Family b = 1** (Legendre quadruple, bordered Goethals–Seidel, as for H(2092), H(2292) and H(2588)): four ±1
  sequences x₀, …, x₃ on G with
  - Σₖ PAFₖ = 4|G| = 4888 at 0, and −4 at every other lag;
  - |G| = v − 1 = 1222.

  The border then gives H(4v) = **H(4892)**. 4892 is open (`data/open_orders_lt10000.json`).
- **The group:** G = G′ × ℤ₁₃, with G′ = ℤ₄₇ × ℤ₂. G ≅ ℤ₁₂₂₂ by the CRT.

**The slice key** (`theory/open_screen_stage1_exact.py` format): `["b1", 1223, [], [47, 2], 13, 4, [[0, 23]]]`. That
is, q = 13, s = 4, f = (q − 1)/s = 3, and S = U × M:
- **U** = the subgroup of order 23 of ℤ₄₇^*, acting on axis 0 of G′ (the squares mod 47). Its orbits on G′ are
  {0} × {b}, QR × {b} and QNR × {b} for b ∈ ℤ₂: **6 orbits**.
- **M** = the subgroup of order 3 of ℤ₁₃^* = {1, 3, 9}, acting on ℤ₁₃. Its cosets are C₀ = {1, 3, 9},
  C₁ = {2, 5, 6}, C₂ = {4, 10, 12} and C₃ = {7, 8, 11}.
- **Orbits of S** on G: (O, 0) and O × C_j for each U-orbit O. So **30 orbits per block** (6 × (1 + 4)), and 120
  signs in all.
- **Sequences:** each xₖ is constant on these orbits. Write αₖ(O) for its value on (O, 0) and c_{k,j}(O) for its value
  on O × C_j.

**Stage 1 (compression along ℤ₁₃).**
- **Definition:** yₖ(b) = Σ_z xₖ(b, z) = αₖ(O) + 3 σₖ(O) for b ∈ O, with σₖ(O) = Σ_j c_{k,j}(O) ∈ {−4, −2, 0, 2, 4}.
  - The alphabet is ±1 + 3·{−4, …, 4}, which is {±1, ±5, ±7, ±11, ±13}.
  - α is determined by y: α ≡ y (mod 3), with (y − α)/3 even.
- **Targets:** Σₖ PAF_{G′}(yₖ) = 4840 at 0 and −52 at every other lag of G′.
- **Index order:** G′ elements are indexed (a, b) ↦ 2a + b, as `itertools.product(range(47), range(2))`.

**Stage 2 (the lift).**
- Given y, the unknowns are the class signs c_{k,j}(O), with Σ_j c_{k,j}(O) = σₖ(O) fixed. So each (block, orbit) has
  C(4, (σ + 4)/2) choices.
- The remaining condition is the full equation on G: −4 at every nonzero lag. It suffices to check one lag per S-orbit
  class.

**What is known** (`logs/open_screen/followup_20261106/`):
- **[M] One method only.** Method B of `theory/open_screen_stage1_exact.py --strong` finds **384 stage-1 solutions**
  (row pattern (0, 0, 0, 2) up to its normalisation) in 55 s, and stores 20.
  - Method A ran out of memory. The planted controls did not finish within 900 s.
- **[M] Independent checks of the stored solutions:**
  - five of them satisfy all 94 stage-1 lag equations, with row sums (−2, 0, 0, 0) and every entry in the alphabet;
  - all 20 are constant on the U-orbits.
- **[E] Screen numbers** (`logs/open_screen/table.jsonl`):
  - stage-1 log₂E = 14.8; E per y = 2^−12.5; C₂ = 58.3 bits; stage-2 fibre per y ≈ 2^45.9;
  - log₂ E_tot = 2.3; cost 2^27.5 per core.
  - If 384 is the complete count, the slice holds about 2^−3.9 ≈ 0.07 expected cores.
  - The estimates have run optimistic at the open orders. **Decide; don't extrapolate.**

## 2. Steps for 4892

**Step 1. Stage 1, complete and confirmed (≤ 2 h).**
- **Two independent complete enumerations.**
  - (a) Method B, with its planted controls run to completion: `--methods B --strong --planted 4` and a long enough
    limit.
  - (b) Your own enumerator, sharing no code with the stage-1 tool. Per block, list the U-invariant y with the
    block's row sum, under the energy bound. Join them 2 + 2 on the 94-lag PAF vector; the six U-orbit values per
    block make this small.
  - Both must agree on the solution set **as a set of explicit (y₀, y₁, y₂, y₃)**. State the symmetries modded out,
    and expand to every ordered quadruple that is a distinct stage-1 solution of the slice.
- **Symmetries.** Normalisations by block order, block negation or translation on G′ preserve the problem only if
  they preserve the slice. List exactly which are used.
- **Planted controls:** recover planted y by both enumerators.
- **Save** `data/h4892_stage1.json`, the complete list, with an independent check of every entry.

**Step 2. Exact stage 2 for every stage-1 solution (≤ 6 h).**
- **The search.** For each y, enumerate each block's class-sign assignments consistent with σₖ(O). That is at most
  6⁶ ≈ 4.7·10⁴ per block, and typically far fewer. Then join 2 + 2 on the exact full PAF vector over G, or on the
  class-sum vector: one representative lag per S-orbit of G, with one equation per orbit pair under negation.
  - Use a hash on integer vectors, with a full recheck of every hit.
  - Alternatively, use the S/D split of `report/paley_523_answer.md` §3 if it is cheaper. The direct join should
    suffice here.
- **The cost bound, as a check.** The per-y fibre is about 2^46 points over four blocks, i.e. about 2^11.5 per block,
  so the 2 + 2 join is about 2^23 pairs per y. 384 or more y give about 2^31–2^32 in all. Measure the real list sizes,
  and if they are far larger, say so and partition.
- **Planted test first:** take a random y from the list, plant random class signs, set the stage-2 target to the
  planted PAFs, and require recovery. Do this for 4 plants.
- **Run all y.** Log per y the list sizes, the pairs visited and the hits.
- **Result:** either a list of cores, or **exactly empty** for the whole slice [M].

**Step 3. If a core is found.**
- Verify it on all 1221 nonzero lags of ℤ₁₂₂₂ (CRT), with independent code.
- Assemble the bordered GS matrix as for H(2092). Reuse the border construction of `pub/H2092/` or
  `theory/q523_finish.py`, adapted to length 1222.
- Verify H Hᵀ = 4892·I with a separate verifier, as in `scripts/verify_h2092.py`.
- Commit at once. Write `pub/H4892/` (sequences, border, normalised matrix, `verify.py`) in the format of
  `pub/H2092/`, and add a README results row.
- Do not post anywhere outside the repository; the owner decides that.

**Step 4. If the slice is empty:** record it [M]. Then check the same group at nearby parameters:
- the other U of order dividing 46, and the other s ∈ {2, 3, 4, 6, 12} with 24–160 orbits;
- the screen's other 4892 slices (`table.jsonl`, v = 1223).

For each, decide stage 1 exactly, then stage 2 if stage 1 is non-empty and the total cost is ≤ 2^40.

## 3. Calibration warnings

- Every positive estimate decided exactly at the open orders so far was zero, except that this slice's stage 1 is
  not empty. The expected count is small (about 0.07, if 384 is complete), so an empty result is the likely
  outcome. It must still be decided, because the decision is cheap.
- A selected representative is not a decision: Step 2 runs every stage-1 solution.

## 4. Second part: clear the remaining cheap slices (the rest of the time box)

Use `theory/open_screen_stage1_exact.py KEY --methods A,B --strong --planted 2` with longer limits (up to 1 h each,
two at a time). Run your Step 1 enumerator as a second method on anything non-empty.

**Undecided at 900 s on 2026-11-06:**
- `["b1", 1795, [], [13, 3, 2], 23, 2, [[0, 12]]]` (7180)
- `["k4", 751, [4], [17], 11, 2, [[1, 8]]]` (3004)
- `["k4", 1291, [4], [7, 2], 23, 2, [[1, 6]]]` (5164)
- `["k4", 1979, [4], [13, 2], 19, 2, [[1, 12]]]` (7916)
- `["b1", 1198, [], [9, 7], 19, 2, [[1, 6]]]` (4792)

**Tool error:** `["b1", 1473, [], [32, 2], 23, 2, [[0, 8, 3]]]` (5892). Here `b_setup` expects (axis, order) pairs.
Fix it minimally (this U entry appears to carry three numbers), with a regression test that the other keys give
unchanged results.

**The tier-0 slices never run** (from `report/open_screen_answer.md` §6):
- `["b1", 947, [], [11, 2], 43, 14, [[0, 10]]]`
- `["k4", 1039, [4], [7], 37, 4, [["neg", [0, 1]]]]`
- `["k4", 1187, [4], [8], 37, 4, [["neg", [0, 1]]]]`
- `["b1", 1259, [], [17, 2], 37, 12, [[0, 16]]]`
- `["b1", 1423, [], [3, 3, 2], 79, 6, [["neg", [0, 1]]]]`
- `["b1", 1427, [], [23, 2], 31, 10, [[0, 22]]]`
- `["b1", 1893, [], [11, 2, 2], 43, 6, [[0, 10]]]`
- `["b1", 2315, [], [13, 2], 89, 8, [[0, 6]]]`

**Any slice with a non-empty stage 1** gets Steps 1–2 if its estimated total stage-2 cost is ≤ 2^40. Otherwise record
its real stage-1 count and its re-read deficit (`analyse(..., real=True)` in the screen code), and stop on it.

**k4 slices with odd |K|** (v ≡ 7 mod 8) complete with the universal border of Theorem B: certificates in
`logs/openai_4p3/border_certificates.json`, assembly in `theory/4p3_review_border.py` and
`theory/4p3_review_hverify.py`.

## 5. Deliverable

Write `report/h4892_answer.md` with:
- the bottom line (4892: H(4892) built, or the slice exactly empty; the second-part table);
- the methods, scripts, logs and planted controls;
- the next computation.

Commit, push, and stop.
