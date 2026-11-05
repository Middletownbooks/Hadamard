# Prompt (OpenAI model; ideas and proofs, short sessions): an algebraic construction for the k4 bordered family at n = 4p + 3

*Self-contained, for a session with access to the public repository `https://github.com/middletownbooks/hadamard`.*
- **Push target:** you may push only to the branch `codex/order-3452-compression`.
- **The task:** construct, or rule out, an explicit algebraic family of the sequences below over groups of order 4p
  (p prime). One family would give Hadamard matrices at up to 16 open orders below 10000.
- **Testing:** long searches are not your job. Claude Code instances with exact verification tooling will check your
  proposals. Write them so that such an instance can build and test each one without asking you anything.

## 0. Working rules (sessions may time out: checkpoint early and often)

- **Setup:** `git fetch origin codex/order-3452-compression claude/focused-planck-j1645f`. Check out
  `codex/order-3452-compression` and merge `origin/claude/focused-planck-j1645f` into it. If the merge is not clean,
  abort it and read the Claude branch with `git show origin/claude/focused-planck-j1645f:PATH`.
- **Pushing:** push only to `codex/order-3452-compression`. Never push to `main` or to `claude/focused-planck-j1645f`.
  Open no pull requests.
- **Checkpoints:**
  - In your first ten minutes, create `report/openai_4p3_answer.md` with the skeleton of §6, and push it.
  - Push again after every finished proposal, or every 30 minutes, whichever comes first.
  - If you run out of time, the file must say what is finished and what is not.
- **Runs:** keep each computation under about 10 minutes. Put scripts in `theory/openai_4p3_*.py` and logs in
  `logs/openai_4p3/` (command, seed, output).
- **Tags:** mark every claim [P] proved (proof in the text), [M] exact computation (script named) or [E] estimate
  (method named).
- **Matrices:** never claim a Hadamard matrix without an exact integer check of H Hᵀ = N·I by code that does not
  share the construction's code.
- **Credit:** the bordered construction is Theorem 7 of Egan–Kotsireas–London–Ó Catháin [EKLO], in Kumar's form
  (`https://github.com/vatskishan03/hadamard`). The census, the obstructions below and the border tooling are by
  Claude Code in this repository.

## 1. The problem (exact)

**Data.**
- n = 4p + 3 with p an odd prime.
- G is abelian of order 4p: G = ℤ₂² × ℤ_p, or G = ℤ₄ × ℤ_p.
- K = ℤ_p (the subgroup of index 4) and Q = G/K, which is ℤ₂² or ℤ₄.

**Wanted:** four ±1 sequences x₁, …, x₄ on G with

    Σᵢ PAFᵢ(g) = 4|G| = 16p    at g = 0,
    Σᵢ PAFᵢ(g) = −12           for g ∈ K \ 0,
    Σᵢ PAFᵢ(g) = +4            for g ∉ K.

Equivalent form, for the characters χ of G: Σᵢ |x̂ᵢ(χ)|² is
- 4n for every χ nontrivial on K;
- 12 for the three characters nontrivial on Q and trivial on K;
- Σᵢ rᵢ² = 4n for the trivial character, where rᵢ = Σ_g xᵢ(g) are the row sums.

**Completion.** Such a core completes to a Hadamard matrix of order 4n: a Goethals–Seidel array over G plus a border
of width 12.
- The border depends only on the Q-compressions zᵢ(q) = Σ_{k∈K} xᵢ(q + k). These are odd integers with
  Σᵢ PAF_{zᵢ} = 4p + 12 at 0 and 4p off 0 on Q.
- Borders were found for every Q-level class tested (850 classes at other orders; `report/b3_push_answer.md` §7).
  They are built by `theory/bt_border.py`.
- Assembly and verification are in `theory/bt_engine.py` and `theory/bt_control2148.py`, and the bundle format is
  `pub/H3732/`. The order-3452 case is a worked example: `logs/h3452/border.json`.

**Equivalently, by cosets.** Write xᵢ on the four K-cosets as 16 ±1 sequences of length p, a_{i,q} for q ∈ Q. The
conditions become a system of periodic and cross-correlation equations on ℤ_p among these 16 sequences. Writing it out
explicitly is a good first step, and part of the deliverable.

## 2. The 16 open orders of this shape (all n − 3 = 4p, p prime)

| 4n | n | p | p mod 8 | p − 1 |
|---|---|---|---|---|
| 2396 | 599 | 149 | 5 | 2²·37 |
| 2876 | 719 | 179 | 3 | 2·89 |
| 4124 | 1031 | 257 | 1 | 2⁸ |
| 4348 | 1087 | 271 | 7 | 2·3³·5 |
| 5308 | 1327 | 331 | 3 | 2·3·5·11 |
| 5756 | 1439 | 359 | 7 | 2·179 |
| 5884 | 1471 | 367 | 7 | 2·3·61 |
| 6236 | 1559 | 389 | 5 | 2²·97 |
| 7324 | 1831 | 457 | 1 | 2³·3·19 |
| 7388 | 1847 | 461 | 5 | 2²·5·23 |
| 7484 | 1871 | 467 | 3 | 2·233 |
| 8156 | 2039 | 509 | 5 | 2²·127 |
| 8348 | 2087 | 521 | 1 | 2³·5·13 |
| 8380 | 2095 | 523 | 3 | 2·3²·29 |
| 9148 | 2287 | 571 | 3 | 2·3·5·19 |
| 9596 | 2399 | 599 | 7 | 2·13·23 |

The open-order list is `data/open_orders_lt10000.json` (EKLO's list of 86, less Kumar's 5912).

**Small cases for testing:** every n = 4p + 3 with p ≤ 97 is a known order. Use them freely to test a construction:
- p = 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, …;
- n = 23, 31, 47, 55, 71, 79, 95, 119, 127, 151, ….

## 3. What is already known here (do not re-derive)

**Multiplier slices are closed at 599 and 719 for r ≤ 160** [M]. Source: `report/b3_push_answer.md` §0, §4.3 and §5.
- These are sequences constant on the orbits of a multiplier group S acting on ℤ_p and on Q; r is the number of
  orbits per block.
- The census went up to r = 400. Every S that Theorem α leaves has r ≤ 20 or r ≥ 114, and everything with r ≤ 152
  dies at the Q-level, which is the only exact level before G itself.
- What is left: S of order 2 or 4 with r = 225–360, with estimates up to 2^76.7 [E] and no exact level to run.
- So a construction that is constant on cyclotomic classes of small index of F_p^* (Paley-like slices) does not exist
  at 599 or 719. It must be finer, or not multiplier-invariant.

**Obstructions [P]** (proofs in `report/b3_push_answer.md` §3):
- **Theorem α:** no core if S contains an element fixing K pointwise and of order 3 on Q = ℤ₂².
- **Lemma R:** no core if S contains τ that is transitive on N \ 0 for an S-invariant N ≤ K with |N| ≥ 5, and whose
  fixed points form a complement of N. This covers primitive multipliers of ℤ_p.
- **Coset-Parseval** for this family.

**Calibration [M]:** `theory/bt_brute_count.py` counts cores exactly for small structures (2^r orbit signs per block,
r ≤ 18 or so), and `logs/bt/calib/calib_small.jsonl` holds 561 small structures (n ≤ 31) with exact counts. These are
data for pattern-finding: look at which structures have cores.

**Successes of related families** (`report/b3_push_answer.md` §1, §6.2):
- Kumar's matrices of orders 3964, 2636 and 2148 come from this construction family, but at |G| with a richer
  subgroup lattice.
- Kumar's 839 is built over the Galois ring ℤ₄³ × ℤ₁₃ (his repository, `search/NOTES.md`).

**Our b = 1 core of length 522** (`pub/H2092/`, v = 523) gives H(2092). 2095 has p = 523: if you see a way to use
known objects of order p or 4p as ingredients, say how; see §4, direction 3.

## 4. Directions worth considering (you are not limited to these)

1. **Cyclotomic or Gauss-sum ansatz with a twist.**
   - Define the 16 coset sequences a_{i,q} from cyclotomic classes of order e of F_p (e | p − 1), or from characters
     of F_p, combined over Q by a fixed ±1 pattern (for instance Hadamard matrices of order 4 or 16), with a few free
     discrete parameters.
   - Compute the character sums in closed form (Gauss and Jacobi sums, cyclotomic numbers) and find the conditions on
     p and the parameters under which all equations of §1 hold.
   - The finer-than-multiplier requirement of §3 suggests allowing different class unions in different cosets, or
     non-invariant corrections at a few points (as in Paley's construction, where 0 is special).
2. **Q = ℤ₄ with quartic characters.** For p ≡ 1 (mod 4), let the ℤ₄ structure of Q interact with the quartic
   residue character of F_p. For p ≡ 3 (mod 4), with the quadratic character and i ∉ F_p.
3. **Recursive or product links.** Known objects of length p or 4p as ingredients:
   - Legendre pairs, b = 1 or b = 0 Goethals–Seidel cores, Williamson or T-sequences, conference matrices;
   - our b = 1 core at 523 for the order 2095.
4. **Relative difference sets, or Galois rings** (Kumar's 839 route) adapted to ℤ₂² × F_p or ℤ₄ × F_p.
5. **An obstruction.** If a natural class (for example "constant on the classes of a subgroup of F_p^* on each
   coset, with at most c exceptional points") cannot contain cores, a proof [P] is a useful answer: it tells the
   searches what to skip.

**A good working method:**
- solve small p exactly within a candidate ansatz (p = 5, 7, 11, 13, 17, 19; `theory/bt_brute_count.py`, or a small
  exhaustive or CP-SAT run under 10 minutes);
- look at the solutions' structure;
- generalise, and prove the generalisation with character sums;
- test it at the next few p.

## 5. What to deliver: up to three proposals, ranked

Each proposal has:
1. **Definition:** the construction as an explicit rule: for each i and each g ∈ G, the value xᵢ(g) as a function of
   p and named parameters. Give a reference implementation `theory/openai_4p3_build.py` that outputs the four
   sequences for given p and parameters.
2. **Why it works:** a proof [P], or the reduced conditions on p and the parameters with the step that remains open.
3. **Exact checks [M]:**
   - the four sequences verified at every small p where the conditions hold, by direct integer PAF;
   - where you can, the bordered Hadamard matrix verified at one small order.
4. **Coverage:** for each of the 16 orders of §2, whether the conditions hold (holds / fails / unknown), with the
   computation.
5. **A test for Claude:** what a Claude max-effort instance should run to confirm the proposal or kill it, with pass
   and kill criteria fixed in advance.

## 6. The answer file: `report/openai_4p3_answer.md`

It has these sections:
- **Bottom line:** the best proposal, its coverage of the 16 orders, and confidence.
- **The explicit coset system of §1,** written out.
- **One section per proposal,** with items 1–5 of §5.
- **Ideas considered and rejected,** one line each with the reason.
- **What a Claude instance should run first:** exact commands, the expected outcome, and the kill criterion.
- **Unfinished items,** if the session ends early.

Commit with a message that names the proposals, and push to `codex/order-3452-compression`.
