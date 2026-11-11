# Working notes, Fable fresh-idea session (state file; survives container restarts)

Started 2026-11-09. Deliverable: `report/fable_fresh_answer.md`. Scripts: `theory/fable_*.py`.

## Orientation done (2026-11-09/10)

- Kumar's repository is still at d885022 (2026-11-06); nothing newer. Clone: /home/user/vatskishan03/hadamard
  (re-clone with GIT_LFS_SKIP_SMUDGE=1 git clone --depth 50 if missing).
- Read: open_screen_answer §0-2, b599_answer §0, alpoge_gap_analysis §0, openai_identity_answer (other branch),
  Kumar FRONTIER.md, skew/README.md, skew/PALEY_SKEW.md, NOTES §34, search/NOTES §102, §107, §108.

## Ideas considered and dropped (reason)

1. Ito / dicyclic (negaperiodic Golay pairs of length 2p): 4p signs, p-1 equations; Gaussian count negative
   for p > ~163 [E]. Dropped.
2. Legendre pairs over a group of order 2p-1 (2 blocks): 4 signs per equation, negative count. Kumar screened
   the Spence (4-block) form over the same groups. Dropped.
3. Norm-equation view: GS over Z_p with multiplier group U <=> beta in O_F^4 (F = fixed field, degree e),
   sum beta_f conj(beta_f) = p, beta_f = s_f/2 mod the prime over p; +-1-ness is automatic [P, short proof:
   Parseval + oddness]. Smith type of the 4x4 similitude at p is (0,0,e,e). Solutions are the minimum
   vectors of p-neighbours of I_4. No algorithmic gain found: bucketing by moments mod p costs as much as it
   saves. Dropped as a search method; keep as language.
4. Torus plug-in: G = F_p^* (order p-1, Theorem 5 group), base block X(x) = chi(x+c) (Sidelnikov). Exact:
   PAF = -1 - kappa(y) off 1; four of them = H(4) x Paley(p+1) seen with the split torus (border 8). Target
   needs Delta PAF = 4 kappa(y): energy 4p must move into the single character kappa; D-part row sums must
   reach sum of squares p, so >= ~25-50 sign changes, and then 299 equations must be restored: no solutions
   near the base point by count [E]. Dropped.
5. Uniform cyclotomy (punctured lines of F_q^2) in Theorem 7 fibres: energy identity forces q <= 5 [P, worked
   in session; same flavour as Theorem T]. Dropped.
6. 10012 = 4*2503 over Z_2^2 x F_625 (EKLO family q = 25): Kumar already has it (search/NOTES line 6466, 7061:
   F_25 x F_25 exactly empty; q = 25 "2% per run, not run"). Dropped.
7. 4792 = 8*599 via Kumar's Spence plug-in Z_7 x Z_19 x F_9, |H| = 36: his calibrated cost 1.5e3 core-hours;
   we have 288 core-hours. Dropped unless a 5x gain appears.
8. Skew T-structures (I, X, -eX^T, Y skew) for Paley blocks, t = 2a^2+1: Kumar exhausted t = 9, 19; count
   heuristic ~1e-3 at t = 33. Dropped.
9. Small skew v (89, 101, 107, 149): no multiplier structure beyond the trivial one (89: only |U| = 11, 31
   bits); unstructured needs 2^311 space. Dropped.

## Idea under test: multiplication of GS quadruples by T-matrices ("octonion" / complex-OD form)

Goal: a theorem "T-matrices of order t  x  GS quadruple of a restricted TYPE over Z_v  =>  GS quadruple over
Z_t x Z_v", for a type that searches can reach at v = 223, 283, 491, 571 (GS cores exist there: AVR, Kumar).
Would settle 2676 = 4*3*223, 4460 = 4*5*223, 5892 = 4*3*491, 6852 = 4*3*571, 7924 = 4*7*283.

Facts so far [P unless marked]:
- Octonion product in complex coordinates is monomial (each z_j = sum of 4 terms +-x_i^(*) y_k^(*)), so
  Z = T * W is a GS quadruple whenever the conjugation pattern on the T side respects supports.
- General T needs supports symmetric => impossible for odd t > 1 (parity: the pair {-s/2, s/2}).
- Complex-linear in x  =>  complex orthogonal design of order 4 in the y's; rate 1 impossible (Hurwitz-Radon
  over C: 6 real dimensions in order 4). So W must lie in a 6-real-dimensional type: e.g. two symmetric
  blocks, or D = C^*, or skew A + ...
- TODO: brute-force which types admit a monomial design (theory/fable_tmult.py), then cost the search for
  that type at 223 / 491 / 571 / 283.

## 2026-11-10 22:00 UTC: T-multiplication idea decided (negative) [M]

Scripts: theory/fable_tmult.py (exhaustive MITM, t = 3), fable_tmult_types.py (all 164 letter-relation types at
t = 3), fable_cod4.py (monomial complex orthogonal designs of order 4), fable_tmult_sat.py / fable_tmult_scan.py
(CP-SAT, agrees with the exhaustive code at t = 3). Logs: logs/fable_fresh/tmult_*.
- t = 3: a formal multiplication (every fibre carries +-X or +-X*, any assignment) exists exactly for the 59
  types of real dimension <= 4 (Williamson s,s,s,s; "propus-like" (A, A or A*, B sym, C sym); "Golay-like"
  (A, A*, C, C*) etc.), and for NO type of real dimension 5, 6 or 8 (105 types, incl. generic, two symmetric
  blocks, one repeated block).
- t = 5 generic: INFEASIBLE (CP-SAT proof, 37 s). AB_sym at t = 5: running.
- Order-4 monomial complex ODs (T (x) W form, any T-matrices): only the all-symmetric type has one.
- Consequence: the multiplied types have 2v sign bits (or fewer), so by the usual count they stop existing near
  v ~ 160; none at 223, 283, 491, 571. Route closed. (Side remark: propus-type x Z_3 gives a GS quadruple of
  order 3v; not checked against the literature; only reaches known orders.)
- Bordered variant (identity only at rho != 1) gives n = 4v(t+1): product orders only.

## Background run (idle cores): skew 361 split torus, fresh lists, seed 2
- Launcher: /home/user/skew_work/run361.sh [threads] (copy: scripts/fable_run361.sh). Work dir
  /home/user/skew_work/p361_s2 (lists cached; kt4 resumes). Log logs/fable_fresh/p361_s2.log.
- Needs: pip install scipy numba ortools; binaries in /home/user/skew_work/bin:
  cc -O3 -march=native -fopenmp theory/{skew_enum_lowe2,skew_lg2,skew_lg,skew_enum_lowe}.c -lm; kt4.c from Kumar.
- Expectation [E]: 0.5 (Metropolis yield 0.07) to 7 (Gaussian) twin classes in the whole lists; the earlier run
  (seed 1, 10,000 of 33,600 split values) found none.

## 2026-11-10 22:05 UTC: octonion-type product GS(t) x GS(v) -> GS(tv) decided (negative) [M + P]

theory/fable_octo.py, logs/fable_fresh/octo.log.
- Z_j = (1/2) sum_i eps V_i(e tau) W_{pi_j(i)}(f g): PSD sum 4tv whenever z = composition with
  sum |z_j|^2 = |v|^2 |w|^2. Exhaustive: 2048 monomial compositions, all on the Klein Latin square; their
  conjugation-pattern sets are, for each factor, either O = {+---, +-++, ++-+, +++-} or
  E = {++++, ++--, +-+-, +--+}.
- Z_j is +-1 valued iff prod_i V_i(e_ji tau) and prod_k W_k(f_jk g) are constant for every row pattern.
  E  <=>  T-matrices with symmetric supports.  O  <=>  8 disjoint ternary sequences, 4 symmetric + 4 skew, all
  supports symmetric, zero total PAF. Both impossible for odd order > 1: at lag s the pair {-s/2, s/2} lies in
  one support and every other pair is matched by g -> -g-s, so sum PAF(s) is odd [P].
- So no product theorem of this monomial kind for two GS quadruples of odd orders > 1.

## 2026-11-10 22:10 UTC: typed multiplication closed at t = 5 as well [M]

theory/fable_tmult_h.py (hashed MITM, exact re-check; reproduces t = 3 counts 384 / 576 / 0 / 0).
- t = 5: the four minimal classes of real dimension 5, (f,s,s,s), (f,f,s,=2), (f,f,s,=0), (f,f,s,*0): 0 solutions
  each (5.0e7 pairs of block classes). Every type of dimension >= 5 specialises to one of them (an assignment
  for a type restricts to its specialisations), so none has a formal multiplication at t = 5.
  Williamson type at t = 5: 18,432.  CP-SAT on (f,f,s,s) had been UNKNOWN at 1500 s; superseded.
- Correction of an earlier line: at t = 3, 59 types have a multiplication, all of dimension <= 4; not every type
  of dimension <= 4 has one (28 of dimension 4 and 20 of dimension 3 have none).

## Large index k with a tiny subgroup S in EKLO's Theorem 7 (worked by hand; [P] where marked)

- Cell-energy identity [P]: with S = Z_3 exactly one of the 4k cells is constant; with Paley-type cells on F_5 two
  cells, on F_7 one cell, must be 'constant with a flipped zero'; larger q impossible (Theorem T flavour).
- S = Z_3, p = 4k - 1: the core is equivalent to (z, y): z = four sequences on Q over {0} u mu_6 with a single
  zero and Hermitian PAF sum 0 (a 4-block conference matrix over sixth roots of unity of order p + 1), and
  y = -z^3 (+-3 at the zero) with PAF sum 12 off 0. More generally an automorphism of order 3 of H(4p) with f
  fixed points gives a Butson weighing matrix over mu_6 of order (4p - f)/3 and weight p.
- Algebraic point: z = chi_6(det) on PG(1,p) (p = 1 mod 6; e.g. 751, 1039, 1087): z satisfies its equations
  exactly, but then y = 3I - (Paley skew conference), PAF sum 0 instead of 12, and no rescaling changes y
  (cubing kills mu_3). Fails by a constant.
- Border for large k [P for the reformulation, E for the conclusion]: the border exists iff there are two +-1
  orthogonal bases E, Z of the block-balanced space W (dim 4k - 4) with all entries of E M^T Z^T equal to
  +-4k (M the compressed core): the bases M E and Z are 'Hadamard-unbiased'. For a generic M that is
  (4k-4)^2 exact conditions against 2 (4k)^2 signs: no solution expected for large k unless M is a scalar
  matrix plus a spike (Kumar's case, which needs |S| >= about sqrt(k)). So large k with |S| <= 7 is closed.

## 2026-11-11 03:05 UTC: state after a container restart (22:1x -> 02:55, nothing ran in between)

- Report draft committed: report/fable_fresh_answer.md (update the RUN361 block and section 7 at the end).
- 361 run is now a multi-seed loop: scripts/fable_run361_loop.sh (installed copy /home/user/skew_work/loop361.sh;
  /home/user/skew_work/run361.sh starts it; flock prevents doubles). Per seed: fresh lists (about 25 min per block
  on 4 threads), then the top 2000 ranked split values. Logs logs/fable_fresh/p361_s<seed>.log and p361_loop.log.
  A hit writes cert_p361torus_*.json in the seed's work dir, a line in logs/skew/HITS.jsonl, and runs
  scripts/verify_skew.py; then publish with theory/skew_publish.py (format of pub/).
- After a restart: `/home/user/skew_work/run361.sh 4` (if /home/user/skew_work vanished: rebuild bin as noted above,
  copy scripts/fable_run361_loop.sh to /home/user/skew_work/loop361.sh and scripts/fable_run361.sh to run361.sh).
- T-matrices with three half-amicability relations: none for t = 3..11 [M] (theory/fable_tamic.py); relation set
  verified numerically.

## 2026-11-11 11:40 UTC: partial character borders (new family) -- first verified instance

(Container was down 03:55 -> 11:09; uptime windows are about 55 minutes; long runs only progress in cached pieces.)
- Framework [P, in report 6.4]: a coset-constant border of a GS core forces PSD sum = n at the characters of G
  outside a set T of characters of Q = G/S and = 4|T| on T, n = 4|G| + 4|T|.  T = {1}: Theorem 5; T = nontrivial:
  EKLO Thm 7; T = all: Kumar's all-character family.  For Q = Z_2^a every character is a +-1 vector, so ANY subset T
  is a candidate shape, and the border depends only on the 4 x 4 integer blocks Gamma_rho (rho in T) of the
  compressed core (Walsh values of the compressions on T), not on S.
- Verified [M]: theory/fable_border4.py 3 3 0,1,2,4 1 380 1 -> Q = Z_2^3, S = Z_3, T = four affinely independent
  characters (not a flat, not a flat minus a point): core by annealing, single-character border by CP-SAT (323 s),
  H(112) assembled, H H^T = 112 I.  Controls: T = nontrivial of Z_2^2 (H(60)), T = {1} (H(52)), T = three
  characters of Z_2^3 = twisted k4 with even |K| (H(108), fable_border2.py).
- Undecided (CP-SAT timeouts, not zeros): general (non single-character) borders for T = 5 characters of Z_2^3
  (b = 20) and for S = Z_3 with Q = Z_6 (k = 6, EKLO shape, non-spike compressions).
- For prime p = 3 mod 4 the new shapes need |T| = 3 mod 4 and start at Q = Z_2^4: p = 16 s + 7 (|T| = 7 not a
  3-flat minus a point) or 16 s + 11; 599 = 16*37 + 7, 907 = 16*56 + 11, 971 = 16*60 + 11, 1019 = 16*63 + 11.
  At 599 with the multiplier group of order 9 on Z_37 the Q-level forces exactly two cells of size 17 and two of
  size 19 (Parseval), and a block with one big cell cannot have 7 small Walsh values (7*12^2 > 16*15): dead for
  every T unless T lies in an affine hyperplane, and then the 8th character of the hyperplane fails.  [P, by hand]
- TODO next window: (1) T of size 7 in Z_2^4 (not flat minus point), S = Z_3: core + border, order 220;
  (2) affine classes of T for k = 16, 32 and the list of (a, |T|, s) per open prime; (3) report section.

## 2026-11-11 11:45 UTC: correction of attribution

Kumar's search/NOTES.md section 18 already states the bordered GS arrays with a character set T (any Galois-stable
T), the orthogonal-matrix form of the border problem, and "for odd q the trivial character must lie in T".  So the
framework in report 6.3/6.4/7 is his.  What this session adds: the converse (a coset-constant border must have that
shape), the corollary for skew 4v (v odd), the affine classification of T in Z_2^a and the observation that only
three classes were ever used, the block-diagonal remark, one verified instance with a non-classical T (H(112)), the
Q-level kills at 599 and 751, and the Parseval table.  Report rewritten accordingly.
- theory/fable_qlevel.py (CP-SAT Q-level gate, all T at once): UNKNOWN at 300 s on the case 599 / |U| = 9 that is
  dead by hand; too weak, a dedicated enumerator is needed.
- Next session TODO: dedicated exact Q-level gate (big-cell configurations x +-1 remainder, Walsh MITM);
  border library per Walsh-quadruple type for |T| = 7, 11 (b = 28, 44).

## 2026-11-11 16:30 UTC: a live Q-level at 599 in an unused shape  (window 4; restarts at 11:55->16:19)

- theory/fable_tcensus2.py: Z_2^4 x Z_s with H = <sigma> x U, sigma of odd order in GL(4,2) acting on Q, T a union of
  sigma-orbits of characters, non-classical; exact Q-level gate.  Partial (599, 907, 971 done; resumes per prime):
  66 dead, 65 alive, rest capped.
- **599 = 16*37 + 7, G = Z_2^4 x Z_37, sigma of order 3 without nonzero fixed points, U = <10> of order 3 on Z_37
  (|H| = 9 odd, 78 orbits per block), T = two 2-flats of characters through 1 (7 characters, e.g.
  {0,1,2,3,4,8,12}; 10 such T): 2688 Q-level solutions [M, theory/fable_q599.py, re-checked from scratch;
  logs/fable_fresh/q599_solutions.json].**  Cell sums are small (1, 5, 7, 11, 13, 17).  At every character of T
  the four Walsh values are (+-2, +-2, +-2, +-4); in 2592 of the 2688 each block's Walsh value is constant on T,
  so the compressed core acts on the border space as Gamma (x) I_7 with Gamma = GS(4,2,2,2) up to signs/order.
- Border (b = 28) for this shape: NOT found yet.  Block functions with Walsh support in T: 30 (the +-1 functions
  of the two coordinates dual to one 2-flat, or to the other).  Tensor bases E = H_4 (x) (7 orthogonal block
  functions) give 0 admissible Z-rows (expected: 7 is not a square); CP-SAT with tables was killed for memory.
  theory/fable_border5.py, fable_border6.py.  Kumar has Singer-invariant borders for the same Gamma with the
  classical T (k = 8); here a sigma-invariant search (order 3) is the analogue to try.
- If a border exists, the slice has the shape of the repository's staged slices (G' x Z_q with q = 37, classes of
  index 12, n_B = 6 orbits of G' = Z_2^4 under sigma, sGeff = 72): cost it with theory/open_screen_deficit.py after
  changing the target, then stage 2.  Fibre per Q-level solution about 2^230, 36 class equations (degree-6 real
  field) [E, rough].

## 2026-11-11 18:45 UTC: wrap-up of the idea phase (user asked to summarise and push; no new ideas until reset)

- 599 slice (report 7.1): Gamma = GS(+-2,+-2,+-2,+-4), constant on T for 2592/2688 Q-level solutions.
  Kumar's 3356 border is not single-character (blocks balanced, mostly not characters): no transplant.
  Cyclic (order 7) single-character borders: exactly none (theory/fable_border9.py).  General border undecided.
- Background, mechanical (safe to continue without new ideas):
  1. 361 skew loop: /home/user/skew_work/run361.sh 4  (resumes; log p361_loop.log, p361_s<seed>.log).
  2. t = 13 relation search: nohup nice -n 19 python3 theory/fable_tamic13.py 13 logs/fable_fresh/tamic13.jsonl &
     (58/64 chunks, 0 so far).
  3. Q-level census: nohup nice -n 10 python3 theory/fable_tcensus2.py 30 120 >> logs/fable_fresh/tcensus2.log &
     (resumes per prime; RAM: do not run together with a CP-SAT border job).
  4. CP-SAT single-character border for the 599 slice, when memory allows (kill the census first):
     nice python3 theory/fable_border7.py logs/fable_fresh/q599_solutions.json 0 37 7200 2
- After a reset, open questions for new thinking: sigma-invariant borders for the 599 slice; costing its stage 2.
- 2026-11-11 20:05 UTC check-in: container had restarted again (6th); nothing finished since 18:30. Relaunched the
  361 loop (block 3 of seed 2 regenerating), the t = 13 search (58/64) and the CP-SAT border job (no verdict yet).
