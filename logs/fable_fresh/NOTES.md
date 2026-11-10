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
