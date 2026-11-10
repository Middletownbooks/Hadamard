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
