r"""Every core of a small c = 4 slice, without staging: the character solver of theory/ab_s1_chars.py run on the whole
group G = Q x H'' (H'' = H' x Z_q, the slice group acting on H''), with letters +-1 on every G-orbit [M].

Used to calibrate the stage-2 Gaussian (theory/ab_s2_s4.py) against an exact core count on a small slice of the
same family.  The admissible row quadruples are all sorted quads of non-negative integers of the parity of |G| with
sum of squares sum_g T(g) (= 4v).  Output: raw int8 rows (4 x r orbit signs per core, canonical up to block order
and block sign), checkpointed, every core verified (sum_k PAF_G = T).
usage: python3 theory/ab_full_chars.py CFG.json OUT.cbin [--plant SEED] [--maxB N]"""
import sys, os, math, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_s1sd as SD, ab_s1_chars as AC


class Unstaged:
    """the attributes theory/ab_s1_chars.solve reads, for the whole group G (G' := G, letters +-1)"""
    def __init__(s, sl):
        s.Gp = sl.G; s.borb = np.asarray(sl.orb); s.breps = list(sl.reps); s.Sp = list(sl.S)
        s.T1 = np.asarray(sl.T, np.int64); s.letters = [[-1, 1] for _ in range(sl.r)]


def row_quads(G, T):
    R = int(np.asarray(T).sum()); par = G.n % 2; out = []
    vals = [x for x in range(par, math.isqrt(R) + 1, 2)]
    for q in itertools.combinations_with_replacement(vals, 4):
        if sum(x * x for x in q) == R: out.append(q)
    return out


def main():
    a = sys.argv[1:]
    cfg = SD.cfg_from_json(open(a[0]).read()); sl = W.build(cfg); U = Unstaged(sl)
    AC.solve(cfg, U, a[1], a, rq=row_quads(sl.G, sl.T))


if __name__ == "__main__":
    main()
