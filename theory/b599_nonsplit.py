r"""Non-split c = 6 slices of 599: G = Z_2 x Z_27 x Z_11, H = 3Z_27 x Z_11 (= Z_9 x Z_11), Q = G/H = Z_2 x Z_3 [P, E].

G/H: (e, x, z) -> (e, x mod 3), a homomorphism onto Z_2 x Z_3 with kernel H (index 6).  The Q-level, the border
classes (the forced corner depends only on the Q-compression and the frames over Q; report/alpoge_design_answer.md
section 3) and the stage-1 target T1 on G' = Z_2 x Z_27 (2176 at 0, -220 on H' \ 0 = 3Z_27 \ 0 in layer 0, +44
elsewhere) are those of the split lead.  The slice: Paley(11) on Z_11 (squares), U = a multiplier group on Z_27:
  x10 (order 3, = 1 mod 9: fixes H' pointwise, permutes each other coset in 3 orbits of 3);  -1 (moves coset 1 to 2);
  x4 (order 9, = 1 mod 3: keeps cosets).
The conditional stage-1 estimator (theory/b599_cond.py with an explicit coset map) and the stage-2 typical-y budget
(theory/ab_typical_levels.level_bits, ab_levels.sd_freedom) are run on these slices.
usage: python3 theory/b599_nonsplit.py --u 10 [--stage2] [--out FILE]"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ab_lib import Slice, AG, primitive_root
import b599_cond as BC


def build(us):
    mods = (2, 27, 11); G = AG(mods)
    Hmask = (G.co[:, 0] == 0) & (G.co[:, 1] % 3 == 0)
    g11 = pow(primitive_root(11), 2, 11)
    gens = [(1, 1, g11)] + [(1, u, 1) for u in us]
    sl = Slice(mods, 599, 5, gens, 2, Hmask)
    Gp = sl.Gp                                   # Z_2 x Z_27
    qc = 3 * Gp.co[:, 0] + (Gp.co[:, 1] % 3)
    cfg = dict(v=599, c=6, b=5, h=99, m=594, Q=(2, 3), H1=('3Z27',), q=11, s=2, f=5, mods=mods, nQ=2,
               U=[('mul27', u) for u in us])
    return cfg, sl, qc


if __name__ == "__main__":
    a = sys.argv[1:]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    us = [int(x) for x in opt('--u', '10').split(',')]
    cfg, sl, qc = build(us)
    print(json.dumps(dict(U=us, orbits=int(sl.r), Gp_orbits=len(sl.breps), T1=sorted(set(sl.T1.tolist())))), flush=True)
    C = BC.Cond(cfg, sl=sl, Qmods=(2, 3), qc=qc)
    outp = opt('--out', ''); fo = open(outp, 'a') if outp else None
    classes = C.qlevel_quads()
    tot = 0.0; nq = 0; t0 = time.time()
    for ci, (key, us_) in enumerate(sorted(classes.items())):
        for zm in C.rows_for(us_):
            Z = np.array([[zm[k] + 2 * us_[k][q] for q in range(C.c)] for k in range(4)], np.int64)
            r = C.budget(Z, N=int(opt('--N', '2000')))
            nq += 1
            if r.get('log2E') is not None and np.isfinite(r['log2E']): tot += 2 ** r['log2E']
            rec = dict(U=us, **{'class': ci}, shapes=[list(u) for u in us_], zmin=list(zm), rows=[int(x) for x in Z.sum(1)], **r)
            if fo: fo.write(json.dumps(rec) + "\n"); fo.flush()
    print(json.dumps(dict(stage="done", U=us, classes=len(classes), zquads=nq, sum_log2E=round(math.log2(tot), 2) if tot else None,
                          secs=round(time.time() - t0, 1))), flush=True)


def zquads_x10():
    """Q-level compressions for the x10 slice [P]: the coset sums off H' are odd multiples of 3 (three orbits of
    weight 3), so z_k = c_k + 2 a_k delta_(0,0) + 2 b_k delta_(1,0) with c_k an odd multiple of 3; the Q-level reads
    sum_k (a_k^2 + b_k^2) = 5 and sum_k a_k b_k = 0 (u_k^(eps, w^j) = a_k + eps b_k for every psi != 1).
    Up to block order, block negation and the swap a <-> b (translation by (1, 0))."""
    import itertools
    pairs = [(a, b) for a in range(-2, 3) for b in range(-2, 3)]
    out = set()
    for P in itertools.combinations_with_replacement(pairs, 4):
        if sum(a * a + b * b for a, b in P) != 5 or sum(a * b for a, b in P) != 0: continue
        out.add(tuple(sorted(P)))
    quads = []
    for P in sorted(out):
        for cs in itertools.product([x for x in range(-21, 22, 6)], repeat=4):     # odd multiples of 3
            rows = [6 * cs[k] + 2 * P[k][0] + 2 * P[k][1] for k in range(4)]
            if sum(r * r for r in rows) != 2396: continue
            Z = np.array([[cs[k] + (2 * P[k][0] if q == 0 else 0) + (2 * P[k][1] if q == 3 else 0) for q in range(6)]
                          for k in range(4)], np.int64)
            quads.append((P, cs, Z))
    return quads


def run_x10(N=2000, outp=''):
    cfg, sl, qc = build([10])
    C = BC.Cond(cfg, sl=sl, Qmods=(2, 3), qc=qc)
    fo = open(outp, 'a') if outp else None
    tot = 0.0; recs = []
    for P, cs, Z in zquads_x10():
        r = C.budget(Z, N=N)
        if r.get('log2E') is not None and np.isfinite(r['log2E']): tot += 2 ** r['log2E']
        rec = dict(U=[10], ab=[list(p) for p in P], c=list(cs), rows=[int(x) for x in Z.sum(1)], **r)
        recs.append(rec)
        if fo: fo.write(json.dumps(rec) + "\n"); fo.flush()
    return tot, recs
