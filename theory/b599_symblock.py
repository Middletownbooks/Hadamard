r"""Stage 1 of the c = 6 family (theory/b599_fam.py) restricted to symmetric blocks [M].

Gamma = (Q-affine maps q -> a q + s, a in Aut(Q) = {+-1 on Z_3}, s in Q) x (Z_n^* / +-1) acts on block vectors
(cells (q, j), j a representative of {j, -j} in Z_n).  Every element maps stage-1 solutions to stage-1 solutions when
applied to one block (it is an automorphism-and-translation of G' preserving H', commuting with -1 on Z_n), and
it maps the block's Q-compression z to z o (affine map).  A block is *symmetric* if a nontrivial element of Gamma
fixes it; then it is fixed by an element of prime order, whose Q-part fixes z.
At the same-shape rung v = 215 the low-count compressions have 93% of their exact solutions with all four blocks
symmetric (report/b599_answer.md section 3), and the conditional Gaussian (generic blocks) undercounts them.
Search (exact for the sub-family): per block, the union over prime-order g fixing z_k of Fix(g) with compression
z_k and energy <= its cap (product over the g-orbits of cosets, exact sums), then the exact join
theory/p647_join.join4 on the stage-1 class signatures; every hit re-checked by the full G' PAF.
usage: python3 theory/b599_symblock.py --n 9 --q 11 [--classes ...] [--out FILE]"""
import sys, os, json, math, itertools, time
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b599_fam import Fam
import b599_flevels as FL
import p647_join as PJ


def gamma(F):
    """list of (label, cell permutation P (image index of each cell), coset map cm, order) for every element"""
    J, n = F.J, F.n
    mus = []
    seen = set()
    for m in range(1, n):
        if math.gcd(m, n) != 1: continue
        key = frozenset([m, n - m])
        if key in seen: continue
        seen.add(key); mus.append(m)
    out = []
    for a in (1, -1):
        for s in range(6):
            se, st = divmod(s, 3)
            cm = [3 * ((e + se) % 2) + ((a * t) + st) % 3 for e in range(2) for t in range(3)]
            for mu in mus:
                jm = [min((mu * j) % n, (-mu * j) % n) for j in range(J)]
                P = np.zeros(6 * J, np.int64)       # P[image cell] = source cell: y'(image) = y(source)
                for q in range(6):
                    for j in range(J): P[J * cm[q] + jm[j]] = J * q + j
                # order
                Q_ = np.arange(6 * J); o = 0
                while True:
                    Q_ = Q_[P]; o += 1
                    if (Q_ == np.arange(6 * J)).all(): break
                out.append(((a, s, mu), P, cm, o))
    return out


def fixed_blocks(F, z, P, cm, Emax):
    """all block vectors y with y[P] == y (g-invariant), Q-compression z and energy <= Emax"""
    J = F.J; L = F.LET; W = F.W
    ncell = 6 * J
    # cell orbits of g
    orb = -np.ones(ncell, np.int64); reps = []
    for x in range(ncell):
        if orb[x] >= 0: continue
        cyc = [x]; y = P[x]
        while y != x: cyc.append(y); y = P[y]
        for c in cyc: orb[c] = len(reps)
        reps.append(cyc)
    # coset orbits under cm (cm[q] = image coset)
    corb = -np.ones(6, np.int64); creps = []
    for q in range(6):
        if corb[q] >= 0: continue
        cyc = [q]; y = cm[q]
        while y != q: cyc.append(y); y = cm[y]
        for c in cyc: corb[c] = len(creps)
        creps.append(cyc)
    if any(len({int(z[q]) for q in cyc}) > 1 for cyc in creps): return np.zeros((0, ncell), np.int64)
    # per coset orbit: the cell orbits contained in it; enumerate letters, check the sum of one coset of the orbit
    parts = []
    for cyc in creps:
        cells = [J * q + j for q in cyc for j in range(J)]
        os_ = sorted({int(orb[c]) for c in cells})
        q0 = cyc[0]
        wq0 = np.array([sum(int(W[c]) for c in reps[o] if c // J == q0) for o in os_], np.int64)   # weight in coset q0
        wt = np.array([sum(int(W[c]) for c in reps[o]) for o in os_], np.int64)                      # total weight
        T = np.array(list(itertools.product(L.tolist(), repeat=len(os_))), np.int64)
        m = (T @ wq0) == int(z[q0])
        T = T[m]; e = (T * T) @ wt
        parts.append((os_, T, e))
    mins = [int(e.min()) if len(e) else 10 ** 9 for _, _, e in parts]
    if max(mins) >= 10 ** 9: return np.zeros((0, ncell), np.int64)
    rem = [sum(mins[i + 1:]) for i in range(len(parts))]
    vals = np.zeros((1, len(reps)), np.int64); E = np.zeros(1, np.int64)
    for i, (os_, T, e) in enumerate(parts):
        ne = E[:, None] + e[None, :]
        a_, b_ = np.nonzero(ne + rem[i] <= Emax)
        nv = vals[a_].copy(); nv[:, os_] = T[b_]
        vals = nv; E = ne[a_, b_]
        if len(vals) > 3 * 10 ** 7: raise RuntimeError("fixed list too large")
    return vals[:, orb]


def fixdim(P):
    seen = set(); d = 0
    for x in range(len(P)):
        if x in seen: continue
        d += 1; y = x
        while y not in seen: seen.add(y); y = P[y]
    return d


def sym_list(F, z, G, Emax, types=None):
    """types: optional set of (order, fixed dimension) pairs of the elements used"""
    out = []
    for lab, P, cm, o in G:
        if lab == (1, 0, 1): continue
        if o not in (2, 3, 5, 7): continue
        if types is not None and (o, fixdim(P)) not in types: continue
        if any(z[cm[q]] != z[q] for q in range(6)): continue
        V = fixed_blocks(F, z, P, cm, Emax)
        if len(V): out.append(V)
    if not out: return np.zeros((0, 6 * F.J), np.int64)
    return np.unique(np.concatenate(out), axis=0)


def solve_z(F, z, G=None, types=None):
    G = gamma(F) if G is None else G
    T0 = int(F.T1SIG[0])
    Ls = [sym_list(F, z[k], G, T0, types) for k in range(4)]
    if min(len(x) for x in Ls) == 0: return dict(sizes=[len(x) for x in Ls], solutions=0), []
    mins = [int(F.energy(x).min()) for x in Ls]
    for k in range(4):
        cap = T0 - (sum(mins) - mins[k]); Ls[k] = Ls[k][F.energy(Ls[k]) <= cap]
    sizes = [len(x) for x in Ls]
    if min(sizes) == 0: return dict(sizes=sizes, solutions=0), []
    S = [F.sig(x) for x in Ls]
    big = max(sizes)
    P = 4093 if big > 2 * 10 ** 5 else (509 if big > 2 * 10 ** 4 else 61)
    sols, st = PJ.join4(S, F.T1SIG.astype(np.int64), P=P)
    Ys = []
    for s_ in sols:
        Y = np.array([Ls[k][s_[k]] for k in range(4)])
        assert F.check_stage1(Y) and (F.qcomp(Y) == np.asarray(z)).all()
        Ys.append(Y)
    return dict(sizes=sizes, solutions=len(Ys), work_log2=st.get('work_log2'), join_secs=st.get('secs')), Ys


if __name__ == "__main__":
    a = sys.argv[1:]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    n, q = int(opt('--n', '9')), int(opt('--q', '11'))
    outp = opt('--out', '')
    cls = None if '--classes' not in a else [int(x) for x in opt('--classes', '').split(',')]
    part, nparts = (int(x) for x in opt('--part', '0/1').split('/'))
    types = None if '--types' not in a else {tuple(int(y) for y in t.split(':')) for t in opt('--types', '').split(',')}
    F = Fam(n, q); G = gamma(F)
    done = set()
    if outp and os.path.exists(outp):
        for l in open(outp):
            try: d = json.loads(l); done.add((d['class'], tuple(d['c'])))
            except Exception: pass
    fo = open(outp, 'a') if outp else None
    tot = 0
    for i, (ci, cl, c, z) in enumerate(FL.class_quads(F)):
        if i % nparts != part: continue
        if cls is not None and ci not in cls: continue
        if (ci, tuple(c)) in done: continue
        t0 = time.time()
        try:
            res, Ys = solve_z(F, z, G, types)
        except (RuntimeError, MemoryError) as e:
            res, Ys = dict(error=repr(e)[:200]), []
        tot += len(Ys)
        rec = dict(n=n, q=q, v=F.v, types=sorted(types) if types else 'all', **{'class': ci}, c=list(c), rows=[int(x) for x in z.sum(1)], **res, secs=round(time.time() - t0, 1))
        if Ys: rec['examples'] = [Y.tolist() for Y in Ys[:50]]
        print(json.dumps({k: v for k, v in rec.items() if k != 'examples'}), flush=True)
        if fo: fo.write(json.dumps(rec) + "\n"); fo.flush()
    print(json.dumps(dict(stage="done", total_solutions=tot)))
