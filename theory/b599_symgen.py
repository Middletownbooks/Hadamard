r"""Symmetric-block stage 1 for any slice handled by theory/b599_cond.py (U trivial on Q) [M].

Gamma = (Q-affine maps q -> A q + s, A in Aut(Q)) x (automorphisms B of H' commuting with U), acting on G' = Q x H'
by (q, h) -> (A q + s, B h).  Each element is an automorphism of G' preserving H' composed with a translation by Q,
so applied to one block it maps stage-1 solutions to stage-1 solutions (T1 depends only on membership in H' and
on 0; the fibre letter sets are preserved because B commutes with U), and it maps the block's Q-compression z to
z o (A q + s).  A block is symmetric if a nontrivial element fixes it; then an element of prime order does.
Search (exact for the family where all four blocks are symmetric): per block, the union over prime-order elements
whose Q-part fixes z_k of their fixed blocks with compression z_k and energy <= cap; then the exact join
theory/p647_join.join4 on the class signatures; every hit verified by the full G' PAF and the letters.
Aut(H') supported: H' = Z_n (units) or Z_3^2 (GL(2, 3)); Aut(Q) for Q = Z_2 x Z_3: {+-1 on Z_3}.
usage: python3 theory/b599_symgen.py 'CFGJSON' [--types 2,3] [--classes ...] [--part i/W] [--out FILE]"""
import sys, os, json, math, itertools, time
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ab_s1sd import cfg_from_json
from ab_lib import AG
import b599_cond as BC
import p647_join as PJ


def haut(mods_h, U_h):
    """automorphisms of H' (as integer matrices acting on coordinate vectors) commuting with every u in U_h"""
    r = len(mods_h)
    if r == 1:
        n = mods_h[0]
        mats = [np.array([[u]]) for u in range(1, n) if math.gcd(u, n) == 1]
    elif tuple(mods_h) == (3, 3):
        mats = []
        for M in itertools.product(range(3), repeat=4):
            A = np.array(M).reshape(2, 2)
            if round(np.linalg.det(A)) % 3 != 0: mats.append(A)
    else:
        raise NotImplementedError(mods_h)
    out = []
    for A in mats:
        ok = True
        for u in U_h:
            D = np.diag(u)
            if not (((A @ D) - (D @ A)) % np.array(mods_h)[:, None] == 0).all(): ok = False; break
        if ok: out.append(A)
    return out


def gamma(C):
    sl = C.sl; Gp = sl.Gp; nQ = C.cfg['nQ']
    mods = Gp.mods; mh = mods[nQ:]
    U_h = [tuple(u[nQ:]) for u in sl.Sp]
    Hs = haut(mh, U_h)
    Qg = C.Q
    # Aut(Q) for Q = Z_2 x Z_3 (or any): unit scalars per axis
    qaut = [np.array(us) for us in itertools.product(*[[u for u in range(1, m) if math.gcd(u, m) == 1] or [0] for m in Qg.mods])]
    els = []
    co = Gp.co
    for qa in qaut:
        for s in range(Qg.n):
            for B in Hs:
                qimg = (co[:, :nQ] * qa + Qg.co[s]) % np.array(Qg.mods)
                himg = (co[:, nQ:] @ B.T) % np.array(mh)
                img = Gp.idx(np.concatenate([qimg, himg], 1))          # element g -> image
                # orbit permutation: P[orbit of image of rep] = orbit of rep  (y'(img) = y(src))
                P = np.zeros(C.L1, np.int64)
                for j, bb in enumerate(sl.breps): P[sl.borb[img[bb]]] = j
                # well-defined on orbits?
                okp = True
                for g in range(0, Gp.n, max(1, Gp.n // 64)):
                    if P[sl.borb[img[g]]] != sl.borb[g]: okp = False; break
                if not okp: continue
                cm = np.zeros(Qg.n, np.int64)
                for q in range(Qg.n): cm[q] = Qg.idx((Qg.co[q] * qa + Qg.co[s]) % np.array(Qg.mods))
                # order
                Q_ = np.arange(C.L1); o = 0
                while True:
                    Q_ = Q_[P]; o += 1
                    if (Q_ == np.arange(C.L1)).all() or o > 100: break
                els.append(((tuple(qa.tolist()), s, tuple(B.ravel().tolist())), P, cm, o))
    # deduplicate identical permutations (with identical coset maps)
    seen = {}; out = []
    for lab, P, cm, o in els:
        key = (tuple(P.tolist()), tuple(cm.tolist()))
        if key in seen: continue
        seen[key] = 1; out.append((lab, P, cm, o))
    return out


def fixed_blocks(C, z, P, cm, Emax):
    """all block vectors y (per G'-orbit values) with y[P] == y, Q-compression z, energy <= Emax"""
    L1 = C.L1; W = C.w
    orb = -np.ones(L1, np.int64); reps = []
    for x in range(L1):
        if orb[x] >= 0: continue
        cyc = [x]; y = P[x]
        while y != x: cyc.append(int(y)); y = P[y]
        for c in cyc: orb[c] = len(reps)
        reps.append(cyc)
    corb = -np.ones(C.c, np.int64); creps = []
    for q in range(C.c):
        if corb[q] >= 0: continue
        cyc = [q]; y = cm[q]
        while y != q: cyc.append(int(y)); y = cm[y]
        for c in cyc: corb[c] = len(creps)
        creps.append(cyc)
    if any(len({int(z[q]) for q in cyc}) > 1 for cyc in creps): return np.zeros((0, L1), np.int64)
    parts = []
    for cyc in creps:
        cells = [j for j in range(L1) if int(C.oq[j]) in cyc]
        os_ = sorted({int(orb[c]) for c in cells}); q0 = cyc[0]
        wq0 = np.array([sum(int(W[c]) for c in reps[o] if int(C.oq[c]) == q0) for o in os_], np.int64)
        wt = np.array([sum(int(W[c]) for c in reps[o]) for o in os_], np.int64)
        letsets = [sorted(set.intersection(*[set(C.lets[c].tolist()) for c in reps[o]])) for o in os_]
        T = np.array(list(itertools.product(*letsets)), np.int64)
        if len(T) == 0: return np.zeros((0, L1), np.int64)
        T = T[(T @ wq0) == int(z[q0])]; e = (T * T) @ wt
        parts.append((os_, T, e))
    mins = [int(e.min()) if len(e) else 10 ** 9 for _, _, e in parts]
    if max(mins) >= 10 ** 9: return np.zeros((0, L1), np.int8)
    rem = [sum(mins[i + 1:]) for i in range(len(parts))]
    vals = np.zeros((1, len(reps)), np.int8); E = np.zeros(1, np.int64)
    for i, (os_, T, e) in enumerate(parts):
        o_ = np.argsort(e, kind='stable'); T = T[o_]; e = e[o_]          # sorted energies: prefix per row
        nv_l = []; ne_l = []; tot = 0
        for c0 in range(0, len(vals), 4096):
            Ec = E[c0:c0 + 4096]
            lim = Emax - rem[i] - Ec                                       # allowed pattern energy per row
            cnt = np.searchsorted(e, lim, side='right')
            if cnt.sum() == 0: continue
            rows = np.repeat(np.arange(len(Ec)), cnt)
            cols = np.concatenate([np.arange(k) for k in cnt]) if len(cnt) else np.zeros(0, np.int64)
            nv = vals[c0:c0 + 4096][rows].copy(); nv[:, os_] = T[cols].astype(np.int8)
            nv_l.append(nv); ne_l.append(Ec[rows] + e[cols]); tot += len(rows)
            if tot > 2 * 10 ** 7: raise RuntimeError("fixed list too large")
        if not nv_l: return np.zeros((0, L1), np.int8)
        vals = np.concatenate(nv_l); E = np.concatenate(ne_l)
    return vals[:, orb]


def count_fixed(C, z, P, cm, Emax):
    """number of blocks fixed by (P, cm) with compression z and energy <= Emax (DP, nothing materialised)"""
    L1 = C.L1; W = C.w
    orb = -np.ones(L1, np.int64); reps = []
    for x in range(L1):
        if orb[x] >= 0: continue
        cyc = [x]; y = P[x]
        while y != x: cyc.append(int(y)); y = P[y]
        for c in cyc: orb[c] = len(reps)
        reps.append(cyc)
    corb = -np.ones(C.c, np.int64); creps = []
    for q in range(C.c):
        if corb[q] >= 0: continue
        cyc = [q]; y = cm[q]
        while y != q: cyc.append(int(y)); y = cm[y]
        for c in cyc: corb[c] = len(creps)
        creps.append(cyc)
    if any(len({int(z[q]) for q in cyc}) > 1 for cyc in creps): return 0
    tot = {0: 1}
    for cyc in creps:
        cells = [j for j in range(L1) if int(C.oq[j]) in cyc]
        os_ = sorted({int(orb[c]) for c in cells}); q0 = cyc[0]
        wq0 = [sum(int(W[c]) for c in reps[o] if int(C.oq[c]) == q0) for o in os_]
        wt = [sum(int(W[c]) for c in reps[o]) for o in os_]
        letsets = [sorted(set.intersection(*[set(C.lets[c].tolist()) for c in reps[o]])) for o in os_]
        st = {(0, 0): 1}
        for i in range(len(os_)):
            nd = defaultdict(int)
            for (sm, e), m in st.items():
                for y in letsets[i]:
                    e2 = e + wt[i] * y * y
                    if e2 <= Emax: nd[(sm + wq0[i] * y, e2)] += m
            st = nd
        per = defaultdict(int)
        for (sm, e), m in st.items():
            if sm == int(z[q0]): per[e] += m
        nt = defaultdict(int)
        for E1, m1 in tot.items():
            for E2, m2 in per.items():
                if E1 + E2 <= Emax: nt[E1 + E2] += m1 * m2
        tot = nt
    return sum(tot.values())


SKIPPED = []


def sym_list(C, z, G, Emax, orders=(2, 3), cap=10 ** 6):
    out = []
    for lab, P, cm, o in G:
        if (P == np.arange(C.L1)).all(): continue
        if o not in orders: continue
        if any(z[cm[q]] != z[q] for q in range(C.c)): continue
        n_ = count_fixed(C, z, P, cm, Emax)
        if n_ == 0: continue
        if n_ > cap:
            SKIPPED.append((lab, n_)); continue
        V = fixed_blocks(C, z, P, cm, Emax)
        if len(V): out.append(V)
    if not out: return np.zeros((0, C.L1), np.int8)
    return np.unique(np.concatenate(out), axis=0)


def check_y(C, Y):
    sl = C.sl
    Yf = np.asarray(Y)[:, sl.borb]
    if not (sl.Gp.paf(Yf).sum(0) == sl.T1).all(): return False
    return all(int(Y[k][j]) in set(C.lets[j].tolist()) for k in range(4) for j in range(C.L1))


def solve_z(C, Z, G, orders=(2, 3), maxlist=3 * 10 ** 6, cap=10 ** 6):
    T0 = C.T0
    del SKIPPED[:]
    Ls = [sym_list(C, Z[k], G, T0, orders, cap) for k in range(4)]
    skipped = len(SKIPPED); skipped_max = max((n_ for _, n_ in SKIPPED), default=0)
    if min(len(x) for x in Ls) == 0: return dict(sizes=[len(x) for x in Ls], solutions=0), []
    En = [((x.astype(np.int64)) ** 2 * C.w).sum(1) for x in Ls]
    mins = [int(e_.min()) for e_ in En]
    for k in range(4):
        cap_ = T0 - (sum(mins) - mins[k]); Ls[k] = Ls[k][En[k] <= cap_]
    sizes = [len(x) for x in Ls]
    if min(sizes) == 0: return dict(sizes=sizes, solutions=0), []
    if max(sizes) > maxlist: return dict(sizes=sizes, skipped=f"list > {maxlist}"), []
    S = [C.sig(x.astype(np.int64)) for x in Ls]
    big = max(sizes)
    P = 4093 if big > 2 * 10 ** 5 else (509 if big > 2 * 10 ** 4 else 61)
    sols, st = PJ.join4(S, C.T1sig, P=P)
    Ys = []
    for s_ in sols:
        Y = np.array([Ls[k][s_[k]] for k in range(4)], np.int64)
        assert check_y(C, Y)
        Ys.append(Y)
    return dict(sizes=sizes, solutions=len(Ys), work_log2=st.get('work_log2'), join_secs=st.get('secs'),
                skipped_elements=skipped, skipped_max=skipped_max), Ys


if __name__ == "__main__":
    a = sys.argv[1:]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    cfg = cfg_from_json(a[0]); C = BC.Cond(cfg); G = gamma(C)
    orders = tuple(int(x) for x in opt('--types', '2,3').split(','))
    outp = opt('--out', ''); part, nparts = (int(x) for x in opt('--part', '0/1').split('/'))
    maxlist = int(opt('--maxlist', str(3 * 10 ** 6)))
    cls = None if '--classes' not in a else [int(x) for x in opt('--classes', '').split(',')]
    print(json.dumps(dict(stage="gamma", elements=len(G), orders=sorted({o for *_, o in G}))), flush=True)
    classes = C.qlevel_quads()
    todo = []
    for ci, (key, us) in enumerate(sorted(classes.items())):
        if cls is not None and ci not in cls: continue
        for zm in C.rows_for(us):
            todo.append((ci, us, zm))
    done = set()
    if outp and os.path.exists(outp):
        for l in open(outp):
            try: d = json.loads(l); done.add((d['class'], tuple(d['zmin'])))
            except Exception: pass
    fo = open(outp, 'a') if outp else None
    tot = 0
    for i, (ci, us, zm) in enumerate(todo):
        if i % nparts != part or (ci, tuple(zm)) in done: continue
        Z = np.array([[zm[k] + 2 * us[k][q] for q in range(C.c)] for k in range(4)], np.int64)
        t0 = time.time()
        try: res, Ys = solve_z(C, Z, G, orders, maxlist)
        except (RuntimeError, MemoryError) as e: res, Ys = dict(error=repr(e)[:200]), []
        tot += len(Ys)
        rec = dict(cfg=dict(v=cfg['v'], H1=list(cfg['H1']), U=[list(u) for u in cfg['U']]), **{'class': ci},
                   shapes=[list(u) for u in us], zmin=list(zm), rows=[int(x) for x in Z.sum(1)], orders=list(orders),
                   **res, secs=round(time.time() - t0, 1))
        if Ys: rec['examples'] = [Y.tolist() for Y in Ys[:50]]
        print(json.dumps({k: v for k, v in rec.items() if k != 'examples'}), flush=True)
        if fo: fo.write(json.dumps(rec) + "\n"); fo.flush()
    print(json.dumps(dict(stage="done", total_solutions=tot)))
