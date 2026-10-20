r"""Controls for the bordered-GS chain (theory/bt_engine, bt_chain, bt_border) against Kumar's published structures
[M].  Reads Kumar's certificates from a local clone (default /home/user/vatskishan03/hadamard, Apache-2.0).

1. 3964 (G = Z_2^2 x Z_247, S of order 27): 46 orbits, 23 equations; Q-level at Kumar's rows (-50, -26, -2, 28);
   the Z_19-level join gives 304 solutions; the planted core (Kumar's) lies at a level-A solution with 160 level-B
   combinations; the final join of that combination returns the planted core (and others); border and assembly.
2. 1916 (G = Z_2^2 x Z_119): the split structure of the prompt (blocks 1, 4 under <(alpha^2, 16), (1, 13)>,
   blocks 2, 3 under <(1, 13)>) has 16 Q-level solutions at rows (6, 6, 38, 20); the published core is, up to a
   translation of blocks 2 and 3 by elements of Z_2^2, invariant under the order-12 group in every block.
3. 2636 (G = Z_4^2 x Z_41, non-split K, S of order 15): 54 orbits, 27 equations, 2954 compressed vectors at the
   Z_4^2 level; border solved for Kumar's core with our frames.
usage: python3 theory/bt_controls.py [--kumar DIR] [--skip-final]
"""
import sys, os, json, time, argparse, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_engine as E, bt_chain as C, bt_border as BB, bt_gates as BG

S = lambda s: np.array([1 if c == '+' else -1 for c in s], np.int64)


def log(**k):
    print(json.dumps(k), flush=True)


def c3964(kdir, final=True):
    st = E.control_3964()
    log(control='3964', r=st.r, e=st.e)
    cert = json.load(open(os.path.join(kdir, 'hadamard/3964/certificate.json')))
    xs = [S(s) for s in cert['sequences']]
    assert st.check_core(xs)
    eps = [np.array([xs[i][st.reps[i][O]] for O in range(st.r[i])], np.int8) for i in range(4)]
    subs = st.invariant_subgroups()
    N19 = [N for N in subs if len(N) == 19][0]; N13 = [N for N in subs if len(N) == 13][0]
    N4 = [N for N in subs if len(N) == 4][0]
    qQ, qA, qB, qC = (st.quotient(N) for N in (st.Kel, N19, N13, N4))
    chQ, chA, chB = (C.char_matrix(q) for q in (qQ, qA, qB))
    rows = [int(x.sum()) for x in xs]
    LQ = [C.q_level_list(qQ, i, rows[i], chQ) for i in range(4)]
    PQ = [C.paf_reduced(qQ, i, LQ[i][0]) for i in range(4)]
    sols, no = C.join4(PQ, qQ.cls_target)
    s = sols[0]
    log(level='Q', rows=rows, solutions=int(no), z=[LQ[i][0][s[i]].tolist() for i in range(4)])
    LA = [C.refine_list(qQ, LQ[i][0][s[i]], qA, i, chars=chA) for i in range(4)]
    PA = [C.paf_reduced(qA, i, LA[i][0]) for i in range(4)]
    solA, noA = C.join4(PA, qA.cls_target)
    yA = [eps[i].astype(np.int64) @ qA.contrib(i) for i in range(4)]
    yB = [eps[i].astype(np.int64) @ qB.contrib(i) for i in range(4)]
    pos = [int(np.flatnonzero((LA[i][0] == yA[i]).all(1))[0]) for i in range(4)]
    log(level='Z19 (G/Z19 = Z2^2 x Z13)', solutions=int(noA), kumar_reports=304,
        planted_at=[k for k, a in enumerate(solA) if list(a) == pos])
    LB = []
    for i in range(4):
        Y, cnt = C.fibre_dp(st, i, qA, yA[i], [], [], qB)
        en = np.abs(Y[:, qB.borb[i]].astype(float) @ chB[0].T) ** 2
        ok = np.all(en <= chB[1][None, :] + 1e-6, axis=1)
        LB.append((Y[ok], cnt[ok]))
    PB = [C.paf_reduced(qB, i, LB[i][0]) for i in range(4)]
    solB, noB = C.join4(PB, qB.cls_target)
    posB = [int(np.flatnonzero((LB[i][0] == yB[i]).all(1))[0]) for i in range(4)]
    log(level='Z13 inside the planted level-A fibre', lists=[len(l[0]) for l in LB], combinations=int(noB),
        kumar_reports=160, planted_at=[k for k, b in enumerate(solB) if list(b) == posB])
    if final:
        t = time.time()
        EPS = [C.expand_fibre(st, i, qA, yA[i], [qB], [yB[i]]) for i in range(4)]
        PHI = [C.char_orbit_matrix(st, i) for i in range(4)]
        EPS = [EPS[i][C.char_filter(st, i, EPS[i], PHI[i])] for i in range(4)]
        fs = C.final_join(st, EPS, qC)
        found = any(all((EPS[i][f[i]] == eps[i]).all() for i in range(4)) for f in fs)
        cores = [[st.expand(i, EPS[i][f[i]]) for i in range(4)] for f in fs]
        assert all(st.check_core(c) for c in cores)
        log(level='final join of the planted combination', fibres=[len(e) for e in EPS], cores=len(fs),
            planted_found=found, secs=round(time.time() - t, 1))
        nb = 0
        for c in cores:
            r = BB.solve_border(st, c)
            if r is not None and BB.is_hadamard(BB.assemble(st, c, *r[:3])):
                nb += 1
        log(level='border + assembly of every core', bordered_and_hadamard=nb, of=len(cores))


def c1916():
    A2 = np.array([[1, 1], [1, 0]])
    g12a = E.block_diag(A2, [[16]]); g12b = E.block_diag(np.eye(2, dtype=int), [[13]])
    mods = (2, 2, 119)
    split = E.BStruct(mods, [(0, 0, 1)], 'k', [[g12a, g12b], [g12b], [g12b], [g12a, g12b]], name='1916-split')
    log(control='1916 split', r=split.r, e=split.e)
    g = BG.Gates(split)
    rows = (6, 6, 38, 20)
    L = [g.q_list(i, rows[i]) for i in range(4)]
    P = [C.paf_reduced(g.qQ, i, L[i][0]) for i in range(4)]
    sol, no = C.join4(P, g.qQ.cls_target)
    log(level='Q (split structure)', rows=rows, solutions=int(no), prompt_reports=16)
    # published core: affine symmetries
    lines = open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data/core_v479_b3.txt')).read().split('\n')
    seqs = {}
    for Ln in lines[1:]:
        if '=' in Ln:
            k, v = Ln.split('=', 1); seqs[k] = S(v.strip())
    G = E.Grp(mods)
    common = E.BStruct(mods, [(0, 0, 1)], 'k', [g12a, g12b], name='1916-common')
    xs = []
    for k in 'ABCD':
        x = seqs[k]
        for c in range(4):                      # translations by Z_2^2 elements
            tr = G.idx(G.co + np.array([c // 2, c % 2, 0]))
            y = x[tr]
            if all((y[P] == y).all() for P in common.bperms[0]):
                xs.append(y); break
        else:
            xs.append(None)
    ok = all(x is not None for x in xs) and common.check_core(xs)
    log(level='published 479 core', all_blocks_invariant_under_order12_after_translation=ok, r=common.r, e=common.e)


def c2636(kdir):
    st = E.control_2636()
    log(control='2636', r=st.r, e=st.e)
    subs = st.invariant_subgroups()
    N41 = [N for N in subs if len(N) == 41][0]
    q1 = st.quotient(N41); ch1 = C.char_matrix(q1)
    Y, cnt = C.level_list(q1, 0, rows=None, chars=ch1)
    log(level='Z4^2 (G/Z41)', compressed_vectors=len(Y), kumar_reports=2954)
    cert = json.load(open(os.path.join(kdir, 'hadamard/2636/certificate.json')))
    xs = [S(s) for s in cert['sequences']]
    r = BB.solve_border(st, xs)
    log(level='border for Kumar core (our frames)', found=r is not None,
        hadamard=bool(r is not None and BB.is_hadamard(BB.assemble(st, xs, *r[:3]))))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--kumar', default='/home/user/vatskishan03/hadamard')
    ap.add_argument('--skip-final', action='store_true')
    a = ap.parse_args()
    c1916()
    c2636(a.kumar)
    c3964(a.kumar, final=not a.skip_final)
