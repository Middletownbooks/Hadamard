r"""971 x469: exhaustive decision of the dual compression (along Z_10 onto Z_97) [P for the reduction, M].

Reduction (theory/nd_971_dual.py, [P]): per block X(0) = R_k = (2, 0, 0, 0) and X = 2Y on the eight classes C_j of
<81> in Z_97^*, Y in [-5, 5]^8, sum_j Y_j = 0; the stage-1A equations are
    lag 0:  sum_{k,j} Y_{k,j}^2 = 80;
    lag h (one per class, 8 classes):  sum_k [ X_k(0) Y_{k,cl(h)} + sum_{j,j'} N_h(j, j') Y_{k,j} Y_{k,j'} ] = -10,
with N_h(j, j') = #{z in C_j : z + h in C_j'} (the PAF divided by 4).
Symmetries of the equations [P]: the class shift j -> pi(j) (multiplication by a generator of Z_97^*; it permutes the
lags and the target is constant), negation of each row-0 block, and permutations of the row-0 blocks.  Every solution
is therefore equivalent to one with e_1 <= e_2 <= e_3 (energies e_k = sum_j Y_{k,j}^2), Y_1, Y_2 negation-canonical
and Y_3 canonical under <pi, -1> (16 images).  Search: theory/nd_ejoin.py with left = (block 0: all Y, block 1:
negation-canonical, e_1 <= ER / 2) and right = (block 2: negation-canonical, block 3: <pi, -1>-canonical, e_2 <= e_3).
Every hit is re-verified by the full cyclic PAF on Z_97.  An empty result decides the dual compression, hence 971 x469,
as empty (all of its stage-1 solutions are one orbit, and every core compresses to a dual solution).
Restart-safe: one JSONL line per energy class EL, run by --workers processes (EL classes dealt round-robin).
usage: python3 theory/nd_971_ejoin.py --dir DIR [--workers 4]  |  --plant N [--E 40]"""
import sys, os, json, math, time, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nd_971_dual import CL, LAGS, expand, q
from nd_ejoin import ejoin4

def class_perm():
    """pi: class index of g * (representative of class j), g a generator of Z_97^*"""
    g = next(a for a in range(2, q) if len({pow(a, k, q) for k in range(q - 1)}) == q - 1)
    rep = [next(z for z in range(1, q) if CL[z] == j) for j in range(8)]
    return np.array([CL[g * rep[j] % q] for j in range(8)]), g

PI, GEN = class_perm()

def nmats():
    N = np.zeros((len(LAGS), 8, 8), np.int64)
    for i, h in enumerate(LAGS):
        for z in range(1, q):
            w = (z + h) % q
            if w != 0: N[i, CL[z], CL[w]] += 1
    return N

NM = nmats()

def all_Y(E=80):
    """every Y in [-5, 5]^8 with sum 0 and sum of squares <= E (vectorised)"""
    v = np.arange(-5, 6, dtype=np.int8)
    G = np.stack(np.meshgrid(*([v] * 7), indexing='ij'), -1).reshape(-1, 7)
    y8 = -G.sum(1, dtype=np.int64)
    ok = np.abs(y8) <= 5
    G = G[ok]; y8 = y8[ok]
    e = (G.astype(np.int64) ** 2).sum(1) + y8 ** 2
    ok = e <= E
    Y = np.concatenate([G[ok].astype(np.int64), y8[ok, None]], 1)
    return Y

def signature(Y, x0):
    """columns: energy, then (PAF / 4) at the 8 lag classes"""
    Y = np.asarray(Y, np.int64)
    out = np.empty((len(Y), 1 + len(LAGS)), np.int64); out[:, 0] = (Y * Y).sum(1)
    for i, h in enumerate(LAGS):
        out[:, 1 + i] = x0 * Y[:, CL[h]] + np.einsum('nj,jk,nk->n', Y, NM[i], Y)
    return out

def code(Y):
    return ((Y + 5) * (11 ** np.arange(8))).sum(1)

def neg_canonical(Y):
    return code(Y) <= code(-Y)

def shiftneg_canonical(Y):
    c = code(Y); best = c.copy(); Z = Y
    for t in range(8):
        for sg in (1, -1):
            best = np.minimum(best, code(sg * Z))
        Z = Z[:, PI]
    return c == best

def check_signature():
    """the signature formula against the direct cyclic PAF of the expanded sequence (nd_971_dual.expand)"""
    rng = np.random.default_rng(3)
    for x0 in (2, 0):
        Y = rng.integers(-5, 6, (50, 8))
        X = expand(x0, Y)
        direct = np.stack([(X * np.roll(X, -h, axis=1)).sum(1) for h in LAGS], 1)
        assert (4 * signature(Y, x0)[:, 1:] == direct).all()
        assert ((X * X).sum(1) == x0 * x0 + 48 * (Y * Y).sum(1)).all()
    # pi is the class map of multiplication by GEN, and the lag classes permute accordingly
    for z in range(1, q): assert CL[GEN * z % q] == PI[CL[z]]
    return True

def lists(E=80):
    Y = all_Y(E)
    S0 = signature(Y, 2); S1 = signature(Y, 0)
    nc = neg_canonical(Y); sc = shiftneg_canonical(Y)
    return Y, [S0, S1[nc], S1[nc], S1[sc]], [np.arange(len(Y)), np.nonzero(nc)[0], np.nonzero(nc)[0], np.nonzero(sc)[0]]

def verify(Ys):
    X = [expand(x0, y)[0] for x0, y in zip((2, 0, 0, 0), Ys)]
    P = sum(np.array([(x * np.roll(x, -h)).sum() for h in range(q)]) for x in X)
    return bool(P[0] == 3844 and (P[1:] == -40).all())

_SH = {}
def _prepare_shared(E, T):
    from nd_ejoin import prepare, hash_weights
    Y, S, idx = lists(E)
    w = hash_weights(len(T))
    P1 = prepare(S[1], E, w)
    _SH.update(Y=Y, idx=idx, prepared=[prepare(S[0], E, w), P1, P1, prepare(S[3], E, w)])

def _worker(args):
    els, path, E, Tlist = args
    Y, idx = _SH['Y'], _SH['idx']
    T = np.array(Tlist, np.int64)
    def cb(EL, found):
        rec = dict(EL=EL, solutions=[[Y[idx[k][t[k]]].tolist() for k in range(4)] for t in found])
        with open(path, "a") as fo: fo.write(json.dumps(rec) + "\n")
    sols, st = ejoin4(None, T, order=(0, 1, 2, 3), bmax_frac=0.5, cled=True, els=els, done_cb=cb, prepared=_SH['prepared'])
    return st

def run(d, workers=4, E=80, log=print):
    os.makedirs(d, exist_ok=True); path = os.path.join(d, "dual971_els.jsonl")
    assert check_signature()
    T = [E] + [-10] * len(LAGS)
    done = set()
    if os.path.exists(path):
        for line in open(path):
            try: done.add(json.loads(line)['EL'])
            except Exception: pass
    todo = [EL for EL in range(0, E + 1, 2) if EL not in done]
    # deal the energy classes round-robin by estimated size (largest first)
    log(json.dumps(dict(stage="start", todo=todo, done=sorted(done))))
    chunks = [todo[i::workers] for i in range(workers)]
    t0 = time.time()
    _prepare_shared(E, T)
    log(json.dumps(dict(stage="lists", sizes=[len(p_['idx']) for p_ in _SH['prepared']], secs=round(time.time() - t0, 1))))
    if workers > 1:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(workers) as pool:
            for st in pool.imap_unordered(_worker, [(c, path, E, T) for c in chunks if c]):
                log(json.dumps(dict(stage="worker_done", **st, wall=round(time.time() - t0, 1))))
    else:
        log(json.dumps(dict(stage="worker_done", **_worker((todo, path, E, T)))))
    sols = []; els = set()
    for line in open(path):
        r = json.loads(line); els.add(r['EL']); sols += r['solutions']
    ver = [verify(s_) for s_ in sols]
    log(json.dumps(dict(stage="done", energy_classes=len(els), solutions=len(sols), verified=sum(ver), secs=round(time.time() - t0, 1))))
    return sols

def plant(seed, E=40):
    """random canonical-form quadruple with total energy E; target = its own signature; must be found"""
    rng = np.random.default_rng(seed)
    Y, S, idx = lists(E)
    e = (Y * Y).sum(1)
    while True:
        k = [int(rng.integers(len(Y))) for _ in range(4)]
        if sum(int(e[i]) for i in k) == E: break
    Ys = [Y[i].copy() for i in k]
    # canonicalise: sort blocks 1..3 by energy, class shift + sign making block 3 canonical, negation of blocks 1, 2
    rest = sorted(Ys[1:], key=lambda y: int((y * y).sum()))
    Ys = [Ys[0]] + rest
    best = None
    for t in range(8):
        Z = [y.copy() for y in Ys]
        for _ in range(t): Z = [z[PI] for z in Z]
        for sg in (1, -1):
            c = int(code((sg * Z[3])[None])[0])
            if best is None or c < best[0]: best = (c, [Z[0], Z[1], Z[2], sg * Z[3]])
    Ys = best[1]
    for i in (1, 2):
        if not neg_canonical(Ys[i][None])[0]: Ys[i] = -Ys[i]
    T = signature(Ys[0][None], 2)[0] + sum(signature(Ys[i][None], 0)[0] for i in (1, 2, 3))
    sols, st = ejoin4(S, T, order=(0, 1, 2, 3), bmax_frac=0.5, cled=True)
    found = any(all((Y[idx[k][t[k]]] == Ys[k]).all() for k in range(4)) for t in sols)
    return dict(seed=seed, E=E, energies=[int((y * y).sum()) for y in Ys], solutions=len(sols), planted_found=found, work_log2=st['work_log2'], secs=st['secs'])

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--dir": "", "--workers": "4", "--plant": "0", "--E": "40"}; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: i += 1
    if int(opt["--plant"]):
        assert check_signature()
        for sd in range(int(opt["--plant"])): print(json.dumps(plant(sd, int(opt["--E"]))), flush=True)
        sys.exit()
    run(opt["--dir"], int(opt["--workers"]), log=lambda s: print(s, flush=True))
