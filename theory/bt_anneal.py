r"""Simulated annealing over the orbit signs of a bordered-GS structure (heuristic production search; every core it
reports is checked exactly on G, bordered and assembled) [M for each core; the search itself proves nothing].

State: orbit signs eps_i (4 blocks, S-invariant).  Energy E = sum_{g != 0} (sum_i PAF_i(g) - F(g))^2, kept per
equation class c (weight |c|): PAF_i(d_c) = const + sum_{o < o'} B_c[o, o'] eps_o eps_o', B_c = A_c + A_c^T,
A_c[o, o'] = #{g in o : g + d_c in o'}.  Fields H_i[c, o] = sum_o' B_c[o, o'] eps_o' make a flip cost O(e) to score
and O(e r) to apply.  Metropolis at a geometric temperature ladder, restarts from random signs; E = 0 is a core.
Optional: --keep-q fixes nothing but starts each restart from a random Q-level solution's fibre (not used by default).
usage: python3 theory/bt_anneal.py (--control 3964|2636 | --census FILE --key KEY | --spec SPEC.json) OUTDIR
       [--hours h] [--seed s] [--steps N] [--t0 x] [--t1 y]
"""
import sys, os, json, time, math, argparse
import numpy as np
import numba as nb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_engine as E


def tables(st):
    """B (4, e, R, R) int32 with zero diagonal, diag constants (4, e), class weights w (e,), targets T (e,)"""
    G = st.G
    reps = st.cls_rep[1:]
    e = len(reps); R = max(st.r)
    B = np.zeros((4, e, R, R), np.int32); D = np.zeros((4, e), np.int64)
    for i in range(4):
        o = st.orb[i]
        for c, d in enumerate(reps):
            A = np.zeros((st.r[i], st.r[i]), np.int64)
            np.add.at(A, (o, o[G.add[:, d]]), 1)
            D[i, c] = int(np.trace(A))
            M = A + A.T
            np.fill_diagonal(M, 0)
            B[i, c, :st.r[i], :st.r[i]] = M
    w = st.cls_size[1:].astype(np.int64)
    T = st.cls_target[1:].astype(np.int64)
    return B, D, w, T


@nb.njit(cache=True)
def _paf(B, D, eps, r):
    nb4, e = D.shape
    P = np.zeros((4, e), np.int64)
    for i in range(4):
        for c in range(e):
            s = 0
            for a in range(r[i]):
                for b in range(a + 1, r[i]):
                    s += B[i, c, a, b] * eps[i, a] * eps[i, b]
            P[i, c] = D[i, c] + s
    return P


@nb.njit(cache=True)
def _fields(B, eps, r):
    e = B.shape[1]; R = B.shape[2]
    H = np.zeros((4, e, R), np.int64)
    for i in range(4):
        for c in range(e):
            for a in range(r[i]):
                s = 0
                for b in range(r[i]):
                    s += B[i, c, a, b] * eps[i, b]
                H[i, c, a] = s
    return H


@nb.njit(cache=True)
def _anneal(B, D, w, T, r, eps, steps, t0, t1, seed):
    """one annealing run from eps (modified in place); returns (best E, best eps, final E)"""
    np.random.seed(seed)
    e = B.shape[1]
    P = _paf(B, D, eps, r)
    S = np.zeros(e, np.int64)
    for c in range(e):
        s = -T[c]
        for i in range(4):
            s += P[i, c]
        S[c] = s
    En = 0
    for c in range(e):
        En += w[c] * S[c] * S[c]
    H = _fields(B, eps, r)
    best = En; beps = eps.copy()
    dP = np.zeros(e, np.int64)
    rtot = r[0] + r[1] + r[2] + r[3]
    lr = math.log(t1 / t0)
    for k in range(steps):
        temp = t0 * math.exp(lr * k / steps)
        u = np.random.randint(rtot)
        i = 0
        while u >= r[i]:
            u -= r[i]; i += 1
        o = u
        ev = eps[i, o]
        dE = 0
        for c in range(e):
            d = -2 * ev * H[i, c, o]
            dP[c] = d
            dE += w[c] * d * (d + 2 * S[c])
        if dE <= 0 or np.random.random() < math.exp(-dE / temp):
            eps[i, o] = -ev
            for c in range(e):
                S[c] += dP[c]
                for a in range(r[i]):
                    H[i, c, a] += -2 * ev * B[i, c, a, o]
            En += dE
            if En < best:
                best = En; beps[:, :] = eps
                if En == 0:
                    return best, beps, En
    return best, beps, En


def build(a):
    if a.control:
        return {'3964': E.control_3964, '2636': E.control_2636}[a.control](), dict(name=a.control)
    if a.census:
        import bt_census_run as R, bt_census2 as C2
        for f in a.census.split(','):
            for L in open(f):
                d = json.loads(L)
                if d.get('key') == a.key:
                    return (C2.build(d) if 'fam' in d else R.build(d)), d
        raise SystemExit('key not found')
    import bt_search as BS
    spec = json.load(open(a.spec))
    return BS.build_spec(spec), spec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('outdir')
    ap.add_argument('--control'); ap.add_argument('--census'); ap.add_argument('--key'); ap.add_argument('--spec')
    ap.add_argument('--hours', type=float, default=1.0); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--steps', type=int, default=2 * 10 ** 6)
    ap.add_argument('--t0', type=float, default=None); ap.add_argument('--t1', type=float, default=None)
    a = ap.parse_args()
    st, meta = build(a)
    os.makedirs(a.outdir, exist_ok=True)
    B, D, w, T = tables(st)
    r = np.array(st.r, np.int64)
    R = B.shape[2]
    # temperature scale: typical |dE| of a random flip from a random state
    rng = np.random.default_rng(a.seed)
    eps = np.zeros((4, R), np.int64)
    for i in range(4):
        eps[i, :r[i]] = rng.choice([-1, 1], r[i])
    P = _paf(B, D, eps, r)
    S = P.sum(0) - T
    H = _fields(B, eps, r)
    dEs = []
    for _ in range(200):
        i = rng.integers(4); o = rng.integers(r[i])
        dP = -2 * eps[i, o] * H[i, :, o]
        dEs.append(abs(int((w * dP * (dP + 2 * S)).sum())))
    scale = float(np.median(dEs))
    t0 = a.t0 if a.t0 else scale
    t1 = a.t1 if a.t1 else scale * 1e-3
    logf = os.path.join(a.outdir, f'anneal_s{a.seed}.jsonl')
    coref = os.path.join(a.outdir, 'cores.jsonl')
    t_end = time.time() + 3600 * a.hours
    run = 0; hist = {}
    with open(logf, 'a') as f:
        f.write(json.dumps(dict(event='start', name=str(meta.get('name', meta.get('key'))), r=list(map(int, r)),
                                e=int(len(T)), t0=t0, t1=t1, steps=a.steps, seed=a.seed, t=time.time())) + '\n')
    while time.time() < t_end:
        run += 1
        for i in range(4):
            eps[i, :] = 0
            eps[i, :r[i]] = rng.choice([-1, 1], r[i])
        t = time.time()
        best, beps, fin = _anneal(B, D, w, T, r, eps, a.steps, t0, t1, int(rng.integers(1 << 30)))
        hist[int(best)] = hist.get(int(best), 0) + 1
        rec = dict(event='run', run=run, best=int(best), final=int(fin), secs=round(time.time() - t, 2))
        if best == 0:
            xs = [st.expand(i, beps[i, :r[i]]) for i in range(4)]
            ok = st.check_core(xs)
            rec['core'] = bool(ok)
            with open(coref, 'a') as fc:
                fc.write(json.dumps(dict(eps=[beps[i, :r[i]].astype(int).tolist() for i in range(4)],
                                         rows=[int(x.sum()) for x in xs], check_core=bool(ok), run=run,
                                         seed=a.seed)) + '\n')
        with open(logf, 'a') as f:
            f.write(json.dumps(rec) + '\n')
        if best == 0:
            print(json.dumps(rec), flush=True)
    with open(logf, 'a') as f:
        f.write(json.dumps(dict(event='stop', runs=run, best_hist=sorted(hist.items())[:12], t=time.time())) + '\n')
    print(json.dumps(dict(runs=run, best_hist=sorted(hist.items())[:12])), flush=True)


if __name__ == '__main__':
    main()
