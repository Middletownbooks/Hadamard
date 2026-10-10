r"""Independent brute-force check of theory/p571_p19gadgets.py on small distributions: enumerate ALL signed gadget
sets per block (one gadget at position 0 per nonempty block, every sign pattern, no excess pruning, no join), apply the
gadget conditions directly (|t_k| <= 1, s(h) in {0, -1}, s(15) = 0, #{h in 1..14: s = -1} = (K - #odd)/2, Lemma 6 at
the trivial character, X-feasibility of every labelling), canonicalise, and compare with the enumerator's classes.
usage: python3 theory/p571_p19gadget_check.py classes.pkl"""
import sys, os, pickle, itertools as it, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p571_p19gadgets as G
n = 30

def block_options(K, block0):
    if K == 0: return [np.zeros(n, np.int64)]
    out = []
    for rest in it.combinations(range(1, n), K - 1):
        pos = (0,) + rest
        for sg in it.product((1, -1), repeat=K):
            if abs(sum(sg)) > 1: continue
            w = np.zeros(n, np.int64); w[list(pos)] = sg; out.append(w)
    return out

def brute(dist, Nm, K):
    opts = [block_options(Kk, k == 0) for k, Kk in enumerate(dist)]
    nspec = (K - sum(x % 2 for x in dist)) // 2
    found = set()
    paf = lambda w: np.array([int((w * np.roll(w, -h)).sum()) for h in range(16)])
    P = [[paf(w) for w in o] for o in opts]
    for sel in it.product(*[range(len(o)) for o in opts]):
        s = sum(P[k][sel[k]] for k in range(4))
        if not all(s[h] in (0, -1) for h in range(1, 15)) or s[15] != 0: continue
        if sum(1 for h in range(1, 15) if s[h] == -1) != nspec: continue
        Ws = [opts[k][sel[k]] for k in range(4)]; ts = [int(w.sum()) for w in Ws]
        if not G.lemma6_trivial(ts): continue
        gpos = [(k, g) for k in range(4) for g in np.nonzero(Ws[k])[0]]
        for jm in it.combinations(range(len(gpos)), Nm):
            eps = [np.zeros(n, np.int64) for _ in range(4)]
            for i_, (k, g) in enumerate(gpos): eps[k][g] = -1 if i_ in jm else 1
            if G.x_feasible(Ws, eps, s, K): found.add(G.canon_config(Ws, eps))
    return found

if __name__ == "__main__":
    classes = pickle.load(open(sys.argv[1], 'rb'))
    for K, Nm, dist in ((6, 1, (0, 2, 2, 2)), (6, 1, (1, 2, 2, 1)), (6, 1, (2, 2, 1, 1)), (7, 6, (1, 2, 2, 2)),
                        (7, 6, (2, 2, 2, 1)), (6, 1, (0, 3, 2, 1)), (7, 6, (0, 3, 2, 2))):
        t0 = time.time(); b = brute(dist, Nm, K)
        e = {c for c, info in classes.items() if info['K'] == K and tuple(info['dist']) == dist}
        print(json.dumps(dict(K=K, dist=dist, brute=len(b), enumerator=len(e), equal=(b == e), secs=round(time.time() - t0, 1))), flush=True)
