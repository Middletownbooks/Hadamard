"""Independent check of theory/ab_s2_hzq.screen_y on a small slice: for real stage-1 y's, enumerate EVERY lift of
every block (product of fibre patterns), form xbar^psi on H' x Z_q directly on group points, take its PAF by FFT,
and count 4-block combinations meeting the H' x Z_q-level target exactly; compare with screen_y (equal counts).
usage: python3 hzq_bruteforce_check.py CFG.json YBIN N"""
import sys, os, json, math, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../../theory'))
import numpy as np, ab_s1sd as SD, ab_gauss as AGs, ab_s2_hzq as HZ
cfg = SD.cfg_from_json(open(sys.argv[1]).read()); S = HZ.setup(cfg); sl = S['sl']; L1 = S['L1']; Hg = S['Hg']
Y = np.fromfile(sys.argv[2], np.int8).reshape(-1, 4, L1); N = int(sys.argv[3])
qg = sl.G.co[:, 0] * 2 + sl.G.co[:, 1]; gidx = Hg.idx(sl.G.co[:, 2:]); reps = S['reps']
for i in np.random.default_rng(3).choice(len(Y), N, replace=False):
    y = Y[i].astype(np.int64); res = []
    for p in range(4):
        sig = []
        for k in range(4):
            pats = [AGs.fibre_patterns(sl, j, int(y[k, j]))[0] for j in range(L1)]
            idx = np.array(list(itertools.product(*[range(len(P)) for P in pats])), np.int64)
            lab = np.zeros((len(idx), sl.r), np.int8)
            for j, P in enumerate(pats):
                for t, (o, _) in enumerate(sl.parts[j]): lab[:, o] = P[idx[:, j], t]
            X = lab[:, sl.orb].astype(np.float64) * np.array(HZ.CHARS[p])[qg][None, :]   # psi-weighted, on G
            Xb = np.zeros((len(X), Hg.n)); np.add.at(Xb.T, gidx, X.T)                     # sum over Q: xbar on H''
            F = np.fft.fftn(Xb.reshape((len(X),) + Hg.mods), axes=tuple(range(1, 1 + len(Hg.mods))))
            R = np.rint(np.fft.ifftn(np.abs(F) ** 2, axes=tuple(range(1, 1 + len(Hg.mods)))).real).reshape(len(X), Hg.n)
            P_ = R[:, reps].astype(np.int64)[:, :-1]
            u, c = np.unique(P_, axis=0, return_counts=True); sig.append((u, c.astype(float)))
        T = S['Tpsi'][p][:-1]; d = len(T); dt = np.dtype([('f%d' % t, np.int64) for t in range(d)])
        sv = lambda R: np.ascontiguousarray(R, np.int64).view(dt).ravel()
        def conv(a, b):
            R = (a[0][:, None, :] + b[0][None, :, :]).reshape(-1, d); v = (a[1][:, None] * b[1][None, :]).ravel()
            u, inv = np.unique(sv(R), return_inverse=True); c = np.zeros(len(u)); np.add.at(c, inv.ravel(), v); return u, c
        A = conv(sig[0], sig[1]); B = conv(sig[2], sig[3])
        q = sv(T[None, :] - B[0].view(np.int64).reshape(-1, d)); pos = np.minimum(np.searchsorted(A[0], q), len(A[0]) - 1); hit = A[0][pos] == q
        brute = float((A[1][pos[hit]] * B[1][hit]).sum()); fast, by = HZ.screen_y(S, y, p)
        res.append(dict(psi=p, brute=brute, screen=fast, equal=brute == fast))
    print(json.dumps(dict(i=int(i), results=res)), flush=True)
