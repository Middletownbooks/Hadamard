r"""Part C, gate C3 (cost): estimates for the index-4 Szolloesi plug-in on K x F_q, K = Z_nK (1291: nK = 19, q = 17)
[E; sizes exact].  Same model as the repository's deficit estimator (report/open_screen_answer.md section 2.3), written
for this family:
  stage 1: y-level count with the lag-0 energy (exact DP over cells, y-level, A and B symmetric), times the lattice
           Gaussian (theory/q523_gauss) of lags 1..nK/2 at the target, for y uniform on the energy shell;
  stage 2: for a synthetic y (uniform on the stage-1 energy shell), H2 = log2 #fibre (class values x_j(k) with sum sigma(k); 1, 4, 6 ways for |sigma| = 4, 2,
           0; A, B by symmetric pairs) and C2 = lattice-Gaussian bits of all G-level lag classes for uniform fibre
           points (the stage-1 combinations are then exact relations, which the lattice handles).
  E_tot = stage-1 log2E + H2 - C2; per-core cost of a staged search >= 2^max(0, C2 - H2) y's, each needing at least a
  2+2 join over fibres of 2^(H2/4) per block.
usage: python3 theory/h2396_c_estimate.py nK q [N] [ys]"""
import sys, os, math, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_gauss import gauss_bits
from ab_gauss import gauss_bits_at_mean
from h2396_c_s1cp import classes

def main(nK, q, N=4000, nys=2, seed=1):
    S, cls, cyc = classes(q); L = (q - 1) // 4; Nn = nK * q; T0 = 4 * Nn - 2
    rng = np.random.default_rng(seed)
    sig_vals = [-4, -2, 0, 2, 4]; mult = {4: 1, 2: 4, 0: 6}
    # cells: (block, list of K-positions, x0 allowed)
    cells = []
    for i in range(4):
        if i < 2:
            cells.append((i, [0], False))
            for k in range(1, nK // 2 + 1): cells.append((i, [k, nK - k], True))
        else:
            for k in range(nK): cells.append((i, [k], True))
    def letters(cell):
        _, pos, has0 = cell
        out = []
        for s in sig_vals:
            for a in ((1, -1) if has0 else (0,)): out.append(a + L * s)
        return out
    # exact y-level energy DP (counts distinct y-quadruples by total energy)
    back = [None] * (len(cells) + 1); back[-1] = np.zeros(T0 + 1); back[-1][0] = 1.0
    for c in range(len(cells) - 1, -1, -1):
        w = len(cells[c][1]); arr = np.zeros(T0 + 1)
        for y in letters(cells[c]):
            e = w * y * y
            if e <= T0: arr[e:] += back[c + 1][:T0 + 1 - e]
        back[c] = arr
    log2Ny = math.log2(back[0][T0])
    def sample_y():
        E = T0; y = [[0] * nK for _ in range(4)]
        for c, cell in enumerate(cells):
            w = len(cell[1]); opts = []; wts = []
            for v in letters(cell):
                e = w * v * v
                if e <= E and back[c + 1][E - e] > 0: opts.append(v); wts.append(back[c + 1][E - e])
            p = np.array(wts) / sum(wts); v = opts[rng.choice(len(opts), p=p)]
            for k in cell[1]: y[cell[0]][k] = v
            E -= w * v * v
        return y
    lags = list(range(1, nK // 2 + 1))
    V = np.zeros((N, len(lags)), np.int64)
    for t in range(N):
        y = np.array(sample_y())
        V[t] = [sum(int((y[i] * np.roll(y[i], -d)).sum()) for i in range(4)) for d in lags]
    g1 = gauss_bits(V, np.zeros(len(lags), np.int64))
    s1 = log2Ny - g1['bits'] if g1.get('in_lattice') else -math.inf
    out = dict(nK=nK, q=q, N=Nn, stage1=dict(log2Ny_energy=round(log2Ny, 2), bits=round(g1['bits'], 2), rank=g1['rank'],
                                             in_lattice=g1['in_lattice'], log2E=round(s1, 2)))
    # stage 2 at synthetic y's
    zrep = [0] + [next(z for z in range(1, q) if cls[z] == l) for l in range(4)]
    lagG = [(d, z) for d in range(0, nK // 2 + 1) for z in zrep if (d, z) != (0, 0)]
    res = []
    for _ in range(nys):
        y = sample_y()
        # decode y -> x0, sigma
        x0 = [[0] * nK for _ in range(4)]; sg = [[0] * nK for _ in range(4)]; H2 = 0.0
        for i in range(4):
            for k in range(nK):
                v = y[i][k]
                a = 0 if (i < 2 and k == 0) else (1 if (v - 1) % (2 * L) == 0 else -1)   # y = x0 + L sigma, sigma even
                s = (v - a) // L; assert a + L * s == v and s in (-4, -2, 0, 2, 4)
                x0[i][k] = a; sg[i][k] = s
        for i in range(4):
            for k in range(nK):
                if i < 2 and k > nK // 2: continue
                H2 += math.log2(mult[abs(sg[i][k])])
        # uniform fibre points -> G-level lag-class sums
        Vs = np.zeros((N, len(lagG)), np.int64)
        for t in range(N):
            tot = np.zeros(len(lagG), np.int64)
            for i in range(4):
                X = np.zeros((nK, q), np.int64)
                for k in range(nK):
                    kk = (nK - k) % nK if (i < 2 and k > nK // 2) else k
                    if kk != k: X[k] = X[kk]; continue
                    s = sg[i][k]; npos = (4 + s) // 2
                    cv = np.array([1] * npos + [-1] * (4 - npos)); rng.shuffle(cv)
                    X[k, 0] = x0[i][k]
                    for z in range(1, q): X[k, z] = cv[cls[z]]
                for j, (d, z) in enumerate(lagG):
                    tot[j] += int((X * np.roll(np.roll(X, -d, axis=0), -z, axis=1)).sum())
            Vs[t] = tot
        # a synthetic y does not solve stage 1, so its own stage-1 combinations miss the target: evaluate the
        # stage-2 constraint at the sample mean, as the repository estimator does for synthetic y (ab_gauss)
        g2 = gauss_bits_at_mean(Vs); b2 = g2.get('bits_at_mean')
        res.append(dict(H2=round(H2, 2), C2=b2, rank=g2.get('rank'), E_per_y=round(H2 - b2, 2) if b2 is not None else None))
    out['stage2_synthetic'] = res
    if res and res[0]['E_per_y'] is not None:
        Ey = float(np.mean([r['E_per_y'] for r in res])); H2m = float(np.mean([r['H2'] for r in res]))
        out['log2E_total'] = round(s1 + Ey, 2); out['log10E_total'] = round((s1 + Ey) * math.log10(2), 2)
        out['ys_per_core_log2'] = round(max(0.0, -Ey), 2)
        out['per_y_join_floor_log2'] = round(H2m / 4, 2)
        out['cost_per_core_floor_log2'] = round(max(0.0, -Ey) + H2m / 4, 2)
    print(json.dumps(out))

if __name__ == "__main__":
    a = sys.argv[1:]
    main(int(a[0]), int(a[1]), int(a[2]) if len(a) > 2 else 4000, int(a[3]) if len(a) > 3 else 2)
