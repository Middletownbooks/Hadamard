r"""Design "SP" inside the Paley-coset ansatz (p = 3 mod 4): one swap pair + two co-symmetric blocks.

Per coset g' of Z_p:  x_k(g', .) = alpha*delta + u*chi + w*eps.  With y = alpha + (p-1) w, beta = alpha - w the reduced
system is (report/paley_design_answer.md section 2, verified in theory/paley_algebra.py):
    S1:  sum_k PAF(y_k)                    = T_J := T0 + (p-1) Ts          (on G')
    S2:  sum_k PAF(beta_k) + p PAF(u_k)    = T0 - Ts  (= 4v delta)
    X :  sum_k CC(beta_k, u_k) - CC(u_k, beta_k) = 0
Swap pair (blocks A, B) from two +-1 sequences a, c on G' (all cosets twisted):
    A: alpha = a, u = c;    B: alpha(g) = -c(s - g), u(g) = a(s - g)
  => X_A + X_B = 0,  S1 gets P := PAF(a) + PAF(c),  S2 gets (p+1) P.
Co-symmetric block (center t): alpha, u, w all satisfy f(t - g) = f(g)  =>  X_k = 0.
Eliminating P:  E_C + E_D = (T0 - Ts) - (p+1) T_J  with  E_k = PAF(beta_k) + p PAF(u_k) - (p+1) PAF(y_k);
then P = T_J - PAF(y_C) - PAF(y_D) must be a "Golay-type" sum PAF(a) + PAF(c).

usage: python3 theory/paley_swap.py search v b p mods... [--Hgen a.b.c] [--cap N]   (exhaustive for small n')"""
import sys, os, json, itertools as it, time
from collections import defaultdict
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from paley_ansatz import LABELS, LC
import paley_ls as LSm

def paf_g(P, f):
    return np.array([(f * f[P.add[:, h]]).sum() for h in range(P.n)])

def paf_batch(P, F):
    """F: (N, n) -> (N, n)"""
    return np.stack([(F * F[:, P.add[:, h]]).sum(1) for h in range(P.n)], 1)

def cosym_orbits(P, t):
    seen, orbs = set(), []
    for g in range(P.n):
        if g in seen: continue
        g2 = int(P.sub[t, g])                         # t - g
        o = sorted({g, g2}); seen.update(o); orbs.append(o)
    return orbs

def cosym_configs(P, t, maxK=None):
    """all co-symmetric label vectors (one label per orbit) -> arrays alpha, u, w of shape (N, n)"""
    orbs = cosym_orbits(P, t); k = len(orbs)
    idx = np.array(list(it.product(range(8), repeat=k)), np.int64)          # (8^k, k)
    if maxK is not None:
        const = np.array([LABELS[L][1] == 0 for L in range(8)])
        nconst = sum(const[idx[:, j]] * len(o) for j, o in enumerate(orbs))
        idx = idx[nconst <= maxK]
    lab = np.zeros((len(idx), P.n), np.int64)
    for j, o in enumerate(orbs): lab[:, o] = idx[:, [j]]
    A = LC[lab]
    return lab, A[..., 0], A[..., 1], A[..., 2]

def search(v, b, p, mods, Hgen=None, cap=None, max_report=5):
    P = LSm.Prob(v, b, p, mods, Hgen); n = P.n; t0 = time.time()
    T0, Ts = P.T[:, 0].astype(np.int64), P.T[:, 1].astype(np.int64)
    TJ = T0 + (p - 1) * Ts; TS2 = T0 - Ts
    TE = TS2 - (p + 1) * TJ
    # norm identity: (p+1)K - 4 N_- = 4 n' - 4 b  (all constants live in C, D)
    out = {"v": v, "b": b, "p": p, "mods": list(mods), "n'": n, "H'": getattr(P, "H", None)}
    Kmax = (4 * n - 4 * b) // max(1, p - 3) if p > 3 else 2 * n
    # +-1 pair table: PAF(a) for all a in {+-1}^n (n <= 20)
    assert n <= 20
    allpm = 1 - 2 * ((np.arange(2 ** n)[:, None] >> np.arange(n)) & 1)
    Pa = paf_batch(P, allpm)
    hashA = defaultdict(list)
    for i, r in enumerate(Pa): hashA[r.tobytes()].append(i)
    sols = []; nCD = 0
    centers = list(range(n))
    tried = set()
    for tC in centers:
        for tD in centers:
            # translation of a whole block by s maps center t -> t - 2s; skip equivalent center classes
            keyc = (int(P.G.idx(P.G.co[tC] % 2)) if False else tC, tD)
            labC, aC, uC, wC = cosym_configs(P, tC, Kmax)
            labD, aD, uD, wD = cosym_configs(P, tD, Kmax)
            def Evec(a, u, w):
                y = a + (p - 1) * w; be = a - w
                return paf_batch(P, be) + p * paf_batch(P, u) - (p + 1) * paf_batch(P, y), paf_batch(P, y), (w != 0).sum(1), ((w != 0) & (w == -a)).sum(1)
            EC, PyC, KC, NmC = Evec(aC, uC, wC); ED, PyD, KD, NmD = Evec(aD, uD, wD)
            H = defaultdict(list)
            for j, r in enumerate(ED): H[r.tobytes()].append(j)
            for i in range(len(EC)):
                need = (TE - EC[i]).tobytes()
                for j in H.get(need, []):
                    K = KC[i] + KD[j]; Nm = NmC[i] + NmD[j]
                    if (p + 1) * K - 4 * Nm != 4 * n - 4 * b: continue
                    nCD += 1
                    Pv = TJ - PyC[i] - PyD[j]
                    # Golay-type pair a, c with PAF(a) + PAF(c) = Pv
                    for ia in range(len(allpm)):
                        hit = hashA.get((Pv - Pa[ia]).tobytes())
                        if not hit: continue
                        for ic in hit[:1]:
                            a = allpm[ia]; c = allpm[ic]
                            for s in range(n):
                                labA = np.array([LABELS.index((int(a[g]), 1, int(c[g]))) for g in range(n)])
                                labB = np.array([LABELS.index((int(-c[P.sub[s, g]]), 1, int(a[P.sub[s, g]]))) for g in range(n)])
                                lab = np.stack([labA, labB, labC[i], labD[j]])
                                ok, X = LSm.verify(P, lab)
                                if ok:
                                    sols.append(dict(lab=lab.tolist(), tC=tC, tD=tD, s=s, K=int(K), Nminus=int(Nm)))
                                    break
                            if len(sols) >= max_report: break
                        if len(sols) >= max_report: break
                    if len(sols) >= max_report: break
                if len(sols) >= max_report: break
            if len(sols) >= max_report: break
        if len(sols) >= max_report: break
    out.update(CD_solutions_scanned=nCD, solutions_found=len(sols), secs=round(time.time() - t0, 1))
    return out, sols

if __name__ == "__main__":
    a = sys.argv[1:]
    cmd = a[0]; opt = {"--Hgen": None, "--cap": None}; pos = []; i = 1
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    v, b, p = pos[:3]; mods = tuple(pos[3:])
    Hgen = None if opt["--Hgen"] is None else [tuple(map(int, e.split('.'))) for e in opt["--Hgen"].split(',')]
    out, sols = search(v, b, p, mods, Hgen)
    print(json.dumps(out), flush=True)
    if sols:
        fn = os.path.join(ROOT, 'theory', 'logs', f"paley_swap_v{v}_b{b}_p{p}_{'x'.join(map(str, mods))}.json")
        json.dump(dict(out, solutions=sols), open(fn, 'w'))
