r"""Templates (core functions F, layout, border K, Kc, corner C, compression) of the published b >= 3 matrices [M].
Generalises scripts/bge3_template347.py to 359 (non-split: Q = G / (2Z_8 x Z_11) = Z_2^3) and 479 (Q = Z_2^2).
Each template is checked by rebuilding the matrix with theory/ab_assemble.py and by the forced-corner formula.
Saves data/template{N}_ab.npz.
usage: python3 theory/ab_template.py"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'scripts'))
from extract_group import load, G as EG, solve
from ab_lib import AG
import ab_assemble as AS

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def block_signs(M, off, m, g):
    E = []
    for c in range(4):
        sub = M[off:off + m, off + c * m:off + (c + 1) * m]; f = solve(sub, g, c > 0)[0]; I = np.arange(m)
        F = f[g.idx(g.coords[I][:, None, :] + g.coords[I][None, :, :])] if c > 0 else f[g.idx(g.coords[I][None, :, :] - g.coords[I][:, None, :])]
        T = sub * F; E.append(T[0, :] * T[0, 0])
    return np.concatenate(E)

def extract(N, b, mods, qfun, Qmods):
    off = 4 * b; v = N // 4; m = v - b; g = EG(mods); Gm = AG(mods)
    M = load(N); E = block_signs(M, off, m, g); D = block_signs(M.T.copy(), off, m, g)
    Mt = M.copy(); Mt[:, off:] *= E[None, :]; Mt[off:, :] *= D[:, None]
    assert np.array_equal(Mt @ Mt.T, N * np.eye(N, dtype=Mt.dtype))
    core = Mt[off:, off:]; F = [core[0, J * m:(J + 1) * m].copy() for J in range(4)]
    sig = {n: AS.sigma_index(Gm, n) for n in ('h-g', 'g+h', 'g-h', '-g-h')}
    layout = {}
    for Ib in range(4):
        for Jb in range(4):
            S = core[Ib * m:(Ib + 1) * m, Jb * m:(Jb + 1) * m]; found = None
            for k in range(4):
                for sn, ix in sig.items():
                    R = S * F[k][ix]
                    rho = R[:, 0] * R[0, 0]; gam = R[0, :]
                    if np.array_equal(R, np.outer(rho, gam)):
                        if (rho == 1).all() and (np.abs(gam) == 1).all() and (gam == gam[0]).all():
                            found = dict(k=k, sigma=sn, rho='1', rho_sign=int(gam[0]), gam='1', gam_sign=1)
                        else:
                            found = dict(k=k, sigma=sn, rho='arr', rho_sign=1, gam='arr', gam_sign=1,
                                         rho_arr=rho.tolist(), gam_arr=gam.tolist())
                        break
                if found: break
            layout[f'{Ib}{Jb}'] = found
    assert all(layout.values()), ('layout not of the simple form', layout)
    qg = qfun(Gm.co); Qm = AG(Qmods); c = Qm.n
    K = np.zeros((off, 4 * c), np.int64); Kc = np.zeros((off, 4 * c), np.int64)
    for p in range(off):
        for J in range(4):
            for q in range(c):
                vals = set(Mt[p, off + J * m:off + (J + 1) * m][qg == q].tolist()); assert len(vals) == 1
                K[p, J * c + q] = vals.pop()
                vals = set(Mt[off + J * m:off + (J + 1) * m, p][qg == q].tolist()); assert len(vals) == 1
                Kc[p, J * c + q] = vals.pop()
    F = np.array(F, np.int64)
    comp = AS.compressions(F, mods, 0, qg=qg, Qmods=Qmods)
    C = Mt[:off, :off].astype(np.int64)
    H = AS.assemble(F, mods, 0, layout, K, Kc, C, qg=qg, Qmods=Qmods)
    ok = AS.is_hadamard(H) and np.array_equal(H.astype(np.int64), Mt)
    Cf = AS.forced_corner(comp, Qmods, layout, K, Kc)
    rec = dict(N=N, v=v, b=b, mods=mods, Qmods=Qmods, rebuilt_hadamard=bool(AS.is_hadamard(H)), equals_corrected=bool(np.array_equal(H.astype(np.int64), Mt)),
               forced_corner_equals_C=bool(Cf is not None and np.array_equal(Cf, C)), comp=comp.tolist(), rows=comp.sum(1).tolist())
    np.savez(os.path.join(ROOT, f'data/template{N}_ab.npz'), C=C, K=K, Kc=Kc, comp=comp, layout=json.dumps(layout), F=F,
             mods=np.array(mods), Qmods=np.array(Qmods), qg=qg)
    return rec

if __name__ == "__main__":
    cases = [(1388, 5, (2, 3, 57), lambda co: co[:, 0] * 3 + co[:, 1], (2, 3)),
             (1916, 3, (2, 2, 119), lambda co: co[:, 0] * 2 + co[:, 1], (2, 2)),
             (1436, 7, (2, 2, 8, 11), lambda co: co[:, 0] * 4 + co[:, 1] * 2 + (co[:, 2] % 2), (2, 2, 2))]
    for c_ in cases:
        print(json.dumps(extract(*c_)), flush=True)
