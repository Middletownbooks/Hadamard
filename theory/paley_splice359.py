r"""Splice a new v = 359 core into the published H(1436) border, and verify H H^T = 1436 I exactly.

The published matrix (data/hadamard_1436.csv.gz, hoa64-normalised) is sign-corrected exactly as in
scripts/bge3_template347.py: every core sub-block becomes rho(g) gam(h) s F_k(sigma(g, h)) with F_k the row-0 functions
of core row-block 0 and sigma in {h-g, g+h, g-h, -g-h}.  The border rows/columns are constant on the cosets of
H = <(0,0,2)> x Z11, so border-core orthogonality only sees the Q-compression (Q = G/H = Z2^3) of the core functions:
a new core with the SAME compression as F (coset by coset, block by block) gives a Hadamard matrix with the same
border and corner.  `selftest` rebuilds the corrected original from F.

usage: python3 theory/paley_splice359.py selftest
       python3 theory/paley_splice359.py ref           (print F in Paley coordinates: alpha/u/w on G' = Z2 x Z2 x Z8)
from python: T = Template359(); H = T.build(Fnew); T.check(H)"""
import sys, os, json, itertools as it
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from extract_group import load, G, solve

N = 1436; B = 7; MODS = (2, 2, 8, 11); OFF = 4 * B; V = N // 4; M_ = V - B

class Template359:
    def __init__(s):
        cache = os.path.join(ROOT, 'theory', 'logs', 'template359.npz')
        s.g = G(MODS); co = s.g.coords
        s.sig = {'h-g': s.g.idx(co[None, :, :] - co[:, None, :]), 'g+h': s.g.idx(co[:, None, :] + co[None, :, :]),
                 'g-h': s.g.idx(co[:, None, :] - co[None, :, :]), '-g-h': s.g.idx(-co[:, None, :] - co[None, :, :])}
        if os.path.exists(cache):
            d = np.load(cache, allow_pickle=False)
            s.Mt = d['Mt']; s.F = d['F']; s.layout = json.loads(str(d['layout'])); return
        M = load(N)
        def block_signs(M):
            E = []
            for c in range(4):
                sub = M[OFF:OFF + M_, OFF + c * M_:OFF + (c + 1) * M_]; f = solve(sub, s.g, c > 0)[0]; I = np.arange(M_)
                F = f[s.g.idx(co[I][:, None, :] + co[I][None, :, :])] if c > 0 else f[s.g.idx(co[I][None, :, :] - co[I][:, None, :])]
                T = sub * F; E.append(T[0, :] * T[0, 0])
            return np.concatenate(E)
        E = block_signs(M); D = block_signs(M.T.copy())
        Mt = M.copy(); Mt[:, OFF:] *= E[None, :]; Mt[OFF:, :] *= D[:, None]
        assert np.array_equal(Mt @ Mt.T, N * np.eye(N, dtype=Mt.dtype))
        core = Mt[OFF:, OFF:]; F = np.array([core[0, J * M_:(J + 1) * M_].copy() for J in range(4)])
        # real characters of G (products of the order-2 characters of the Z2, Z2, Z8 factors)
        chars = {}
        for e in it.product((0, 1), repeat=3):
            chars[''.join(map(str, e))] = np.where(((co[:, 0] * e[0] + co[:, 1] * e[1] + co[:, 2] * e[2]) % 2) == 0, 1, -1)
        layout = {}
        for Ib in range(4):
            for Jb in range(4):
                S = core[Ib * M_:(Ib + 1) * M_, Jb * M_:(Jb + 1) * M_]; found = None
                for k in range(4):
                    for sn, ix in s.sig.items():
                        R = S * F[k][ix]; rho = R[:, 0] * R[0, 0]; gam = R[0, :]
                        if np.array_equal(R, np.outer(rho, gam)):
                            cr = [n for n, ch in chars.items() if np.array_equal(rho, ch) or np.array_equal(rho, -ch)]
                            cg = [n for n, ch in chars.items() if np.array_equal(gam, ch) or np.array_equal(gam, -ch)]
                            if cr and cg:
                                found = dict(k=k, sigma=sn, rho=cr[0], rho_sign=int(rho[0] * chars[cr[0]][0]),
                                             gam=cg[0], gam_sign=int(gam[0] * chars[cg[0]][0])); break
                    if found: break
                assert found, (Ib, Jb)
                layout[f'{Ib}{Jb}'] = found
        s.Mt, s.F, s.layout = Mt, F, layout
        np.savez(cache, Mt=Mt, F=F, layout=json.dumps(layout))

    def chars(s):
        co = s.g.coords
        return {''.join(map(str, e)): np.where(((co[:, 0] * e[0] + co[:, 1] * e[1] + co[:, 2] * e[2]) % 2) == 0, 1, -1)
                for e in it.product((0, 1), repeat=3)}

    def build(s, Fnew):
        """Fnew: (4, 352) +-1 functions in the SAME coordinates as s.F (mixed radix 2,2,8,11, last fastest)"""
        H = s.Mt.copy(); ch = s.chars()
        for key, L in s.layout.items():
            Ib, Jb = int(key[0]), int(key[1])
            rho = L['rho_sign'] * ch[L['rho']]; gam = L['gam_sign'] * ch[L['gam']]
            H[OFF + Ib * M_:OFF + (Ib + 1) * M_, OFF + Jb * M_:OFF + (Jb + 1) * M_] = np.outer(rho, gam) * np.asarray(Fnew[L['k']])[s.sig[L['sigma']]]
        return H

    @staticmethod
    def check(H):
        return bool(np.array_equal(H @ H.T, N * np.eye(N, dtype=H.dtype)))

    def compression(s, Fk):
        """coset sums over H = <(0,0,2)> x Z11, indexed by (a, b, e mod 2)"""
        f = np.asarray(Fk).reshape(MODS); out = np.zeros((2, 2, 2), np.int64)
        for a, b, e in it.product(range(2), range(2), range(8)): out[a, b, e % 2] += f[a, b, e].sum()
        return out

def selftest():
    T = Template359()
    H = T.build(T.F)
    print(json.dumps({"rebuild_equals_corrected_original": bool(np.array_equal(H, T.Mt)), "hadamard": T.check(H),
                      "layout": T.layout, "compressions": [T.compression(T.F[k]).reshape(-1).tolist() for k in range(4)]}))

def ref():
    sys.path.insert(0, os.path.join(ROOT, 'theory'))
    from paley_algebra import components
    T = Template359()
    for k in range(4):
        a, u, w, ok = components(T.F[k].reshape(MODS), 11)
        print('ABCD'[k], 'in ansatz', ok, 'K', int((w != 0).sum()), 'rowsum', int(T.F[k].sum()))

if __name__ == "__main__":
    {"selftest": selftest, "ref": ref}[sys.argv[1]]()
