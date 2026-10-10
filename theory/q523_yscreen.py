r"""y-screen for the S/D route (report/paley_523_answer.md, Q3) [E, lattice-exact Gaussian; parity screen exact on the
sample lattice].

For each stage-1 solution y: stage-2 entropy H2, constraint bits C2 (all 36 coordinates), the S-level (C_S), the Z3
level (C_Z3), C_mid = C_S - C_Z3 (the Wagner second merge), C(D|S), h_D = E[log2 #D-options | uniform config] and
the D-deficit D_def = C(D|S) - h_D (the measured value at 523 agrees within 0.3 bit, report section 3.4), plus the
predicted cost of the route in pair operations:
   S-step  2^{(C_mid + D_def)/2 + 1}   (two-level Wagner tree, one middle value),
   D-check 2^{D_def} * t_D             (t_D = 2^{17}: measured D3-join size per S-solution at 523).
`in_lattice` false means the stage-2 target is outside the lattice spanned by sampled stage-2 differences: no
completion exists if that lattice is the full one (checked exactly by theory/q523_lattice.py).
usage: python3 theory/q523_yscreen.py F q s y1.json [y2.json ...] [--N 30000]"""
import sys, os, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_lib import Cyclo, pz, target_P
from q523_levels import sample_P, level_matrix, block_entropies
from q523_gauss import gauss_bits

def screen(F, q, s, Y, N=30000, seed=1):
    cy = Cyclo(q, s); v = F * q + 1; rng = np.random.default_rng(seed)
    P = sample_P(cy, Y, N, rng).reshape(N, -1); T = target_P(F, v, s).ravel()
    gal = [d for d in range(1, s + 1) if s % d == 0]
    C = lambda ds: gauss_bits(P @ level_matrix(F, s, ds, gal).T, level_matrix(F, s, ds, gal) @ T)
    divs = [d for d in range(1, F + 1) if F % d == 0]; F2 = F // 2
    ca = C(divs); cs = C([d for d in divs if F2 % d == 0]); c3 = C([1, 3])
    be = [block_entropies(cy, Y[k], F) for k in range(4)]
    H2 = sum(b["full"]["H"] for b in be); hD = sum(b["H(D|S)"] for b in be)
    out = dict(v=v, H2=round(H2, 2), C2=ca['bits'], yield_=round(H2 - ca['bits'], 2), in_lattice=ca['in_lattice'],
               C_S=cs['bits'], C_Z3=c3['bits'], C_mid=round(cs['bits'] - c3['bits'], 2), CDS=round(ca['bits'] - cs['bits'], 2),
               h_D=round(hD, 2), letters=sorted([abs(int(t)) for r in Y for t in r if abs(t) > 1]))
    out["D_def"] = round(out["CDS"] - hD, 2)
    out["cost_S_log2"] = round((out["C_mid"] + max(0, out["D_def"])) / 2 + 1, 2)
    out["cost_D_log2"] = round(max(0, out["D_def"]) + 17, 2)
    out["cost_log2"] = round(math.log2(2 ** out["cost_S_log2"] + 2 ** out["cost_D_log2"]), 2)
    return out

if __name__ == "__main__":
    a = sys.argv[1:]; N = 30000
    if '--N' in a: i = a.index('--N'); N = int(a[i + 1]); a = a[:i] + a[i + 2:]
    F, q, s = map(int, a[:3])
    for f in a[3:]:
        try:
            d = json.load(open(f))
        except Exception:
            continue
        Y = d['y'] if isinstance(d, dict) else d
        if Y is None: continue
        r = screen(F, q, s, Y, N); r['file'] = os.path.basename(f)
        print(json.dumps(r), flush=True)
