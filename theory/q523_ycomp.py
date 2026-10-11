r"""Composition screen without stage-1 solutions (report/paley_523_answer.md, section 5) [E].

The D-deficit of a y is D_def = C(D|S) - h_D.  h_D is exact from the pair types (y(g), y(g + 9)) (section 5); C(D|S) =
C2 - C_S is a Gaussian quantity whose covariance part depends only on the coset alphabets and positions.  For a
composition without a known stage-1 solution we place its letters at random (random signs; stage 1 not imposed) and
evaluate the Gaussian bits at the sample mean (Mahalanobis term dropped), "CDS0"; on real stage-1 solutions the same
quantity is compared with the true C(D|S) to calibrate the offset.
usage: python3 theory/q523_ycomp.py F q s [--real y1.json ...] [--comps "13:6,15:4,1:62;..."] [--fakes 3] [--N 20000]"""
import sys, os, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_lib import Cyclo, target_P
from q523_levels import sample_P, level_matrix, block_entropies
from q523_gauss import gauss_bits

def bits0(g):
    """Gaussian bits with the Mahalanobis term removed"""
    return g['bits'] - 0.5 * g['mahal'] ** 2 / math.log(2)

def cds(F, q, s, Y, N, seed=1):
    cy = Cyclo(q, s); v = F * q + 1; rng = np.random.default_rng(seed)
    P = sample_P(cy, Y, N, rng).reshape(N, -1); T = target_P(F, v, s).ravel()
    gal = [d for d in range(1, s + 1) if s % d == 0]; divs = [d for d in range(1, F + 1) if F % d == 0]; F2 = F // 2
    def C(ds, at_sample=False):
        Mx = level_matrix(F, s, ds, gal); V = P @ Mx.T
        return gauss_bits(V, V[0] if at_sample else Mx @ T)          # a sample point is always in the lattice
    ca = C(divs); cs = C([d for d in divs if F2 % d == 0])
    ca0 = C(divs, True); cs0 = C([d for d in divs if F2 % d == 0], True)
    hD = sum(block_entropies(cy, Y[k], F)["H(D|S)"] for k in range(4))
    out = dict(h_D=round(hD, 2), CDS0=round(bits0(ca0) - bits0(cs0), 2))
    if ca.get('in_lattice') and math.isfinite(ca['bits']) and math.isfinite(cs['bits']):
        out['CDS'] = round(ca['bits'] - cs['bits'], 2)
    return out

def fake_y(F, comp, rng, pair_big=0):
    """random placement of the composition (random signs), optionally with pair_big same-sign big-big (g, g+9) pairs"""
    mags = []
    for a, n in comp.items(): mags += [a] * n
    mags = np.array(mags); rng.shuffle(mags); Y = (mags * rng.choice([-1, 1], len(mags))).reshape(4, F)
    return Y.tolist()

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--comps": "", "--fakes": "3", "--N": "20000"}; real = []; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        elif a[i] == '--real':
            i += 1
            while i < len(a) and not a[i].startswith('--'): real.append(a[i]); i += 1
        else: pos.append(int(a[i])); i += 1
    F, q, s = pos; N = int(opt["--N"])
    for f in real:
        Y = json.load(open(f)); Y = Y['y'] if isinstance(Y, dict) else Y
        r = cds(F, q, s, Y, N); r.update(kind="real", file=os.path.basename(f)); print(json.dumps(r), flush=True)
    rng = np.random.default_rng(7)
    for cs_ in [c for c in opt["--comps"].split(';') if c]:
        comp = {int(t.split(':')[0]): int(t.split(':')[1]) for t in cs_.split(',')}
        for t in range(int(opt["--fakes"])):
            r = cds(F, q, s, fake_y(F, comp, rng), N); r.update(kind="fake", comp=cs_); print(json.dumps(r), flush=True)
