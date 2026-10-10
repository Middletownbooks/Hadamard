"""Planted test of stage 2 (theory/p571_hsolve.stage2): take a real b = 1 core with a Z6 factor (v = 43 on Z42), compute
its y and the S-class of each block, run stage 2 on random S-class sublists forced to contain the true classes, and
require the original core to be recovered (S-level join, D-level check, assembly, exact core check).
usage: python3 theory/p571_test_stage2.py"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p571_hsolve import Inst, stage2, Stage2Block
from p571_ternary import load_core

W = np.exp(2j * np.pi / 3)
def key(z): return (int(round(z.real * 2)), int(round(z.imag * 2 / np.sqrt(3))))

if __name__ == "__main__":
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    x, m = load_core(os.path.join(root, 'data', 'core_v43_b1_paley.txt'))
    inst = Inst(m, []); S = inst.S; reps = [o[0] for o in S.orbs]
    Xc = S.to_cube(x); y = Xc.sum(2)[:, :, reps]
    Zp = ((Xc * W ** np.arange(3)[:, None]).sum(2) / -2)[:, :, reps]          # Z' (4, 2, R)
    force = []
    for k in range(4):
        b = Stage2Block(inst, y[k]); digit = []
        for i in range(inst.R):
            s = Zp[k, 0, i] + Zp[k, 1, i]
            sa = int(round((s.real + s.imag / np.sqrt(3))))              # s = a + b w: b = 2 Im/sqrt3, a = Re + b/2
            sb = int(round(2 * s.imag / np.sqrt(3)))
            digit.append(b.Sval[i].index((sa, sb)))
        idx = 0
        for i in range(inst.R): idx = idx * int(b.nS[i]) + digit[i]
        force.append(idx)
    rec, cores = stage2(inst, y, sub=300, seed=1, tbits=16, force=force)
    found = any(np.array_equal(c, x) for c in cores)
    print(json.dumps(dict(**rec, recovered_original=bool(found))))
