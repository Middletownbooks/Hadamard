r"""Stage-2 budget of the Z3^2-a slice (H' = Z3^2, -1 on one axis) for y's selected by stage-2 freedom [E].

H2 (stage-2 freedom given y) = number of (block, G'-orbit) with letter +-1 (two fibre patterns each; +-9, +-11 have
one).  Stage 1 of this slice is rich (report/b599_answer.md section 5), so y can be chosen; this draws K typical y
(theory/ab_typical_levels.typical_y_sampler), ranks them by H2 and budgets the top and bottom ones
(level_bits, sd_freedom): C2, C_S, h_D, D_def and the cost model 2^(D_def + C_S - H_S/2).
usage: python3 theory/b599_ysel.py [--K 40] [--out FILE]"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W
from ab_s1sd import cfg_from_json
from ab_typical_levels import typical_y_sampler, level_bits
from ab_levels import sd_freedom, subgroup_mask_axes

CFG = '{"v":599,"c":6,"b":5,"h":99,"m":594,"Q":[2,3],"H1":[3,3],"q":11,"s":2,"f":5,"mods":[2,3,3,3,11],"nQ":2,"U":[[2,2,2]]}'

if __name__ == "__main__":
    a = sys.argv[1:]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    K = int(opt('--K', '40')); outp = opt('--out', '')
    sl = W.build(cfg_from_json(CFG))
    cands = []
    for seed in range(K):
        Y = typical_y_sampler(sl, 100 + seed)()
        H2 = int(sum((np.abs(Y[k, bb]) == 1) for k in range(4) for bb in sl.breps))
        cands.append((H2, seed, Y))
    cands.sort(key=lambda t: -t[0])
    fo = open(outp, 'a') if outp else None
    print(json.dumps(dict(H2_list=[c[0] for c in cands])), flush=True)
    for H2, seed, Y in [cands[0], cands[1], cands[-1]]:
        r = level_bits(sl, Y, {'S(axis0)': subgroup_mask_axes(sl.Gp, [0])}, N=4000, seed=seed)
        sd = sd_freedom(sl, Y)
        CS = r['levels']['S(axis0)']['bits']; HS = r['H2'] - sd['h_D']; Dd = r['C2'] - CS - sd['h_D']
        rec = dict(seed=100 + seed, H2=r['H2'], C2=r['C2'], per_y=round(r['H2'] - r['C2'], 2), C_S=CS, h_D=sd['h_D'],
                   D_def=round(Dd, 2), cost=round(Dd + CS - HS / 2, 2), y=Y.tolist())
        print(json.dumps({k: v for k, v in rec.items() if k != 'y'}), flush=True)
        if fo: fo.write(json.dumps(rec) + "\n"); fo.flush()
