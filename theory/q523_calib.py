r"""Calibration of the lattice-exact Gaussian on the S-step and D-check at full size [M vs E].

For a stage-1 y and a chosen Z3 parent, build the S-step lists (theory/q523_sdsolve.py), predict the number of
S-solutions as  prod |L_k| * 2^{-bits},  bits = lattice-Gaussian of the sum of signatures of random 4-tuples drawn from
the lists, and count the S-solutions exactly with the 2+2 join.  Optionally runs the D-check on every S-solution and
compares the number of D3-solutions with its prediction.
usage: python3 theory/q523_calib.py F q s y.json [--K 32] [--beta 1] [--top 24] [--seed 0] [--parent i] [--dcheck n]"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_sdsolve import SDInst, choose_parent, fibre_list, seq_sig, join_all, d_check
from q523_gauss import gauss_bits

def predict(sigs, target, M=40000, seed=0):
    rng = np.random.default_rng(seed)
    V = sum(s[rng.integers(len(s), size=M)] for s in sigs)
    g = gauss_bits(V, target)
    return g, sum(math.log2(len(s)) for s in sigs) - g['bits']

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--K": "32", "--beta": "1", "--top": "24", "--seed": "0", "--parent": "0", "--dcheck": "0"}
    pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(a[i]); i += 1
    F, q, s = map(int, pos[:3]); Y = json.load(open(pos[3])); Y = Y['y'] if isinstance(Y, dict) else Y
    K, beta, top, seed, pi = int(opt["--K"]), float(opt["--beta"]), int(opt["--top"]), int(opt["--seed"]), int(opt["--parent"])
    inst = SDInst(F, q, s, Y); rng = np.random.default_rng(seed)
    log = lambda d: print(json.dumps(d), flush=True)
    groups = [[inst.row_groups(k, r) for r in range(3)] for k in range(4)]
    parents = choose_parent(inst, groups, beta, top, seed, want=pi + 1, log=log)
    par, sc = parents[pi]
    lists = [fibre_list(inst, groups, k, par[k], K, beta, rng) for k in range(4)]
    sigs = [seq_sig(L[0], inst.cy.T, inst.nl) for L in lists]
    g, pred = predict(sigs, inst.tsigS(), seed=seed)
    t0 = time.time(); hits, st = join_all(sigs, inst.tsigS(), seed=seed, tbits=24)
    log(dict(stage="S-step calibration", F=F, q=q, K=K, beta=beta, list_log2=[round(math.log2(len(x)), 2) for x in sigs],
             gauss_bits=g['bits'], mahal=g.get('mahal'), predicted_log2=round(pred, 2), found=int(len(hits)),
             found_log2=round(math.log2(len(hits)), 2) if len(hits) else None, secs=round(time.time() - t0, 1), probes=st['probes']))
    nd = int(opt["--dcheck"])
    if nd:
        tot3 = 0; hD = []; t1 = time.time(); nsol = 0
        for hq in hits[:nd]:
            pidx4 = [lists[k][1][hq[k]] for k in range(4)]
            sols, info = d_check(inst, pidx4, seed=seed)
            tot3 += info['d3_solutions']; hD.append(sum(info['d_options_log2'])); nsol += len(sols)
        log(dict(stage="D-check calibration", checked=min(nd, len(hits)), d3_solutions=tot3, full_solutions=nsol,
                 mean_hD=round(float(np.mean(hD)), 2), secs=round(time.time() - t1, 1)))
