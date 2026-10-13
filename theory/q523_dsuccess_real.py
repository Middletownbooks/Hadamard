r"""D-success per S-solution measured on production S-solutions (report/paley_523_answer.md, sections 3.4 and 6) [E
from measured residuals, as theory/q523_dfine_calib.py].

Takes a random sample of real S-solutions from a production checkpoint (array (n, 4, F2) of pair-S indices written by
theory/q523_wagner.run), and for each one enumerates its D3-solutions and D-candidates with the reference D-check
(theory/q523_dcheck.py).  The full D-level signature sums of sampled candidates are fitted by the lattice-exact
Gaussian: log2 P(full D-level match | D3 match) = -fine_bits.  The D-success per S-solution is
log2(mean D-candidates per S-solution) - fine_bits, i.e. -D_def measured for the production y and lists.
usage: python3 theory/q523_dsuccess_real.py F q s y.json P4.npy [--n 400] [--seed 1]"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_dcheck import DChecker
from q523_dfine_calib import candidates_sums
from q523_gauss import gauss_bits
from q523_sdsolve import SDInst

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--n": "400", "--seed": "1"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(a[i]); i += 1
    F, q, s = map(int, pos[:3]); Y = json.load(open(pos[3])); Y = Y['y'] if isinstance(Y, dict) else Y
    inst = SDInst(F, q, s, Y); P4 = np.load(pos[4], mmap_mode='r')
    rng = np.random.default_rng(int(opt["--seed"])); n = min(int(opt["--n"]), len(P4))
    pick = np.sort(rng.choice(len(P4), size=n, replace=False))
    dc = DChecker(inst, int(opt["--seed"])); allv = []; nd3s = []; ncs = []; hds = []; t0 = time.time()
    for j in pick:
        p4 = [np.asarray(P4[j, k], np.int64) for k in range(4)]
        v, nd3, nc, hd = candidates_sums(dc, p4); allv += v; nd3s.append(nd3); ncs.append(nc); hds.append(hd)
    V = np.array(allv); g = gauss_bits(V, dc.TF)
    print(json.dumps(dict(stage="D-success on production S-solutions", file=os.path.basename(pos[4]), s_solutions=n,
                          mean_hD=round(float(np.mean(hds)), 2), d3_per_S_log2=round(math.log2(max(1e-9, np.mean(nd3s))), 2),
                          candidates_per_S_log2=round(math.log2(max(1e-9, np.mean(ncs))), 2), fine_bits=g['bits'],
                          mahal=g.get('mahal'), samples=len(V),
                          D_success_per_S_log2=round(math.log2(max(1e-9, np.mean(ncs))) - g['bits'], 2),
                          secs=round(time.time() - t0, 1))), flush=True)
