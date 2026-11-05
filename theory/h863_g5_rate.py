r"""[M] Throughput of the G/Z5-level join (theory/h863_g5_solve.py) on the Codex y: w-solutions per second and per
pair visit (left + right list entries of both stages), uncapped stage beta.  One core.
usage: python3 theory/h863_g5_rate.py [secs] [na] [nb]"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h863_g5_solve import Problem, join22
secs = float(sys.argv[1]) if len(sys.argv) > 1 else 60; na = int(sys.argv[2]) if len(sys.argv) > 2 else 1024
nb = int(sys.argv[3]) if len(sys.argv) > 3 else 512
pb = Problem('data/h3452_stage1.json', 'q2'); Ta, Tb = pb.targets(); rng = np.random.default_rng(31)
t0 = time.time(); nw = 0; nal = 0; visits = 0
while time.time() - t0 < secs:
    H = [pb.sample_halves(i, na, rng) for i in range(4)]
    sig = [pb.alpha_sig(i, H[i][0], H[i][1]) for i in range(4)]
    hits = join22(*sig, Ta, cap=10 ** 7); visits += 2 * na * na
    for (a, b, c, d) in hits:
        nal += 1; pk = [a, b, c, d]; bs = []
        for i in range(4):
            Wa, Wb, pa = pb.splits(i, 0, H[i][0][pk[i]]); Wc, Wd, pc = pb.splits(i, 1, H[i][1][pk[i]])
            ka = rng.choice(len(Wa), size=nb, p=pa / pa.sum()); kc = rng.choice(len(Wc), size=nb, p=pc / pc.sum())
            bs.append(pb.beta_sig(i, Wa[ka] - Wb[ka], Wc[kc] - Wd[kc]))
        nw += len(join22(*bs, Tb, cap=10 ** 7)); visits += 2 * nb * nb
dt = time.time() - t0
print(json.dumps(dict(secs=round(dt, 1), na=na, nb=nb, alpha_solutions=nal, w_solutions=nw,
                      ms_per_w=round(1000 * dt / max(nw, 1), 3), log2_visits_per_w=round(np.log2(visits / max(nw, 1)), 2))))
