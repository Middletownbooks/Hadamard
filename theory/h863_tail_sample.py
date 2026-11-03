r"""[M] Samples for the tail calibration of the mixed-level Gaussian at v = 863: uniform configurations in the fibre
of a G/Z5 solution (y, w) (Codex certificate y; w from data/h863_g5_solutions.json), class PAF vectors saved as int16.
usage: python3 theory/h863_tail_sample.py W_INDEX N SEED OUT.npy"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h863_lib as L, h863_yw_estimate as YW
k, N, seed, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
Y, E, S = L.load_cert()
W = json.load(open('data/h863_g5_solutions.json'))['w']
ts, lf = YW.yw_fibre(S, W[k])
rng = np.random.default_rng(seed)
V = np.empty((N, 72), np.int16)
B = 16384; t0 = time.time()
for a in range(0, N, B):
    n = min(B, N - a)
    x = YW.sample_yw(E, ts, n, rng)
    V[a:a + n] = L.class_paf(x).sum(1)
    print(a + n, round(time.time() - t0, 1), flush=True)
np.save(out, V)
