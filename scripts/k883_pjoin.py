"""Partial joins on the production flat pools (packed uint64 orbit signs): count 4-tuples (A,B,C,D) from the pools of
row sums (a1,a2,a3,a4) whose class sums vanish on the first k classes of a random order, k = 1..kmax.
usage: pjoin.py v d a1,a2,a3,a4 cache_dir P kmax seed"""
import sys, json, numpy as np, resource
resource.setrlimit(resource.RLIMIT_AS, (10 << 30, 10 << 30))
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from redlib import Red
v, d = int(sys.argv[1]), int(sys.argv[2]); prof = [int(x) for x in sys.argv[3].split(",")]
cache, P, kmax, seed = sys.argv[4], int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
R = Red(v, d); r = R.r; rng = np.random.default_rng(seed)
paf = {}
for a in set(prof):
    u = np.load(f"{cache}/pool_{a}.npy", mmap_mode="r")
    u = np.asarray(u[np.sort(rng.choice(len(u), P, replace=False))])
    Y = np.where(((u[:, None] >> np.arange(r, dtype=np.uint64)[None, :]) & np.uint64(1)) == 1, 1, -1).astype(np.int64)
    paf[a] = R.paf(R.expand(Y))[:, R.lags].astype(np.int64)
order = rng.permutation(R.L)
for k in range(1, kmax + 1):
    idx = order[:k]; W = rng.integers(1, 2 ** 62, size=k).astype(np.uint64)
    h = [(paf[a][:, idx].astype(np.uint64) * W).sum(1) for a in prof]
    AB = (h[0][:, None] + h[1][None, :]).ravel()
    CD = np.sort((np.uint64(0) - (h[2][:, None] + h[3][None, :])).ravel())
    lo = np.searchsorted(CD, AB, "left"); hi = np.searchsorted(CD, AB, "right"); cnt = int((hi - lo).sum())
    del AB, CD
    tot = 4 * np.log2(P)
    print(json.dumps(dict(k=k, matches=cnt, tuples_log2=round(tot, 1), bits=round(float(tot - np.log2(cnt)), 2) if cnt else None)), flush=True)
