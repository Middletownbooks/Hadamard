r"""Calibration of the stage-2 Gaussian against exact counts [M vs E].

Inputs: the cfg; the complete list of cores of the slice (raw int8, 4 x r orbit signs per core, canonical up to block
order and sign; theory/ab_full_chars.py); the stage-1 solutions (raw int8, 4 x L1, canonical; theory/ab_s1_chars.py);
the per-y Gaussian records (JSON lines of theory/ab_s2_s4.py, keyed by the y's index).
For every stage-1 y the exact number of ORDERED lifts (x_1..x_4) with compression exactly y (the quantity the
Gaussian estimates) is the number of distinct orderings of each core's blocks whose compressions equal y in order,
summed over the cores; it is compared with 2^(H2 - C2) per y and in total.
usage: python3 theory/ab_s2_exactcmp.py CFG.json CORES.cbin YBIN GAUSS.jsonl"""
import sys, os, json, math, itertools
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_s1sd as SD


def main():
    a = sys.argv[1:]
    cfg = SD.cfg_from_json(open(a[0]).read()); sl = W.build(cfg)
    r = sl.r; L1 = len(sl.breps)
    X = np.fromfile(a[1], np.int8).reshape(-1, 4, r)
    Ys = np.fromfile(a[2], np.int8).reshape(-1, 4, L1)
    yidx = {Ys[i].tobytes(): i for i in range(len(Ys))}
    # compression of every block: y(b) = sum over the fibre of b; orbit values at the G'-orbit representatives
    Mb = np.zeros((sl.G.n, L1), np.int64)                                       # G-point -> its G'-point, if a rep
    for jj, bb in enumerate(sl.breps): Mb[np.nonzero(sl.bidx == bb)[0], jj] = 1
    T = np.asarray(sl.T); Yo = np.zeros((len(X), 4, L1), np.int64)
    for c0 in range(0, len(X), 8192):
        full = X[c0:c0 + 8192].astype(np.int64)[:, :, sl.orb]                  # (n, 4, |G|)
        assert (sl.G.paf(full).sum(1) == T[None, :]).all(), "a listed core fails sum_k PAF = T"
        Yo[c0:c0 + 8192] = full @ Mb                                            # y at the G'-orbit representatives
    exact = defaultdict(int); unmatched = 0
    for n in range(len(X)):
        blocks = [tuple(X[n, k].tolist()) for k in range(4)]
        ys = [Yo[n, k].astype(np.int8) for k in range(4)]
        sgn = [1 if int(Yo[n, k] @ sl.bsz) > 0 or (int(Yo[n, k] @ sl.bsz) == 0 and Yo[n, k, 0] > 0) else -1 for k in range(4)]
        ysg = [ys[k] * sgn[k] for k in range(4)]                                # block signs as in the stage-1 file
        xsg = [tuple((np.array(blocks[k]) * sgn[k]).tolist()) for k in range(4)]
        order = sorted(range(4), key=lambda k: ysg[k].tobytes())
        key = np.stack([ysg[k] for k in order]).tobytes()
        i = yidx.get(key)
        if i is None: unmatched += 1; continue
        seen = set()
        for p in itertools.permutations(range(4)):
            if all((ysg[p[t]] == ysg[order[t]]).all() for t in range(4)):
                seen.add(tuple(xsg[p[t]] for t in range(4)))
        exact[i] += len(seen)
    G = {}
    for l in open(a[3]):
        z = json.loads(l)
        if 'idx' in z and z.get('log2E') is not None: G[z['idx']] = z
    both = [i for i in G]
    d = [math.log2(exact[i]) - G[i]['log2E'] for i in both if exact.get(i, 0) > 0]
    zero = [i for i in both if exact.get(i, 0) == 0]
    tot_ex = sum(exact[i] for i in both); tot_g = sum(2.0 ** G[i]['log2E'] for i in both)
    out = dict(cores=int(len(X)), y_with_cores=len(exact), unmatched_cores=unmatched, y_compared=len(both),
               log2_exact_total=round(math.log2(tot_ex), 3) if tot_ex else None, log2_gauss_total=round(math.log2(tot_g), 3),
               exact_minus_gauss_total_bits=round(math.log2(tot_ex) - math.log2(tot_g), 3) if tot_ex else None,
               per_y_bits=dict(n=len(d), mean=round(float(np.mean(d)), 3), sd=round(float(np.std(d)), 3),
                               min=round(float(np.min(d)), 3), max=round(float(np.max(d)), 3)) if d else None,
               y_with_zero_exact=len(zero))
    print(json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
