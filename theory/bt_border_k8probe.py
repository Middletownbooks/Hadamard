r"""Probe of the k = 8 border problem at 751 [M for each found border; a failure is a failure of the frame family]:
for the first Q-level classes of the best 751 k8 structure, 300 balanced Kronecker frames H_4 (x) H_8[1:] (columns
permuted per block) with the admissible set by theory/bt_border.admissible_balanced and a clique of 28 orthogonal
vectors.  Output: logs/bt/border_k8_751.log.
usage: python3 theory/bt_border_k8probe.py
"""
import sys, time, numpy as np, json
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_border as BB, bt_gates as BG, bt_ladder as LD, bt_census2 as C2
recs = [json.loads(L) for f in ('logs/bt/census2/k8_w0.jsonl', 'logs/bt/census2/k8_w1.jsonl') for L in open(f)]
d = sorted([d for d in recs if d['n'] == 751 and d['verdict'].startswith('passes')], key=lambda d: -(d.get('log2E_Q') or -1e9))[0]
st = C2.build(d)
g = BG.Gates(st); qs = g.q_level(patterns=LD.patterns_for(st))
out = []
for rows, sol, L in qs:
    for s0 in sol[:2]:
        z = np.array([g.qQ.full_from_reduced(i, L[i][0][s0[i]]) for i in range(4)])
        Q = BB.QData(st); Mc = BB.gs_array([zz for zz in z], Q.add, Q.neg)
        sizes = []
        t = time.time(); found = None
        for k, Rh in enumerate(BB.frames_balanced(8, count=300, seed=7)):
            S = BB.admissible_balanced(Mc, Rh, 7, 8)
            sizes.append(len(S))
            if len(S) >= 28:
                Xc = BB.clique(S, 28, tmax=30)
                if Xc is not None:
                    found = k; break
        rec = dict(key=d['key'], rows=list(map(int, rows)), z=z.tolist(), frames=len(sizes), max_admissible=int(max(sizes)),
                   frames_with_28=int(sum(x >= 28 for x in sizes)), border_frame=found, secs=round(time.time() - t, 1))
        print(json.dumps(rec), flush=True)
        out.append(rec)
