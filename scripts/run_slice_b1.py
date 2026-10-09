"""Worker: alternate designed-gadget stage 1 (theory/paley_slicepipe.guided, stage 2 off, records y's) with the exact
stage-2 engine (scripts/paley_s2join.solve_y) on every new y.  Stops at the first verified Hadamard core.
usage: run_slice_b1.py v p n mu K Nminus seed0 [rounds] [round_secs] [kmax] [dcomp] [thresh_x10]"""
import sys, os, json, time
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'theory')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import paley_slicepipe as SP, paley_s2join as J
a = list(map(int, sys.argv[1:]))
v, p, n, mu, K, Nm, seed0 = a[:7]
rounds = a[7] if len(a) > 7 else 1000; rsecs = a[8] if len(a) > 8 else 600; kmax = a[9] if len(a) > 9 else K
dcomp = a[10] if len(a) > 10 else 0
thresh = (a[11] if len(a) > 11 else 0) / 10.0
t0 = time.time(); tot_y = 0
for r in range(rounds):
    seed = seed0 + 1000 * r
    SP.guided(v, p, n, mu, K, Nm, max_runs=10 ** 6, seed=seed, kmax=kmax, stage2=False, time_cap=rsecs, thresh=thresh)
    f = os.path.join(ROOT, 'theory', 'logs', f'paley_slicepipe_stage1_v{v}_mu{mu}_K{K}_s{seed}.npz')
    if not os.path.exists(f):
        print(json.dumps(dict(round=r, seed=seed, ys=0, total_y=tot_y, secs=round(time.time() - t0))), flush=True); continue
    d = np.load(f); As, Ws = d['A'], d['W']; tot_y += len(As); nsol = 0
    for i in range(len(As)):
        rec, sols = J.solve_y(As[i], Ws[i], v, p, n, mu, dcomp=dcomp, maxhits=4)
        if sols:
            print(json.dumps(dict(FOUND=True, v=v, seed=seed, y=i, npz=f, **{k: rec[k] for k in ('solutions', 'core_ok', 'hadamard') if k in rec})), flush=True)
            if rec.get('hadamard'): sys.exit(0)
    print(json.dumps(dict(round=r, seed=seed, ys=len(As), total_y=tot_y, secs=round(time.time() - t0))), flush=True)
