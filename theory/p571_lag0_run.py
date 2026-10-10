"""Run the exact lag-0 (J, rho) filter (theory/p571_z10profiles.Lag0) on the top-weight Z10 profiles and save survivors.
usage: python3 theory/p571_lag0_run.py m gen [--top 400] [--out file.pkl]"""
import sys, os, json, math, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p571_hsolve import Inst
from p571_z10profiles import profiles, Lag0

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--top": "400", "--out": ""}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    inst = Inst(pos[0], pos[1:]); t0 = time.time()
    sols, _ = profiles(inst); sols.sort(key=lambda t: -t[1])
    Jtot = (4 * inst.m - 8 - 8 * inst.q) // 8; L = Lag0(inst); keep = []
    for n, (P, w) in enumerate(sols[:int(opt["--top"])]):
        if L.ok(P, Jtot): keep.append((P, w))
        if (n + 1) % 50 == 0:
            print(json.dumps(dict(checked=n + 1, survivors=len(keep), secs=round(time.time() - t0, 1))), flush=True)
    print(json.dumps(dict(summary=True, profiles=len(sols), checked=min(len(sols), int(opt["--top"])), Jtot=Jtot,
                          survivors=len(keep), log2w_top=(round(math.log2(keep[0][1]), 1) if keep else None))), flush=True)
    if opt["--out"]: pickle.dump(keep, open(opt["--out"], 'wb'))
