"""Independent second implementation of stage A (for W4): recursive generation, PSD from the PAF,
exact tuple-keyed dict join, no pair-PSD pruning. Compares its quadruple sets with wil_search.comp_quads.
Usage: wil_stageA_alt.py n p q"""
import sys, json, math, numpy as np
from wil_common import rowsum_cases

def gen(m, k):
    hm = (m + 1) // 2; out = []
    def rec(pre):
        if len(pre) == hm: out.append(tuple(pre)); return
        for v in range(-k, k + 1, 2): rec(pre + [v])
    rec([]); return out

def full(h, m):
    return [h[min(i, m - i)] for i in range(m)]

def paf(x):
    m = len(x); return [sum(x[i] * x[(i + s) % m] for i in range(m)) for s in range(m)]

def psd_ok(P, n):
    m = len(P)
    for f in range(1, m):
        val = sum(P[s] * math.cos(2 * math.pi * s * f / m) for s in range(m))   # PSD = DFT of PAF (real, exact up to fp)
        if val > 4 * n + 1e-6: return False
    return True

def quads(n, m, case):
    k = n // m; hm = (m + 1) // 2; lists = {}
    for h in gen(m, k):
        x = full(h, m); s = sum(x)
        if s not in case: continue
        P = paf(x)
        if psd_ok(P, n): lists.setdefault(s, []).append((h, tuple(P[1:hm])))
    A, B, C, D = (lists.get(s, []) for s in case)
    tab = {}
    for c, pc in C:
        for d, pd in D: tab.setdefault(tuple(-(u + v) for u, v in zip(pc, pd)), []).append((c, d))
    res = set()
    for a, pa in A:
        for b, pb in B:
            for c, d in tab.get(tuple(u + v for u, v in zip(pa, pb)), []): res.add((a, b, c, d))
    return res

if __name__ == "__main__":
    import wil_search as ws
    n, p, q = map(int, sys.argv[1:4]); rng = np.random.default_rng(0); out = []
    for case in rowsum_cases(n):
        for m in (p, q):
            alt = quads(n, m, case)
            _, Q = ws.comp_quads(n, m, case, rng)
            main = {tuple(tuple(int(t) for t in blk) for blk in qd) for qd in Q}
            out.append({"case": case, "m": m, "alt": len(alt), "main": len(main), "equal": alt == main})
            print(json.dumps(out[-1]), flush=True)
    print(json.dumps({"n": n, "all_equal": all(o["equal"] for o in out)}))
