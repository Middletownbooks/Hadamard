#!/usr/bin/env python3
"""Chunked, restart-safe version of theory/fable_tamic.py for t = 13 (or any t): the top two base-8 digits of the
assignment code select one of 64 chunks; finished chunks are appended to the log file and skipped on restart.
Usage: fable_tamic13.py t logfile"""
import sys, json, os, time
import numpy as np
from numba import njit

@njit(cache=True)
def chunk(t, top):
    ntm = 0; nout = 0
    blk = np.zeros(t, dtype=np.int8); sg = np.zeros(t, dtype=np.int8)
    T = np.zeros((4, t), dtype=np.int8)
    half = (t - 1) // 2
    sub = 8 ** (t - 3)
    first = -1
    for low in range(sub):
        c = top * sub + low
        blk[0] = 0; sg[0] = 1
        for i in range(1, t):
            d = c % 8; c //= 8
            blk[i] = d // 2; sg[i] = 1 - 2 * (d % 2)
        for j in range(4):
            for i in range(t): T[j, i] = 0
        for i in range(t): T[blk[i], i] = sg[i]
        ok = True
        for s in range(1, half + 1):
            acc = 0
            for i in range(t):
                ip = (i + s) % t
                if blk[i] == blk[ip]: acc += sg[i] * sg[ip]
            if acc != 0:
                ok = False; break
        if not ok: continue
        ntm += 1
        good = True
        for s in range(1, t):
            a12 = 0; a34 = 0; a13 = 0; a24 = 0; a14 = 0; a23 = 0
            for i in range(t):
                ip = (i + s) % t; im = (i - s) % t
                a12 += T[0, ip] * T[1, i] - T[0, im] * T[1, i]
                a34 += T[2, ip] * T[3, i] - T[2, im] * T[3, i]
                a13 += T[0, ip] * T[2, i] - T[0, im] * T[2, i]
                a24 += T[1, ip] * T[3, i] - T[1, im] * T[3, i]
                a14 += T[0, ip] * T[3, i] - T[0, im] * T[3, i]
                a23 += T[1, ip] * T[2, i] - T[1, im] * T[2, i]
            if a12 != a34 or a13 != -a24 or a14 != a23:
                good = False; break
        if good:
            nout += 1
            if first < 0: first = top * sub + low
    return ntm, nout, first

if __name__ == '__main__':
    t = int(sys.argv[1]); logf = sys.argv[2]
    done = set()
    if os.path.exists(logf):
        for l in open(logf):
            try: done.add(json.loads(l)['chunk'])
            except Exception: pass
    for top in range(64):
        if top in done: continue
        t0 = time.time()
        ntm, nout, first = chunk(t, top)
        with open(logf, 'a') as f:
            f.write(json.dumps(dict(t=t, chunk=top, tmatrices=int(ntm), with_relations=int(nout), first=int(first), secs=round(time.time() - t0))) + '\n')
