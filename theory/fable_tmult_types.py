#!/usr/bin/env python3
"""Scan all letter-relation types of GS quadruples for a formal multiplication by Z_t (t = 3 exhaustive).

A type assigns to each of A, B, C, D one of: 'f' (free), 's' (symmetric: X* = X), or a reference
'=K' / '*K' (equal to an earlier block K or to its reverse).  For each type we run the exhaustive
meet-in-the-middle of theory/fable_tmult.py at t = 3.  [M]
"""
import sys, itertools
sys.path.insert(0, 'theory')
import fable_tmult as F

def make(typ):
    # typ: list of 4 entries
    sym = {}; conj = {}; nxt = 0; letters = []
    for k, e in enumerate(typ):
        if e == 'f':
            sym[(k, 0)] = nxt; sym[(k, 1)] = nxt + 1; conj[nxt] = nxt + 1; conj[nxt + 1] = nxt
            letters += [(k, 0), (k, 1)]; nxt += 2
        elif e == 's':
            sym[(k, 0)] = nxt; sym[(k, 1)] = nxt; conj[nxt] = nxt
            letters += [(k, 0)]; nxt += 1
        else:
            op, K = e[0], int(e[1])
            if op == '=':
                sym[(k, 0)] = sym[(K, 0)]; sym[(k, 1)] = sym[(K, 1)]
            else:
                sym[(k, 0)] = sym[(K, 1)]; sym[(k, 1)] = sym[(K, 0)]
    return sym, conj, letters

def types():
    out = []
    opts = ['f', 's']
    for a in opts:
        for b in opts + ['=0', '*0']:
            for c in opts + ['=0', '*0', '=1', '*1']:
                for d in opts + ['=0', '*0', '=1', '*1', '=2', '*2']:
                    t = [a, b, c, d]
                    # drop references to a referencing block (keep canonical: refs point to f/s blocks only)
                    ok = all(e in ('f', 's') or t[int(e[1])] in ('f', 's') for e in t)
                    # '*K' with K symmetric is the same as '=K'
                    ok = ok and not any(e[0] == '*' and t[int(e[1])] == 's' for e in t if e not in ('f', 's'))
                    if ok:
                        out.append(t)
    return out

if __name__ == '__main__':
    t = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    orig = F.symbols_for_type
    res = []
    for ty in types():
        F.symbols_for_type = lambda _x, ty=ty: make(ty)
        nsol, sols, cfgs, letters = F.run(t, 'x', verbose=False)
        # real dimension of the type
        dim = sum(2 if e == 'f' else 1 if e == 's' else 0 for e in ty)
        print(' '.join(f"{e:>2}" for e in ty), f" dim={dim}  solutions={nsol}", flush=True)
        res.append((ty, nsol))
