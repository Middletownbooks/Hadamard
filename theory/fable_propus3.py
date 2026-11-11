#!/usr/bin/env python3
"""Among the formal multiplications by Z_t of the propus-like type (A free, B and C symmetric, D = A), list those
whose OUTPUT is again of that type: two blocks symmetric (up to a shift) and the other two equal up to shift,
sign and reversal.  Such an assignment turns a propus quadruple of order v into one of order t*v.  [M]
Usage: fable_propus3.py t
"""
import sys, itertools
import numpy as np
sys.path.insert(0, 'theory')
import fable_tmult as F, fable_tmult_types as TT

def main(t):
    ty = ['f', 's', 's', '=0']
    F.symbols_for_type = lambda _x: TT.make(ty)
    sym, conj, letters = TT.make(ty)     # letters: (0,0)=A, (0,1)=A*, (1,0)=B, (2,0)=C
    # use the exhaustive code to get ALL solutions: re-implement the join here to collect every solution
    nsol, sols, cfgs, letters = F.run(t, 'x', verbose=False, want=10**9)
    star = {0: 1, 1: 0, 2: 2, 3: 3}      # letter index -> starred letter index
    def variants(c):
        """all images of a block under shift, sign, reversal (as tuples)"""
        out = set()
        for rev in (0, 1):
            d = c
            if rev:
                d = tuple((c[(-i) % t][0], star[c[(-i) % t][1]]) for i in range(t))
            for sh in range(t):
                e = tuple(d[(i + sh) % t] for i in range(t))
                for sg in (1, -1):
                    out.add(tuple((x * sg, l) for x, l in e))
        return out
    def is_sym(c):
        for sh in range(t):
            d = tuple(c[(i + sh) % t] for i in range(t))
            if all(d[(-i) % t] == (d[i][0], star[d[i][1]]) for i in range(t)):
                return True
        return False
    good = []
    for (p, q) in sols:
        idx = list(p) + list(q)
        blocks = [cfgs[i] for i in idx]
        syms = [is_sym(b) for b in blocks]
        ok = False
        for i, j in itertools.combinations(range(4), 2):
            rest = [k for k in range(4) if k not in (i, j)]
            if blocks[j] in variants(blocks[i]) and syms[rest[0]] and syms[rest[1]]:
                ok = True
        if ok: good.append(blocks)
    print(f"t={t}: {len(sols)} multiplications of the propus-like type; output again propus-like: {len(good)}")
    names = {0: 'A ', 1: 'A*', 2: 'B ', 3: 'C '}
    for blocks in good[:3]:
        for b in blocks:
            print('   ', ' '.join(('+' if e > 0 else '-') + names[l] for e, l in b))
        print('    --')

if __name__ == '__main__':
    main(int(sys.argv[1]))
