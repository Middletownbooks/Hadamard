#!/usr/bin/env python3
"""Independent audit of the final compression-conditioned 512-table gate.

Imports no builder or screening module. Expands all tables, computes subgroup
PAF sums by literal integer shifts, and uses an 8+8 meet-in-the-middle join instead
of the screening script's 16-step dynamic program. Fixed short audit cases.
"""
import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path
import numpy as np


def audit(p,r):
    m=p*r
    # Representatives 0, positive half, negative half; no character identities.
    lp=np.array([0]+[1]*((p-1)//2)+[2]*((p-1)//2))
    lr=np.array([0]+[1]*((r-1)//2)+[2]*((r-1)//2))
    bits=np.arange(512,dtype=np.int64)[:,None]>>np.arange(9)
    tables=(2*(bits&1)-1).reshape(512,3,3)
    cells=tables[:,lp[:,None],lr[None,:]]
    sums=cells.sum(axis=(1,2))
    SP=np.zeros(512,dtype=np.int64)
    SR=np.zeros(512,dtype=np.int64)
    for lag in range(p):
        SP+=(cells*np.roll(cells,lag,axis=1)).sum(axis=(1,2))
    for lag in range(r):
        SR+=(cells*np.roll(cells,lag,axis=2)).sum(axis=(1,2))
    assert np.array_equal(SP,(cells.sum(axis=1)**2).sum(axis=1))
    assert np.array_equal(SR,(cells.sum(axis=2)**2).sum(axis=1))
    assert np.all((SP-r)%(p-1)==0) and np.all((SR-p)%(r-1)==0)
    cp=(SP-r)//(p-1)
    cr=(SR-p)//(r-1)
    options={}
    target=(16*r-12,16*p-12)
    for z,c,d in zip(sums,cp,cr):
        if c<=target[0] and d<=target[1]:
            options.setdefault(int(z),set()).add((int(c),int(d)))
    Z=set(int(z) for z in sums)
    ts=sorted(t for t in Z if t+2 in Z and t*(t+1)<=m)
    normal=[]
    for a in sorted(z for z in Z if 0<z and z*z<=m):
        for ts3 in itertools.combinations_with_replacement(ts,3):
            if a*a+sum(t*(t+1) for t in ts3)==m:
                normal.append((a,*ts3))

    @lru_cache(None)
    def half(zs):
        if len(zs)==1:
            return frozenset(options.get(zs[0],set()))
        cut=len(zs)//2
        left,right=half(zs[:cut]),half(zs[cut:])
        return frozenset((a+c,b+d) for a,b in left for c,d in right
                         if a+c<=target[0] and b+d<=target[1])

    profiles=[]
    for a,*ts3 in normal:
        zs=[a]*4+[z for t in ts3 for z in (t+2,t,t,t)]
        left=half(tuple(zs[:8]))
        right=half(tuple(zs[8:]))
        matches=[(a,b) for a,b in left if (target[0]-a,target[1]-b) in right]
        profiles.append(dict(normal_form=[a,*ts3],left_states=len(left),right_states=len(right),
                             matching_sums=len(matches)))
    return dict(p=p,r=r,m=m,all_512_tables_integer_subgroup_paf_checked=True,
                normal_forms=profiles,passing_profiles=sum(p['matching_sums']>0 for p in profiles))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',default='logs/openai_identity/conditioned_gate_audit.json')
    a=ap.parse_args()
    # Small calibrations followed by the three final targets and both real-core sizes.
    cases=[(3,3),(3,5),(3,7),(5,5),(3,11),(5,7),(3,13),(7,7),(3,17),(5,11),
           (11,17),(17,23),(17,35),(7,17),(13,19)]
    results=[audit(p,r) for p,r in cases]
    for rec in results[10:]:
        assert rec['passing_profiles']==0,rec
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    Path(a.out).write_text(json.dumps(dict(script='theory/openai_identity_gatecheck.py',results=results),indent=2)+'\n')
    print(json.dumps(dict(cases=len(results),tables_checked=512*len(results),
                          final_targets_and_known_sizes_excluded=True)))


if __name__=='__main__':
    main()
