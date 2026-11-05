#!/usr/bin/env python3
"""Exact finite arithmetic screens for the obstructions in the identity report.

No random search, numerical solver, or floating-point decisions. The universal
conclusions come from the proofs; finite tests and necessary gates are labelled.
"""
import argparse
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import numpy as np

TARGETS = [149,179,187,215,257,259,271,275,305,329,331,341,355,359,367,371,385,389,391,
           395,415,429,445,457,461,467,469,509,515,521,523,551,553,571,581,595,599,605]


def prime_divisors(m):
    return [q for q in range(3,m+1,2) if m%q==0 and all(q%d for d in range(2,math.isqrt(q)+1))]


def fibre_gate(m,q):
    assert m%q==0 and q>=3 and q%2
    r=m//q
    rhs=16*r-12
    pairs=[]
    for A in range(min(16*r, rhs//(q+1))+1):
        remain=rhs-A*(q+1)
        if q==3:
            if remain==0:
                pairs.extend((A,B) for B in range(16*r-A+1))
        elif remain%(q-3)==0:
            B=remain//(q-3)
            if 0<=B<=16*r-A:
                pairs.append((A,B))
    return dict(m=m,q=q,r=r,rhs=rhs,solutions=len(pairs),first=None if not pairs else pairs[0],
                construction_status='fails' if not pairs else 'unknown (necessary gate passes)')


def three_stratum_gate(p,r):
    """All 512 sign tables, compressed only; both Q share the normal form."""
    m=p*r
    weights=[a*b for a in (1,(p-1)//2,(p-1)//2) for b in (1,(r-1)//2,(r-1)//2)]
    sums={sum(s*w for s,w in zip(signs,weights)) for signs in itertools.product((-1,1),repeat=9)}
    sums={z for z in sums if z*z<=4*m+12}
    ts=sorted(t for t in sums if t+2 in sums and t*(t+1)<=m)
    normal=[]
    for a in sorted(z for z in sums if z>0 and z*z<=m):
        for t in itertools.combinations_with_replacement(ts,3):
            if a*a+sum(v*(v+1) for v in t)==m:
                normal.append([a,*t])
    gates=[fibre_gate(m,p),fibre_gate(m,r)]
    reasons=[]
    if not normal:
        reasons.append('normal-form compression gate empty')
    if any(g['solutions']==0 for g in gates):
        reasons.append('quadratic/three-stratum fibre equation empty')
    if min(p,r)>=37:
        reasons.append('Theorem C forces the forbidden balanced tensor interior')
    # Keep the association between a cell's compression and its two fibre costs.
    # Separate attainable-sum and global fibre-count gates lose this information.
    options={}
    wp=np.array([1,(p-1)//2,(p-1)//2],np.int64)
    wr=np.array([1,(r-1)//2,(r-1)//2],np.int64)
    target=(16*r-12,16*p-12)
    for signs in itertools.product((-1,1),repeat=9):
        f=np.array(signs,np.int64).reshape(3,3)
        z=int(wp@f@wr)
        cp=(int(((wp@f)**2)@wr)-r)//(p-1)
        cr=(int(wp@((f@wr)**2))-p)//(r-1)
        if cp<=target[0] and cr<=target[1]:
            options.setdefault(z,set()).add((cp,cr))
    profiles=[]
    for a,*ts in normal:
        zs=[a]*4+[v for t in ts for v in [t+2,t,t,t]]
        dp={(0,0)}
        layers=[1]
        for z in zs:
            dp={(a+c,b+d) for a,b in dp for c,d in options.get(z,set())
                if a+c<=target[0] and b+d<=target[1]}
            layers.append(len(dp))
        profiles.append(dict(normal_form=[a,*ts],cell_compressions=zs,
                             cell_cost_options={str(z):sorted(options.get(z,set())) for z in sorted(set(zs))},
                             reachable_layers=layers,final_state_count=len(dp),
                             final_states_first_16=sorted(dp)[:16],passes=target in dp))
    if normal and not any(pf['passes'] for pf in profiles):
        reasons.append('compression-conditioned fibre-cost gate empty')
    return dict(p=p,r=r,m=m,cell_tables=512,attainable_compressions=sorted(sums),
                normal_forms=normal,fibre_gates=gates,conditioned_fibre_target=target,
                conditioned_profiles=profiles,status='fails' if reasons else 'unknown',reasons=reasons)


def boundary_enumeration(factors):
    """Exhaust every binary boundary filling, with interior B identically zero."""
    co=list(itertools.product(*(range(d) for d in factors)))
    boundary=[i for i,c in enumerate(co) if 0 in c]
    M=np.array([[1]],np.int64)
    for d in factors:
        M=np.kron(M,d*np.eye(d,dtype=np.int64)-np.ones((d,d),np.int64))
    M=M[np.ix_(boundary,boundary)]
    L=math.prod(d-1 for d in factors)
    b=len(boundary)
    minimum=None
    minimizers=0
    for start in range(0,1<<b,4096):
        codes=np.arange(start,min(start+4096,1<<b),dtype=np.int64)
        B=1-2*((codes[:,None]>>np.arange(b))&1)
        vals=np.sum((B@M)*B,axis=1)
        v=int(vals.min())
        if minimum is None or v<minimum:
            minimum=v
            minimizers=int((vals==v).sum())
        elif v==minimum:
            minimizers+=int((vals==v).sum())
    assert minimum>=L
    return dict(factors=factors,boundary_points=b,fillings=1<<b,
                minimum_primitive_boundary_energy=str(Fraction(minimum,L)),
                minimizing_fillings=minimizers)


def transversal_q_gate(m,qt):
    add=np.array([[(u^v) if qt=='v4' else (u+v)%4 for v in range(4)] for u in range(4)])
    rows=np.array([[m-2*((mask>>q)&1) for q in range(4)] for mask in range(16)],np.int64)
    pafs=np.array([[int(row@row[add[:,g]]) for g in range(4)] for row in rows])
    target=np.array([4*m+12,4*m,4*m,4*m])
    survivors=[]
    for ix in itertools.combinations_with_replacement(range(16),4):
        if np.array_equal(pafs[list(ix)].sum(0),target):
            survivors.append(ix)
    assert not survivors
    return dict(m=m,q=qt,row_occupancy_patterns=16,multisets=3876,passing_multisets=0)


def norm_product(a,b):
    # A shared scalar must be 1/12 from the nonprincipal quotient characters.
    principal=Fraction((16*a+12)*(16*b+12),12)
    desired=16*a*b+12
    gap=principal-desired
    assert gap==Fraction(16*a*b,3)+16*(a+b) and gap>0
    return dict(m1=a,m2=b,m=a*b,forced_scalar='1/12',
                principal_energy=str(principal),required=desired,excess=str(gap),status='fails')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',default='logs/openai_identity/obstructions.json')
    a=ap.parse_args()
    coverage=[]
    for m in TARGETS:
        gates=[fibre_gate(m,p) for p in prime_divisors(m)]
        extra=[]
        if m==275:
            extra=[dict(host='C5 x C5 x C11',gate=fibre_gate(m,25))]
        if m==605:
            extra=[dict(host='C11 x C11 x C5',gate=fibre_gate(m,121))]
        pairs=[three_stratum_gate(p,m//p) for p in range(3,math.isqrt(m)+1,2) if m%p==0]
        # For prime m the single-factor class is tested by the fibre gate.
        c_status=('fails' if all(p['status']=='fails' for p in pairs) else 'unknown') if pairs else gates[0]['construction_status']
        coverage.append(dict(order=16*m+12,n=4*m+3,m=m,tensor_class='fails',
                             scalar_norm_product_class='fails',transversal_class='fails',
                             quadratic_fibre_gates=gates,additional_elementary_abelian_hosts=extra,
                             three_stratum_status=c_status,all_arithmetic_two_factor_presentations=pairs))
    out=dict(script='theory/openai_identity_obstructions.py',
             boundary_exhaustions=[boundary_enumeration(fs) for fs in ([3],[3,3],[3,5],[5,5],[3,3,3])],
             scalar_product_first_ten=[norm_product(a,b) for a,b in itertools.combinations_with_replacement([1,3,5,7],2)],
             transversal_small=[transversal_q_gate(m,qt) for m in range(3,100,2) for qt in ('v4','z4')],
             single_field_small=[fibre_gate(m,m) for m in range(3,100,2)
                                 if len(prime_divisors(m))==1],coverage=coverage)
    out['three_stratum_small']=[three_stratum_gate(p,r) for p in range(3,34,2)
                               for r in range(p,100//p+1,2)]
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    Path(a.out).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(boundaries=out['boundary_exhaustions'],norm_pairs=10,transversal_cases=98,
                          targets=len(coverage),targets_killed_by_a_prime_fibre=sum(any(g['solutions']==0 for g in r['quadratic_fibre_gates']) for r in coverage),
                          three_stratum_targets_killed=sum(r['three_stratum_status']=='fails' for r in coverage))))


if __name__=='__main__':
    main()
