"""Independent all-lag CP-SAT check of a fixed S14 Z43 compression fibre.

No imports from the structure/search code. Q = Z2^2, H = Q x Z5.
S14 acts as (q1,q2,t) -> (q2,q1,-t). Each compressed value is
epsilon + 7 sum of six signs. The target is 2936 at zero, -516 at
nonzero elements with q=0, and 172 otherwise. Scope: the input z only.
"""
import argparse, json, time
import numpy as np
from ortools.sat.python import cp_model


def check(z, secs, max_abs=None):
    values=sorted({a+7*b for a in (-1,1) for b in range(-6,7,2)})
    if max_abs is not None:values=[v for v in values if abs(v)<=max_abs]
    model=cp_model.CpModel()
    # Independent enumeration of quotient involution orbits.
    reps=[];labels=[]
    for q in range(4):
        for t in range(5):
            swap=((q&1)<<1)|(q>>1)
            key=min(q*5+t,swap*5+(-t)%5)
            if key not in reps:reps.append(key)
            labels.append(reps.index(key))
    y=[[model.new_int_var_from_domain(cp_model.Domain.from_values(values),f'y{i}_{r}')
        for r in reps] for i in range(4)]
    for i in range(4):
        for q in range(4):model.add(sum(y[i][labels[q*5+t]] for t in range(5))==z[i][q])
    products={}
    for q in range(4):
        for t in range(5):
            terms=[]
            for i in range(4):
                for a in range(4):
                    for b in range(5):
                        u,v=sorted((labels[a*5+b],labels[(a^q)*5+(b+t)%5]))
                        key=i,u,v
                        if key not in products:
                            pv=sorted({x*w for x in values for w in values})
                            p=model.new_int_var_from_domain(cp_model.Domain.from_values(pv),f'p{i}_{u}_{v}')
                            model.add_multiplication_equality(p,[y[i][u],y[i][v]])
                            products[key]=p
                        terms.append(products[key])
            target=2936 if (q,t)==(0,0) else (-516 if q==0 else 172)
            model.add(sum(terms)==target)
    s=cp_model.CpSolver();s.parameters.max_time_in_seconds=secs;s.parameters.num_workers=1
    status=s.solve(model)
    out=dict(status=s.status_name(status),seconds=s.wall_time,branches=s.num_branches,
             conflicts=s.num_conflicts,z=z,max_abs=max_abs,domain=values,
             scope='fixed S14 Z43 compression fibre; independent all-lag model')
    if status in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        Y=np.array([[s.value(y[i][labels[u]]) for u in range(20)] for i in range(4)])
        for q in range(4):
            for t in range(5):
                corr=sum(int(Y[i,a*5+b])*int(Y[i,(a^q)*5+(b+t)%5]) for i in range(4) for a in range(4) for b in range(5))
                assert corr==(2936 if (q,t)==(0,0) else (-516 if q==0 else 172))
        out['y']=Y.tolist()
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('output');p.add_argument('--secs',type=float,default=180)
    a=p.parse_args();r=json.load(open(a.input))
    out=check(r['z'],a.secs,r.get('max_abs'));out['input']=a.input
    open(a.output,'w').write(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
