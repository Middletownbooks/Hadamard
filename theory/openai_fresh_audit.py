#!/usr/bin/env python3
"""Independent integer source audit, arithmetic gates, and projector benchmarks."""
import argparse
import gzip
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np


def seqs(p):
    root=Path('forensics/originals')
    if p in (223,283,487):
        lines=(root/f'seq_v{p}_g.txt').read_text().splitlines()[1:5]
        rows=[line.split('=',1)[1] for line in lines]
    else:
        matches=list(root.glob(f'matrix_H*_v{p}_*.txt.gz'))
        assert len(matches)==1
        border=4 if p in (167,179) else 0
        m=p-1 if border else p
        with gzip.open(matches[0],'rt') as f:
            header=f.readline()
            for _ in range(border+1):
                row=f.readline().strip()
        rows=[row[border+b*m:border+(b+1)*m] for b in range(4)]
    x=np.array([[1 if c=='+' else -1 for c in row] for row in rows],np.int64)
    m=x.shape[1]
    x[1:]=x[1:,(-1-np.arange(m))%m]
    return x


def pafs(x):
    return np.array([np.sum(x*np.roll(x,-t,axis=1),axis=1,dtype=np.int64)
                     for t in range(x.shape[1])],np.int64).T


def source_audit():
    flat=[]; half=[]
    for p in (223,283,311,487,491):
        x=seqs(p); pa=pafs(x)
        assert np.array_equal(pa.sum(0),np.array([4*p]+[0]*(p-1)))
        energy=(pa[:,1:]**2).sum(1)
        rows=x.sum(1)
        numerator=p*energy-(rows**2-p)**2
        denominator=(p-1)*p*p
        flat.append(dict(p=p,integer_core_verified=True,rowsums=rows.tolist(),
                         phi_numerator=numerator.tolist(),phi_denominator=denominator,
                         phi=(numerator/denominator).tolist()))
    for p in (167,179):
        x=seqs(p);m=p-1;h=m//2;pa=pafs(x)
        assert np.array_equal(pa.sum(0),np.array([4*m]+[-4]*(m-1)))
        sig=x[:,:h]*x[:,h:]
        changes=(sig!=np.roll(sig,1,axis=1)).sum(1)
        fits=[]
        for row in sig:
            opts=[]
            for boundary in range(1,h):
                pattern=np.where(np.arange(h)>=boundary,1,-1)
                ex=np.flatnonzero(row!=pattern).tolist()
                opts.append((len(ex),boundary,ex))
            fits.append(min(opts))
        c=(x[:,:h]+x[:,h:])//2;d=(x[:,:h]-x[:,h:])//2
        pc=pafs(c).sum(0)
        nd=[]
        for t in range(h):
            pos=np.arange(h)+t
            nd.append(int(np.sum(d*d[:,pos%h]*np.where(pos<h,1,-1))))
        assert np.array_equal(pc,np.array([2*h-1]+[-2]*(h-1)))
        assert nd==[2*h+1]+[0]*(h-1)
        half.append(dict(p=p,sign_changes=changes.tolist(),anchored_arc_fits=fits,
                         exceptions=sum(f[0] for f in fits),
                         periodic_and_negacyclic_integer_checks=True))
    return dict(flatness=flat,halfshift=half)


def row_gate(p,d):
    sums=sorted({abs(e+d*(2*k-(p-1)//d)) for e in (-1,1)
                 for k in range((p-1)//d+1)} & set(range(math.isqrt(4*p)+1)))
    profiles=[list(r) for r in itertools.combinations_with_replacement(sums,4)
              if sum(a*a for a in r)==4*p]
    return dict(p=p,d=d,allowed_absolute_rowsums=sums,profiles=profiles)


def benchmark():
    from openai_fresh_sphere import Model
    from openai_fresh_cycle import CycleModel
    from openai_fresh_parametric import ParametricProjector
    rng=np.random.default_rng(991)
    out=[]
    for kind,p,d,profile in [('b0',1103,19,row_gate(1103,19)['profiles'][0]),
                            ('k4_863',863,7,None)]:
        model=Model(kind,p,d,profile,'data/h3452_stage1.json')
        r=CycleModel(model,True,False)
        fast=ParametricProjector()
        x=rng.normal(size=r.ncopies)
        # Same rounded integer fields for the two projection implementations.
        x=np.rint(x*1024.)
        t0=time.perf_counter();a=r.pa(x);slow_seconds=time.perf_counter()-t0
        r.assert_local(a)
        t0=time.perf_counter();b=np.empty_like(a);calls=[]
        for lo,hi,nc,length,k,lag in r.factors:
            b[lo:hi],score=fast.integer(x[lo:hi].astype(np.int64),nc,length,k)
            assert score==int(a[lo:hi]@x[lo:hi])
            calls.append(fast.last_calls)
        for lo,hi,k in r.rows:
            b[lo:hi]=a[lo:hi]
        fast_seconds=time.perf_counter()-t0
        r.assert_local(b)
        dp_updates=sum(4*nc*(length-1)*(k+1) for lo,hi,nc,length,k,lag in r.factors)
        fast_updates=sum(4*nc*(length-1)*ncall for (lo,hi,nc,length,k,lag),ncall in zip(r.factors,calls))
        # Complete relaxed iteration, with the real-field quantization used by
        # the search. Timings exclude model construction and compilation.
        r.dp=fast
        x=rng.normal(size=r.ncopies)
        t0=time.perf_counter()
        for _ in range(3):
            aa=r.pa(x);bb=r.pb(2*aa-x);x+=bb-aa
        iteration_seconds=(time.perf_counter()-t0)/3
        out.append(dict(kind=kind,p=p,d=d,profile=profile,variables=r.nvars,
                        copies=r.ncopies,lag_factors=len(r.factors),
                        quadratic_DP_state_updates=dp_updates,
                        parametric_state_updates=fast_updates,
                        integer_field_oracle_calls=calls,
                        quadratic_projection_seconds=slow_seconds,
                        parametric_projection_seconds=fast_seconds,
                        quantized_iteration_seconds=iteration_seconds,
                        every_factor_optimum_agrees=True))
    return out


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--benchmark',action='store_true')
    ap.add_argument('--out',default='logs/openai_fresh/audit.json')
    args=ap.parse_args()
    t0=time.perf_counter()
    result=dict(sources=source_audit(),row_gates=[row_gate(599,13),row_gate(599,23),row_gate(1103,19)])
    if args.benchmark:
        result['benchmark']=benchmark()
    result['seconds']=time.perf_counter()-t0
    Path(args.out).parent.mkdir(parents=True,exist_ok=True)
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
