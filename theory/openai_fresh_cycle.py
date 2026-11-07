#!/usr/bin/env python3
"""Cycle-exact per-lag projections, with a compiled cardinality Ising DP.

A lag is a union of translation cycles. The DP projects onto binary sequences
with exactly (4|G|-target)/2 domain walls across those cycles. Replicas are per
lag and vertex, not per edge endpoint. A C++17 compiler is required. The C++
source below is compiled into a temporary directory, not into the repository.
"""
import argparse
import ctypes
import itertools
import json
import math
from pathlib import Path
import subprocess
import tempfile
import time

import numpy as np
from openai_fresh_sphere import Model
from openai_fresh_sphere import exact_paf


CPP = r'''
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint>
extern "C" double project(const double* h, int nc, int len, int K, int8_t* out) {
    const double neg = -1e200;
    const int Q = K+1;
    std::vector<double> carry(Q,neg), nextcarry(Q,neg), cur(2*Q), nxt(2*Q);
    std::vector<int8_t> parent((size_t)nc*2*len*Q*2, -1);
    std::vector<int8_t> starts((size_t)nc*Q, -1), ends((size_t)nc*Q,-1);
    auto pi = [=](int c,int f,int j,int k,int s)->size_t {
        return (((((size_t)c*2+f)*len+j)*Q+k)*2+s);
    };
    carry[0]=0;
    for(int c=0;c<nc;c++) {
        std::fill(nextcarry.begin(),nextcarry.end(),neg);
        for(int f=0;f<2;f++) {
            std::fill(cur.begin(),cur.end(),neg);
            for(int k=0;k<=K;k++) if(carry[k]>neg/2)
                cur[2*k+f]=carry[k]+(2*f-1)*h[c*len];
            for(int j=1;j<len;j++) {
                std::fill(nxt.begin(),nxt.end(),neg);
                for(int k=0;k<=K;k++) for(int s=0;s<2;s++) {
                    double a=cur[2*k+s];
                    double b=k?cur[2*(k-1)+1-s]:neg;
                    int prev=(a>=b)?s:1-s;
                    double val=std::max(a,b);
                    if(val>neg/2) {
                        nxt[2*k+s]=val+(2*s-1)*h[c*len+j];
                        parent[pi(c,f,j,k,s)]=prev;
                    }
                }
                cur.swap(nxt);
            }
            for(int k=0;k<=K;k++) for(int s=0;s<2;s++) {
                int kk=k+(s!=f);
                if(kk<=K && cur[2*k+s]>nextcarry[kk]) {
                    nextcarry[kk]=cur[2*k+s];
                    starts[c*Q+kk]=f;
                    ends[c*Q+kk]=s;
                }
            }
        }
        carry.swap(nextcarry);
    }
    double score=carry[K];
    if(score<neg/2) return score;
    int k=K;
    for(int c=nc-1;c>=0;c--) {
        int f=starts[c*Q+k], s=ends[c*Q+k];
        k-=(s!=f);
        for(int j=len-1;j>=1;j--) {
            out[c*len+j]=2*s-1;
            int prev=parent[pi(c,f,j,k,s)];
            k-=(prev!=s);
            s=prev;
        }
        out[c*len]=2*f-1;
    }
    return score;
}
'''


class CycleProjector:
    def __init__(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='openai_fresh_cycle_')
        src = Path(self.tmp.name)/'project.cpp'
        lib = Path(self.tmp.name)/'project.so'
        src.write_text(CPP)
        subprocess.run(['g++', '-std=c++17', '-O3', '-shared', '-fPIC', str(src), '-o', str(lib)], check=True)
        self.lib = ctypes.CDLL(str(lib))
        self.fn = self.lib.project
        self.fn.argtypes = [np.ctypeslib.ndpointer(np.float64, flags='C_CONTIGUOUS'),
                            ctypes.c_int, ctypes.c_int, ctypes.c_int,
                            np.ctypeslib.ndpointer(np.int8, flags='C_CONTIGUOUS')]
        self.fn.restype = ctypes.c_double

    def __call__(self, x, nc, length, k):
        h = np.ascontiguousarray(x, np.float64)
        out = np.empty(h.shape, np.int8)
        val = self.fn(h, nc, length, k, out)
        if val < -1e100:
            raise ValueError('infeasible domain-wall count')
        return out, val


class CycleModel:
    def __init__(self, model, spectral_concur=False, parametric=False):
        self.model = model
        self.spectral_concur = spectral_concur
        if parametric:
            from openai_fresh_parametric import ParametricProjector
            self.dp = ParametricProjector()
        else:
            self.dp = CycleProjector()
        mods = model.mods
        points = list(itertools.product(*(range(m) for m in mods)))
        index = {g:i for i,g in enumerate(points)}
        if model.kind == 'b0':
            vmap = model.labels
            step = len(model.orbits)
        else:
            vmap = (np.arange(20)[:, None]*7 + model.labels).ravel()
            step = 140
        self.vertex = np.stack([vmap+b*step for b in range(4)])
        self.nvars = 4*step
        seen = {tuple(0 for _ in mods)}
        ids = []
        self.factors = []
        cursor = 0
        for lag in points:
            if lag in seen:
                continue
            for u in model.orbits[1]:
                for sign in (-1, 1):
                    a = [(sign*t) % m for t,m in zip(lag,mods)]
                    a[-1] = a[-1]*u % mods[-1]
                    seen.add(tuple(a))
            visited = set()
            cycles = []
            for g in points:
                if g in visited:
                    continue
                cyc = []
                cur = g
                while cur not in visited:
                    visited.add(cur)
                    cyc.append(index[cur])
                    cur = tuple((a+b) % m for a,b,m in zip(cur,lag,mods))
                cycles.append(cyc)
            length = len(cycles[0])
            assert all(len(c)==length for c in cycles)
            order = np.array(cycles).ravel()
            allids = self.vertex[:, order].ravel()
            target = int(model.target[lag])
            assert (4*model.n-target) % 2 == 0
            walls = (4*model.n-target)//2
            nc = 4*len(cycles)
            self.factors.append((cursor, cursor+len(allids), nc, length, walls, lag))
            ids.extend(allids)
            cursor += len(allids)
        self.rows = []
        if model.profile:
            for b,target in enumerate(model.profile):
                v = self.vertex[b]
                self.rows.append((cursor,cursor+len(v),(len(v)+target)//2))
                ids.extend(v)
                cursor += len(v)
        if model.kind == 'k4_863':
            for v,target in zip(self.vertex.reshape(80,43), model.y.ravel()):
                self.rows.append((cursor,cursor+43,(43+int(target))//2))
                ids.extend(v)
                cursor += 43
        self.ids = np.array(ids, np.int32)
        self.counts = np.bincount(self.ids, minlength=self.nvars)
        self.ncopies = cursor
        if spectral_concur:
            # Every full coordinate has the same number of replicas. Thus the
            # induced metric is ordinary full-coordinate Euclidean distance.
            sizes = np.bincount(self.vertex.ravel(), minlength=self.nvars)
            assert np.all(self.counts/sizes == self.counts[0]/sizes[0])

    def pa(self,x):
        a = np.empty_like(x)
        for lo,hi,nc,length,k,lag in self.factors:
            a[lo:hi], _ = self.dp(x[lo:hi],nc,length,k)
        for lo,hi,k in self.rows:
            a[lo:hi] = -1
            if k:
                j = np.argpartition(x[lo:hi],hi-lo-k)[-k:]
                a[lo+j] = 1
        return a

    def means(self,x):
        return np.bincount(self.ids, weights=x, minlength=self.nvars)/self.counts

    def pb(self,x):
        means = self.means(x)
        if self.spectral_concur:
            full = means[self.vertex].reshape((4,)+self.model.mods)
            ff = np.fft.fftn(full, axes=self.model.axes)
            norms = np.sqrt(np.sum(np.abs(ff)**2, axis=0))
            scale = np.divide(self.model.radius,norms,out=np.zeros_like(norms),where=norms>1e-15)
            ff *= scale[None]
            ff[0] = np.where(norms<=1e-15,self.model.radius,ff[0])
            if self.model.profile:
                ff[:,0] = self.model.profile
            elif self.model.kind == 'k4_863':
                ff[...,0] = np.fft.fftn(self.model.y,axes=(1,2,3))
            full = np.fft.ifftn(ff, axes=self.model.axes).real.reshape(4,-1)
            sizes = np.bincount(self.vertex.ravel(),minlength=self.nvars)
            means = np.bincount(self.vertex.ravel(),weights=full.ravel(),minlength=self.nvars)/sizes
            return means[self.ids]
        # Exact weighted projection onto concurrence AND the prescribed linear
        # sums. The metric weights are the occurrence counts of each variable.
        if self.model.profile:
            sizes = np.bincount(self.model.labels)
            for b, target in enumerate(self.model.profile):
                lo = b*len(sizes)
                hi = lo+len(sizes)
                direction = sizes/self.counts[lo:hi]
                lam = (target-sizes@means[lo:hi])/(sizes@direction)
                means[lo:hi] += lam*direction
        elif self.model.kind == 'k4_863':
            sizes = np.array([1]+[7]*6)
            for cell, target in enumerate(self.model.y.ravel()):
                lo=7*cell
                hi=lo+7
                direction=sizes/self.counts[lo:hi]
                lam=(target-sizes@means[lo:hi])/(sizes@direction)
                means[lo:hi]+=lam*direction
        return means[self.ids]

    def assert_local(self,a):
        assert np.all(np.abs(a)==1)
        for lo,hi,nc,length,k,lag in self.factors:
            aa = a[lo:hi].reshape(nc,length)
            assert np.sum(aa != np.roll(aa,1,axis=-1)) == k


def controls():
    dp = CycleProjector()
    rng = np.random.default_rng(5)
    ncases = 0
    for nc,length in [(1,3),(2,3),(2,4),(4,3)]:
        allx = np.array(list(itertools.product((-1,1),repeat=nc*length)),np.int8)
        xx = allx.reshape(-1,nc,length)
        walls = np.sum(xx != np.roll(xx,1,axis=-1),axis=(1,2))
        for rep in range(3):
            h = rng.normal(size=nc*length)
            for k in sorted(set(walls.tolist())):
                a,val = dp(h,nc,length,k)
                optimum = np.max(allx[walls==k]@h)
                assert abs(val-optimum)<1e-9
                assert abs(a@h-optimum)<1e-9
                aa=a.reshape(nc,length)
                assert np.sum(aa!=np.roll(aa,1,axis=-1))==k
                ncases += 1
    return dict(exhaustive_projection_cases=ncases,passed=True)


def run(args):
    model = Model(args.kind,args.p,args.d,args.profile,args.cert)
    r = CycleModel(model,args.spectral_concur,args.parametric)
    rng = np.random.default_rng(args.seed)
    t0 = time.perf_counter()
    records=[]
    for trial in range(args.trials):
        x=rng.normal(size=r.ncopies)
        r.assert_local(r.pa(x))
        best=math.inf
        hit=None
        outside_profile=None
        trial_start=time.perf_counter()
        for it in range(1,args.iters+1):
            a=r.pa(x)
            b=r.pb(2*a-x)
            x+=args.beta*(b-a)
            if it%args.check==0:
                z=np.where(r.means(a)>=0,1,-1)
                seq=z[r.vertex].reshape((4,)+model.mods)
                f=np.fft.fftn(seq,axes=model.axes)
                err=np.rint(np.fft.ifftn(np.sum(np.abs(f)**2,axis=0)).real-model.target).astype(np.int64)
                energy=int(np.sum(err**2))
                best=min(best,energy)
                profile_ok=model.profile is None or np.array_equal(seq.sum(-1),model.profile)
                fibre_ok=model.kind!='k4_863' or np.array_equal(seq.sum(-1),model.y)
                if energy==0 and not (profile_ok and fibre_ok) and outside_profile is None:
                    assert np.array_equal(exact_paf(seq,model.mods),model.target)
                    outside_profile=dict(rowsums=seq.reshape(4,-1).sum(1).tolist(),
                                         sequences=seq.reshape(4,-1).tolist(),iterations=it)
                if energy==0 and profile_ok and fibre_ok:
                    assert model.integer_check(seq)
                    hit=seq.reshape(4,-1).tolist()
                    break
                if time.perf_counter()-t0>=args.seconds:
                    break
        records.append(dict(trial=trial,iterations=it,seconds=time.perf_counter()-trial_start,
                            best_integer_paf_energy=best,solved=hit is not None,sequences=hit,
                            outside_requested_sums=outside_profile))
        print(json.dumps({k:v for k,v in records[-1].items()
                          if k not in ('sequences','outside_requested_sums')}
                         | {'other_core_recorded':outside_profile is not None}),flush=True)
        if time.perf_counter()-t0>=args.seconds:
            break
    result=dict(algorithm='cycle-exact-lag-projection',parametric_integer_fields=args.parametric,affine_sums=True,spectral_concur=args.spectral_concur,kind=args.kind,p=args.p,d=args.d,
                profile=args.profile,beta=args.beta,seed=args.seed,trials_requested=args.trials,
                variables=r.nvars,copies=r.ncopies,lag_factors=len(r.factors),
                seconds=time.perf_counter()-t0,solved=sum(r['solved'] for r in records),
                trials=records,certificate='direct integer PAF only; no matrix claim')
    if args.out:
        Path(args.out).parent.mkdir(parents=True,exist_ok=True)
        Path(args.out).write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--kind',choices=['b0','k4_863'],default='b0')
    p.add_argument('--p',type=int,default=43)
    p.add_argument('--d',type=int,default=3)
    p.add_argument('--profile',type=lambda x:list(map(int,x.split(','))))
    p.add_argument('--cert',default='data/h3452_stage1.json')
    p.add_argument('--trials',type=int,default=4)
    p.add_argument('--seed',type=int,default=0)
    p.add_argument('--beta',type=float,default=.5)
    p.add_argument('--iters',type=int,default=2000)
    p.add_argument('--seconds',type=float,default=240)
    p.add_argument('--check',type=int,default=10)
    p.add_argument('--controls',action='store_true')
    p.add_argument('--spectral-concur',action='store_true')
    p.add_argument('--parametric',action='store_true')
    p.add_argument('--out')
    args=p.parse_args()
    if args.controls:
        result=controls()
        if args.out:
            Path(args.out).parent.mkdir(parents=True,exist_ok=True)
            Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result))
    else:
        run(args)
