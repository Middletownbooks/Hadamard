#!/usr/bin/env python3
"""Parametric exact projection onto a prescribed number of cyclic domain walls.

Integer fields, exact int64 arithmetic. The concavity/exchange proof and the
bound for quantizing real fields are in report/openai_fresh_answer.md.
The independent quadratic DP oracle is openai_fresh_cycle.CycleProjector.
"""
import argparse
import ctypes
import itertools
import json
from pathlib import Path
import subprocess
import tempfile
import time

import numpy as np


CPP = r'''
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
using ll = int64_t;
const ll BAD = -(1LL<<61);
struct State { ll value; int count; };
bool better(State a, State b, int mode) {
    return a.value>b.value || (a.value==b.value && mode*a.count>mode*b.count);
}
// Maximize 2 sum h_i x_i - lambda * (# walls), with minimum/maximum
// wall count among ties. Returns one optimizer on each translation cycle.
ll oracle(const ll* h, int nc, int len, ll lambda, int mode, int8_t* out, int& total) {
    ll score=0; total=0;
    std::vector<int8_t> parent(2*len*2), temp(len), chosen(len);
    for(int c=0;c<nc;c++) {
        State best={BAD,0};
        for(int first=0;first<2;first++) {
            State cur[2]={{BAD,0},{BAD,0}};
            cur[first]={2*h[c*len]*(2*first-1),0};
            for(int j=1;j<len;j++) {
                State nxt[2];
                for(int s=0;s<2;s++) {
                    State keep=cur[s];
                    State flip={cur[1-s].value-lambda,cur[1-s].count+1};
                    int prev=better(keep,flip,mode)?s:1-s;
                    State z=(prev==s)?keep:flip;
                    z.value+=2*h[c*len+j]*(2*s-1);
                    nxt[s]=z;
                    parent[(first*len+j)*2+s]=prev;
                }
                cur[0]=nxt[0];cur[1]=nxt[1];
            }
            for(int s=0;s<2;s++) {
                State z={cur[s].value-lambda*(s!=first),cur[s].count+(s!=first)};
                if(better(z,best,mode)) {
                    best=z;
                    int bit=s;
                    for(int j=len-1;j>=1;j--) {
                        temp[j]=2*bit-1;
                        bit=parent[(first*len+j)*2+bit];
                    }
                    temp[0]=2*first-1;
                    chosen=temp;
                }
            }
        }
        score+=best.value; total+=best.count;
        std::copy(chosen.begin(),chosen.end(),out+c*len);
    }
    return score;
}
extern "C" int project_int(const ll* h, int nc, int len, int walls, int8_t* out,
                            ll* scoreout, int* callsout) {
    int n=nc*len;
    if(walls<0 || walls>nc*2*(len/2) || walls%2) return 1;
    ll sum=0;
    for(int j=0;j<n;j++) sum+=std::llabs(h[j]);
    // Keep every accumulated Lagrangian score far from both int64 overflow
    // and the unreachable-state sentinel, including extreme lambda values.
    if(sum>(1LL<<58)/(2LL*n+4)) return 2;
    ll lo=-2*sum-1, hi=2*sum+1;
    std::vector<int8_t> a(n),b(n);
    int na=0,nb=0,calls=0;
    bool found=false;
    ll lam=0, optimum=0;
    while(lo<=hi) {
        lam=lo+(hi-lo)/2;
        ll va=oracle(h,nc,len,lam,-1,a.data(),na);
        ll vb=oracle(h,nc,len,lam, 1,b.data(),nb);
        calls+=2;
        if(va!=vb || na>nb) return 3;
        if(na<=walls && walls<=nb) {found=true;optimum=va;break;}
        if(walls<na) lo=lam+1;
        else hi=lam-1;
    }
    if(!found) return 4;
    // Exchange disagreement components. Every +2-wall component can be
    // swapped independently between two Lagrangian optimizers, preserving
    // the total score. Both resulting sequences are still optimizers.
    for(int c=0;c<nc && na<walls;c++) {
        int off=c*len, common=-1;
        for(int j=0;j<len;j++) if(a[off+j]==b[off+j]) {common=j;break;}
        if(common<0) continue; // whole-cycle disagreement has wall delta 0
        int t=1;
        while(t<len && na<walls) {
            int j=(common+t)%len;
            if(a[off+j]==b[off+j]) {t++;continue;}
            int start=t;
            while(t<len && a[off+(common+t)%len]!=b[off+(common+t)%len]) t++;
            int end=t-1;
            int left=(common+start-1)%len;
            int first=(common+start)%len;
            int last=(common+end)%len;
            int right=(common+end+1)%len;
            int old=(a[off+left]!=a[off+first])+(a[off+last]!=a[off+right]);
            int alt=(a[off+left]!=b[off+first])+(b[off+last]!=a[off+right]);
            if(alt-old==2) {
                for(int k=start;k<=end;k++) std::swap(a[off+(common+k)%len],b[off+(common+k)%len]);
                na+=2;
            }
        }
    }
    if(na!=walls) return 5;
    ll dot=0;int check=0;
    for(int c=0;c<nc;c++) for(int j=0;j<len;j++) {
        int i=c*len+j;
        dot+=h[i]*a[i];
        check+=(a[i]!=a[c*len+(j+1)%len]);
    }
    if(check!=walls || 2*dot-lam*walls!=optimum) return 6;
    std::copy(a.begin(),a.end(),out);
    *scoreout=dot;*callsout=calls;
    return 0;
}
'''


class ParametricProjector:
    def __init__(self, quantization=1048576.):
        self.quantization=float(quantization)
        self.tmp=tempfile.TemporaryDirectory(prefix='openai_fresh_parametric_')
        src=Path(self.tmp.name)/'project.cpp'
        lib=Path(self.tmp.name)/'project.so'
        src.write_text(CPP)
        subprocess.run(['g++','-std=c++17','-O3','-shared','-fPIC',str(src),'-o',str(lib)],check=True)
        self.lib=ctypes.CDLL(str(lib))
        self.fn=self.lib.project_int
        self.fn.argtypes=[np.ctypeslib.ndpointer(np.int64,flags='C_CONTIGUOUS'),
                          ctypes.c_int,ctypes.c_int,ctypes.c_int,
                          np.ctypeslib.ndpointer(np.int8,flags='C_CONTIGUOUS'),
                          ctypes.POINTER(ctypes.c_int64),ctypes.POINTER(ctypes.c_int)]
        self.fn.restype=ctypes.c_int
        self.last_calls=0

    def integer(self,h,nc,length,k):
        h=np.ascontiguousarray(h,np.int64)
        n=nc*length
        if (nc<1 or length<1 or h.ndim!=1 or len(h)!=n
                or np.any(h==np.iinfo(np.int64).min)
                or int(np.max(np.abs(h),initial=0))*n>2**58//(2*n+4)):
            raise ValueError('invalid input size or unsafe integer range')
        out=np.empty(h.shape,np.int8)
        score=ctypes.c_int64()
        calls=ctypes.c_int()
        status=self.fn(h,nc,length,k,out,ctypes.byref(score),ctypes.byref(calls))
        if status:
            raise ValueError(f'parametric projector status {status}')
        self.last_calls=calls.value
        return out,score.value

    def __call__(self,x,nc,length,k):
        scale=self.quantization/max(1.,float(np.max(np.abs(x)))/1024.)
        h=np.rint(np.asarray(x)*scale).astype(np.int64)
        out,score=self.integer(h,nc,length,k)
        # Exact for h; only an approximate Euclidean projection for real x.
        return out,float(np.asarray(x)@out)


def controls(seed=19):
    from openai_fresh_cycle import CycleProjector
    dp=CycleProjector()
    fast=ParametricProjector()
    rng=np.random.default_rng(seed)
    brute_cases=0; dp_cases=0; concavity_cases=0; maximum_calls=0
    for nc,length in [(1,3),(1,4),(1,5),(2,3),(2,4),(3,3),(4,3)]:
        allx=np.array(list(itertools.product((-1,1),repeat=nc*length)),np.int8)
        xx=allx.reshape(-1,nc,length)
        walls=np.sum(xx!=np.roll(xx,1,axis=-1),axis=(1,2))
        for rep in range(12):
            h=rng.integers(-4,5,nc*length,dtype=np.int64)
            vals=allx@h
            f=[]
            for k in sorted(set(walls.tolist())):
                a,score=fast.integer(h,nc,length,k)
                optimum=int(np.max(vals[walls==k]))
                assert score==optimum==int(a@h)
                assert np.sum(a.reshape(nc,length)!=np.roll(a.reshape(nc,length),1,axis=-1))==k
                f.append(optimum)
                maximum_calls=max(maximum_calls,fast.last_calls)
                brute_cases+=1
            assert all(2*f[k]>=f[k-1]+f[k+1] for k in range(1,len(f)-1))
            concavity_cases+=1
    for nc,length in [(4,17),(4,43),(4,61),(16,5),(8,43),(4,131)]:
        for rep in range(10):
            h=rng.integers(-10000,10001,nc*length,dtype=np.int64)
            k=2*int(rng.integers(nc*(length//2)+1))
            a,score=fast.integer(h,nc,length,k)
            b,reference=dp(h.astype(float),nc,length,k)
            assert score==int(reference)==int(a@h)==int(b@h)
            maximum_calls=max(maximum_calls,fast.last_calls)
            dp_cases+=1
    # A dimensionally large field must be rejected before score arithmetic
    # can overflow, even when its unweighted field sum would fit int64.
    try:
        fast.integer(np.full(100,2**48,dtype=np.int64),4,25,100)
    except ValueError:
        range_guard=True
    else:
        raise AssertionError('unsafe integer field was accepted')
    return dict(exhaustive_integer_projections=brute_cases,
                integer_concavity_tables=concavity_cases,
                independent_quadratic_DP_cases=dp_cases,
                maximum_oracle_calls=maximum_calls,unsafe_range_rejected=range_guard,
                passed=True,seed=seed)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--controls',action='store_true')
    ap.add_argument('--out')
    args=ap.parse_args()
    t0=time.perf_counter()
    result=controls()
    result['seconds']=time.perf_counter()-t0
    if args.out:
        Path(args.out).parent.mkdir(parents=True,exist_ok=True)
        Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
