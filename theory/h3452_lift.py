"""Bounded lift search inside a certified S7 Z43 compression fibre.

Each of 80 cells has six size-seven cyclotomic orbits. Its compression
fixes the sign at zero and the number of positive cyclotomic signs.
Swapping opposite signs in one cell preserves every compressed equation.
Exact quadratic PAF updates; no approximate Fourier acceptance checks.
Failure is not an obstruction. A zero residual is rechecked on all 860 lags.
"""
import argparse,json,time,math
from pathlib import Path
import numpy as np
from numba import njit
from h3452_compress import structure,verify


def coefficients(st):
    r=st.r[0];K=np.zeros((st.ncls,r,r),np.int64)
    for j,g in enumerate(st.cls_rep):
        np.add.at(K[j],(st.orb[0],st.orb[0][st.G.add[:,g]]),1)
    return K,K+K.transpose(0,2,1)


@njit
def evaluate(eps,K):
    out=np.zeros(K.shape[0],np.int64)
    for j in range(K.shape[0]):
        for a in range(len(eps)):
            for b in range(len(eps)):out[j]+=K[j,a,b]*eps[a]*eps[b]
    return out


@njit
def delta(eps,fields,K,S,a,b):
    out=np.empty(K.shape[0],np.int64)
    for j in range(len(out)):
        out[j]=-2*eps[a]*fields[j,a]-2*eps[b]*fields[j,b]+4*K[j,a,a]+4*K[j,b,b]+4*eps[a]*eps[b]*S[j,a,b]
    return out


@njit
def search(y,K,S,target,seed,restarts,steps):
    np.random.seed(seed);nc=K.shape[0];r=K.shape[1]
    best=1e100;best_eps=np.zeros((4,r),np.int64);attempts=0
    for restart in range(restarts):
        eps=np.ones((4,r),np.int64)
        for i in range(4):
            for cell in range(20):
                yy=y[i,cell];a=1 if (yy-1)%14==0 else -1
                k=((yy-a)//7+6)//2
                eps[i,7*cell]=a
                order=np.random.permutation(6)
                for j in range(6):eps[i,7*cell+1+order[j]]=1 if j<k else -1
        fields=np.zeros((4,nc,r),np.int64)
        res=-target.copy()
        for i in range(4):
            res+=evaluate(eps[i],K)
            for j in range(nc):
                for a in range(r):
                    for b in range(r):fields[i,j,a]+=S[j,a,b]*eps[i,b]
        cost=(res*res).sum()/256.0
        for step in range(steps):
            i=np.random.randint(4);cell=np.random.randint(20)
            a=cell*7+1+np.random.randint(6);b=cell*7+1+np.random.randint(6)
            if eps[i,a]==eps[i,b]:continue
            dd=delta(eps[i],fields[i],K,S,a,b);rr=res+dd
            cc=(rr*rr).sum()/256.0
            temp=0.25+100.0*(1-step/steps)**3
            if cc<=cost or np.random.random()<np.exp((cost-cc)/temp):
                da=-2*eps[i,a];db=-2*eps[i,b]
                for j in range(nc):
                    for o in range(r):fields[i,j,o]+=S[j,o,a]*da+S[j,o,b]*db
                eps[i,a]=-eps[i,a];eps[i,b]=-eps[i,b];res=rr;cost=cc
            attempts+=1
            if cost<best:
                best=cost;best_eps=eps.copy()
                if best==0:return best,best_eps,attempts
    return best,best_eps,attempts


def validate_updates(st,K,S):
    rng=np.random.default_rng(3452)
    for _ in range(100):
        eps=rng.choice([-1,1],st.r[0]).astype(np.int64)
        fields=np.einsum('jab,b->ja',S,eps)
        a,b=rng.choice(st.r[0],2,replace=False)
        dd=delta(eps,fields,K,S,a,b)
        changed=eps.copy();changed[[a,b]]*=-1
        direct0=st.paf_classes(st.expand(0,eps));direct1=st.paf_classes(st.expand(0,changed))
        assert np.array_equal(evaluate(eps,K),direct0)
        assert np.array_equal(dd,direct1-direct0)


def main():
    p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('--out',required=True)
    p.add_argument('--restarts',type=int,default=10);p.add_argument('--steps',type=int,default=100000)
    p.add_argument('--seed',type=int,default=1);a=p.parse_args()
    cert=json.load(open(a.input));st,rec=structure('s7');ng=[[0,0,0,1]]
    y=np.array(cert['y'],np.int64);verify(y,rec,ng)
    # Check the assumed cell/orbit layout explicitly rather than relying on numbering.
    assert st.r[0]==140
    for cell in range(20):
        assert st.orb[0][43*cell]==7*cell
        assert set(st.orb[0][43*cell+1:43*(cell+1)])==set(range(7*cell+1,7*cell+7))
    K,S=coefficients(st);validate_updates(st,K,S)
    t=time.monotonic();cost,eps,attempts=search(y,K,S,st.cls_target,a.seed,a.restarts,a.steps)
    xs=np.array([st.expand(i,eps[i]) for i in range(4)])
    assert np.array_equal(xs.reshape(4,20,43).sum(2),y)
    # Final residual from direct sums, independent of the incremental updates.
    full=sum(st.paf(row) for row in xs)-st.F
    assert cost==float(np.sum(full[st.cls_rep]**2))/256
    entropy=sum(math.log2(math.comb(6,((int(v)-(1 if (v-1)%14==0 else -1))//7+6)//2)) for v in y.flat)
    out=dict(stage='full core lift; no matrix unless core passes and border is assembled',
             order=3452,input=a.input,seed=a.seed,restarts=a.restarts,steps=a.steps,
             attempts=attempts,seconds=time.monotonic()-t,cost=cost,
             status='CORE_FOUND' if not full.any() else 'NO_CORE_FOUND',
             max_abs_residual=int(abs(full).max()),nonzero_lags=int(np.count_nonzero(full)),
             class_residual=full[st.cls_rep].tolist(),orbit_signs=eps.tolist(),
             log2_fibre_size=entropy,delta_tests=100,compression_preserved=True)
    dest=Path(a.out);dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('orbit_signs','class_residual')}),flush=True)


if __name__=='__main__':main()
