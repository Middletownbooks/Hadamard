"""Whole-cell tabu search inside the certified Z43 fibre. Failure is inconclusive."""
import argparse,json,time
from pathlib import Path
import numpy as np
from numba import njit
from h3452_compress import structure,verify
from h3452_lift import coefficients,evaluate


@njit
def rebuild(eps,K,S,target):
    fields=np.zeros((4,K.shape[0],K.shape[1]),np.int64);res=-target.copy()
    for i in range(4):
        res+=evaluate(eps[i],K)
        for j in range(K.shape[0]):
            for a in range(K.shape[1]):
                for b in range(K.shape[1]):fields[i,j,a]+=S[j,a,b]*eps[i,b]
    return fields,res


@njit
def search(start,K,S,target,patterns,seed,steps,tenure):
    np.random.seed(seed);eps=start.copy();nc=K.shape[0]
    fields,res=rebuild(eps,K,S,target)
    cost=(res*res).sum();best=cost;best_eps=eps.copy()
    taboo=np.zeros((4,20,64),np.int64);codes=np.zeros((4,20),np.int64)
    count=np.zeros((4,20),np.int64)
    for i in range(4):
        for c in range(20):
            for a in range(6):
                if eps[i,7*c+1+a]==1:codes[i,c]+=1<<a;count[i,c]+=1
    evaluated=0;lastbest=0
    for step in range(steps):
        bestmove=1<<60;bi=-1;bc=-1;bp=-1
        dd=np.empty(nc,np.int64);bestdelta=np.empty(nc,np.int64)
        change=np.empty(6,np.int64)
        # Scan every cardinality-preserving replacement, not only two-sign swaps.
        for i in range(4):
            for c in range(20):
                base=7*c+1
                for p in range(64):
                    if p==codes[i,c]:continue
                    npos=0
                    for a in range(6):npos+=patterns[p,a]==1
                    if npos!=count[i,c]:continue
                    for a in range(6):change[a]=patterns[p,a]-eps[i,base+a]
                    cc=0
                    for j in range(nc):
                        v=0
                        for a in range(6):v+=change[a]*fields[i,j,base+a]
                        # Same-cell quadratic terms are supported only at h=0.
                        if j<4:
                            for a in range(6):
                                for b in range(6):v+=change[a]*K[j,base+a,base+b]*change[b]
                        dd[j]=v;cc+=(res[j]+v)**2
                    evaluated+=1
                    if taboo[i,c,p]>step and cc>=best:continue
                    if cc<bestmove or (cc==bestmove and np.random.random()<0.01):
                        bestmove=cc;bi=i;bc=c;bp=p;bestdelta=dd.copy()
        if bi<0:break
        old=codes[bi,bc];taboo[bi,bc,old]=step+tenure+np.random.randint(5)
        base=7*bc+1
        for a in range(6):change[a]=patterns[bp,a]-eps[bi,base+a]
        for j in range(nc):
            for o in range(K.shape[1]):
                for a in range(6):fields[bi,j,o]+=S[j,o,base+a]*change[a]
        eps[bi,base:base+6]=patterns[bp];codes[bi,bc]=bp;res+=bestdelta;cost=bestmove
        if cost<best:
            best=cost;best_eps=eps.copy();lastbest=step
            if best==0:break
        if step-lastbest>=400:
            # Diversify with eight fibre-preserving random cell replacements.
            for _ in range(8):
                i=np.random.randint(4);c=np.random.randint(20)
                pp=np.random.permutation(eps[i,7*c+1:7*c+7])
                eps[i,7*c+1:7*c+7]=pp;code=0
                for a in range(6):
                    if pp[a]==1:code+=1<<a
                codes[i,c]=code
            fields,res=rebuild(eps,K,S,target);cost=(res*res).sum();lastbest=step
    return best,best_eps,evaluated


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',default='logs/h3452/lift_seed1.json')
    p.add_argument('--out',required=True);p.add_argument('--seed',type=int,default=1)
    p.add_argument('--steps',type=int,default=1000);p.add_argument('--batches',type=int,default=4)
    p.add_argument('--tenure',type=int,default=9);a=p.parse_args()
    st,rec=structure('s7');y=np.array(json.load(open('data/h3452_stage1.json'))['y'])
    verify(y,rec,[[0,0,0,1]])
    eps=np.array(json.load(open(a.input))['orbit_signs'],np.int64)
    K,S=coefficients(st)
    for c in range(20):assert not K[4:,7*c+1:7*c+7,7*c+1:7*c+7].any()
    patterns=np.array([[1 if code&(1<<j) else -1 for j in range(6)] for code in range(64)],np.int64)
    runs=[]
    for batch in range(a.batches):
        t=time.monotonic();score,eps,neval=search(eps,K,S,st.cls_target,patterns,a.seed+batch,a.steps,a.tenure)
        xs=np.array([st.expand(i,e) for i,e in enumerate(eps)])
        assert np.array_equal(xs.reshape(4,20,43).sum(2),y)
        full=sum(st.paf(x) for x in xs)-st.F
        assert np.sum(full[st.cls_rep]**2)==score
        runs.append(dict(seed=a.seed+batch,steps=a.steps,tenure=a.tenure,evaluated=neval,
                         seconds=time.monotonic()-t,cost=score/256))
        out=dict(order=3452,stage='full core lift',input=a.input,status='CORE_FOUND' if score==0 else 'NO_CORE_FOUND',
                 orbit_signs=eps.tolist(),class_residual=full[st.cls_rep].tolist(),cost=score/256,
                 max_abs_residual=int(abs(full).max()),nonzero_lags=int(np.count_nonzero(full)),
                 compression_preserved=True,all_lags_rechecked=True,runs=runs)
        dest=Path(a.out);dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps(runs[-1]),flush=True)
        if score==0:break


if __name__=='__main__':main()
