"""Rectangle-switch tabu search preserving both certified complementary projections.

The original Z43 compression is unchanged. An additional H-compression onto
Z43 is fixed to prime_gate.json. Feasibility of both projections is not a core.
"""
import argparse,itertools,json,time
from pathlib import Path
import numpy as np
from numba import njit
from ortools.sat.python import cp_model
from h3452_compress import structure,verify
from h3452_lift import coefficients
from h3452_tabu import rebuild


def seed(y,prime,hint):
    m=cp_model.CpModel();b={};terms=[]
    for i in range(4):
        for c in range(20):
            for j in range(6):
                b[i,c,j]=m.new_bool_var(f'b{i}_{c}_{j}')
                terms.append(int(hint[i,7*c+1+j])*b[i,c,j])
            e=1 if (y[i,c]-1)%14==0 else -1
            m.add(sum(b[i,c,j] for j in range(6))==int(((y[i,c]-e)//7+6)//2))
        for j in range(6):m.add(sum(b[i,c,j] for c in range(20))==int(10+prime[i,j+1]))
    m.maximize(sum(terms));s=cp_model.CpSolver();s.parameters.max_time_in_seconds=10;s.parameters.num_workers=1
    status=s.solve(m);assert status in (cp_model.OPTIMAL,cp_model.FEASIBLE)
    eps=hint.copy()
    for (i,c,j),v in b.items():eps[i,7*c+1+j]=2*s.value(v)-1
    return eps,s.status_name(status)


def moves(K):
    ids=[]
    for a,b in itertools.combinations(range(20),2):
        for u,v in itertools.combinations(range(6),2):ids.append([7*a+1+u,7*a+1+v,7*b+1+u,7*b+1+v])
    ids=np.array(ids,np.int64);quad=np.zeros((len(ids),K.shape[0]),np.int64)
    signs=np.array([-2,2,2,-2])
    for n,ix in enumerate(ids):quad[n]=np.einsum('a,jab,b->j',signs,K[:,ix[:,None],ix[None,:]],signs)
    return ids,quad


@njit
def search(start,K,S,target,ids,quad,seednum,steps,tenure):
    np.random.seed(seednum);eps=start.copy();fields,res=rebuild(eps,K,S,target)
    best=(res*res).sum();best_eps=eps.copy();taboo=np.zeros((4,len(ids)),np.int64)
    evaluated=0;lastbest=0
    for step in range(steps):
        bm=1<<60;bi=-1;bn=-1;bd=np.zeros(len(res),np.int64)
        for i in range(4):
            for n in range(len(ids)):
                a,b,c,d=ids[n];sg=eps[i,a]
                if eps[i,b]!=-sg or eps[i,c]!=-sg or eps[i,d]!=sg:continue
                cc=0;dd=np.empty(len(res),np.int64)
                for j in range(len(res)):
                    v=2*sg*(-fields[i,j,a]+fields[i,j,b]+fields[i,j,c]-fields[i,j,d])+quad[n,j]
                    dd[j]=v;cc+=(res[j]+v)**2
                evaluated+=1
                if taboo[i,n]>step and cc>=best:continue
                if cc<bm or (cc==bm and np.random.random()<0.01):bm=cc;bi=i;bn=n;bd=dd
        if bi<0:break
        a,b,c,d=ids[bn];sg=eps[bi,a]
        for j in range(len(res)):
            for o in range(K.shape[1]):fields[bi,j,o]+=2*sg*(-S[j,o,a]+S[j,o,b]+S[j,o,c]-S[j,o,d])
        for ix in ids[bn]:eps[bi,ix]*=-1
        res+=bd;taboo[bi,bn]=step+tenure+np.random.randint(5)
        if bm<best:
            best=bm;best_eps=eps.copy();lastbest=step
            if best==0:break
        if step-lastbest>=400:
            for _ in range(8):
                for attempt in range(100):
                    i=np.random.randint(4);n=np.random.randint(len(ids));a,b,c,d=ids[n];sg=eps[i,a]
                    if eps[i,b]==-sg and eps[i,c]==-sg and eps[i,d]==sg:
                        for ix in ids[n]:eps[i,ix]*=-1
                        break
            fields,res=rebuild(eps,K,S,target);lastbest=step
    return best,best_eps,evaluated


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',default='logs/h3452/tabu_seed1.json')
    p.add_argument('--out',required=True);p.add_argument('--seed',type=int,default=1)
    p.add_argument('--steps',type=int,default=2000);p.add_argument('--batches',type=int,default=5)
    p.add_argument('--tenure',type=int,default=15);a=p.parse_args()
    st,rec=structure('s7');y=np.array(json.load(open('data/h3452_stage1.json'))['y'])
    verify(y,rec,[[0,0,0,1]])
    prime=np.array(json.load(open('logs/h3452/prime_gate.json'))['orbit_values'])
    hint=np.array(json.load(open(a.input))['orbit_signs'],np.int64)
    eps,status=seed(y,prime,hint);K,S=coefficients(st);ix,quad=moves(K)
    # Check 100 arbitrary legal rectangle updates against full quadratic evaluation.
    fields,res=rebuild(eps,K,S,st.cls_target);checked=0
    for n,ids in enumerate(ix):
        aa,bb,cc,dd=ids;sg=eps[0,aa]
        if eps[0,bb]!=-sg or eps[0,cc]!=-sg or eps[0,dd]!=sg:continue
        changed=eps.copy();changed[0,ids]*=-1
        delta=2*sg*(-fields[0,:,aa]+fields[0,:,bb]+fields[0,:,cc]-fields[0,:,dd])+quad[n]
        direct=sum(st.paf_classes(st.expand(i,e)) for i,e in enumerate(changed))-st.cls_target
        assert np.array_equal(res+delta,direct)
        checked+=1
        if checked==100:break
    assert checked==100
    runs=[]
    for batch in range(a.batches):
        t=time.monotonic();cost,eps,neval=search(eps,K,S,st.cls_target,ix,quad,a.seed+batch,a.steps,a.tenure)
        xs=np.array([st.expand(i,e) for i,e in enumerate(eps)]);xx=xs.reshape(4,20,43)
        assert np.array_equal(xx.sum(2),y)
        compressed=xx.sum(1)
        assert np.array_equal(eps.reshape(4,20,7).sum(1),2*prime)
        paf=np.array([sum(int(row@np.roll(row,-g)) for row in compressed) for g in range(43)])
        assert np.array_equal(paf,np.array([3452]+[0]*42))
        full=sum(st.paf(row) for row in xs)-st.F
        assert np.sum(full[st.cls_rep]**2)==cost
        runs.append(dict(seed=a.seed+batch,steps=a.steps,tenure=a.tenure,evaluated=neval,seconds=time.monotonic()-t,cost=cost/256))
        out=dict(order=3452,stage='full core lift with two fixed projections',input=a.input,
                 status='CORE_FOUND' if cost==0 else 'NO_CORE_FOUND',orbit_signs=eps.tolist(),
                 class_residual=full[st.cls_rep].tolist(),cost=cost/256,
                 max_abs_residual=int(abs(full).max()),nonzero_lags=int(np.count_nonzero(full)),
                 original_compression_preserved=True,prime_compression_preserved=True,
                 prime_paf_exact=True,rectangle_tests=checked,runs=runs)
        dest=Path(a.out);dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps(runs[-1]),flush=True)
        if cost==0:break


if __name__=='__main__':main()
