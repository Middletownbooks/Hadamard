"""Heuristic search of fixed-Q Z43 compressions; zero residuals verified exactly.

Uses small per-cell pools, not exponentially large full-block lists. Restriction:
S7 or S14, with |compressed value| <= --max-abs.
Failure proves nothing. Even success is only a compressed solution.
"""
import argparse, itertools, json, time
from pathlib import Path
import numpy as np
from numba import njit
from h3452_compress import structure, domains, verify


def pools(qt, z, max_abs, kind='s14'):
    ds = domains(qt)
    allp, offsets, lengths = [], [], []
    for block in range(4):
        for cells in (([0], [1,2], [3]) if kind == 's14' else ([0], [1], [2], [3])):
            pos = np.concatenate([np.arange(5*c,5*c+5) for c in cells])
            ids = np.unique(qt.borb[0][pos])
            vals = np.array(list(itertools.product(*[[v for v in ds[u] if abs(v)<=max_abs] for u in ids])), dtype=np.int64)
            sub = vals[:,np.searchsorted(ids,qt.borb[0][pos])]
            ok = (sub.reshape(-1,len(cells),5).sum(2)==np.array(z[block])[cells]).all(1)
            # Other positions in the four blocks have odd values, so at least one unit of energy each.
            ok &= (sub*sub).sum(1) <= int(qt.FN[0])-(80-len(pos))
            sub=sub[ok]
            a=np.zeros((len(sub),20),np.int64);a[:,pos]=sub
            assert len(a)>0
            offsets.append(sum(lengths));lengths.append(len(a));allp.append(a)
    return np.concatenate(allp),np.array(offsets),np.array(lengths)


@njit
def paf(row,add):
    p=np.zeros(20,np.int64)
    for w in range(20):
        for u in range(20):p[w]+=row[u]*row[add[u,w]]
    return p


@njit
def run(pool,offsets,lengths,add,target,seed,restarts,steps,neighbors,use_swaps):
    np.random.seed(seed)
    best=1e100;besty=np.zeros((4,20),np.int64)
    attempts=0
    groups=len(lengths);per_block=groups//4
    for restart in range(restarts):
        choices=np.zeros(groups,np.int64);y=np.zeros((4,20),np.int64)
        for g in range(groups):
            choices[g]=offsets[g]+np.random.randint(lengths[g])
            y[g//per_block]+=pool[choices[g]]
        ps=np.zeros((4,20),np.int64)
        for i in range(4):ps[i]=paf(y[i],add)
        res=ps.sum(0)-target
        cost=(res*res).sum()/784.0
        for step in range(steps):
            g=np.random.randint(groups);i=g//per_block
            nxt=offsets[g]+np.random.randint(lengths[g])
            if use_swaps and per_block==4 and np.random.random()<0.8:
                nxt=neighbors[choices[g],np.random.randint(10)]
            row=y[i]-pool[choices[g]]+pool[nxt]
            pp=paf(row,add)
            rr=res+pp-ps[i]
            cc=(rr*rr).sum()/784.0
            temp=0.15+150.0*(1-step/steps)**3
            if cc<=cost or np.random.random()<np.exp((cost-cc)/temp):
                y[i]=row;ps[i]=pp;res=rr;cost=cc;choices[g]=nxt
            attempts+=1
            if cost<best:
                best=cost;besty=y.copy()
                if best==0:return best,besty,attempts
    return best,besty,attempts


def main():
    p=argparse.ArgumentParser();p.add_argument('--z-index',type=int,default=1)
    p.add_argument('--kind',choices=['s14','s7','cyclic'],default='s14')
    p.add_argument('--proposal',choices=['uniform','swap'],default='swap')
    p.add_argument('--max-abs',type=int,default=15);p.add_argument('--seed',type=int,default=1)
    p.add_argument('--restarts',type=int,default=50);p.add_argument('--steps',type=int,default=100000)
    p.add_argument('--out',required=True);a=p.parse_args()
    st,rec=structure(a.kind);ng=[[0]*(st.G.k-1)+[1]];qt=st.quotient(st.G.subgroup(ng))
    entry=next(json.loads(l) for l in open('logs/bt/border_census.jsonl') if json.loads(l)['key']=='fbf4f510b2566552')
    z=entry['detail'][a.z_index]['z']
    pool,off,lens=pools(qt,z,a.max_abs,a.kind)
    neighbors=np.repeat(np.arange(len(pool))[:,None],10,axis=1)
    if a.kind!='s14':
        for g,(start,length) in enumerate(zip(off,lens)):
            lookup={tuple(pool[j]):j for j in range(start,start+length)}
            for j in range(start,start+length):
                for k,(u,v) in enumerate(itertools.combinations(range(5*(g%4),5*(g%4)+5),2)):
                    row=pool[j].copy();row[u],row[v]=row[v],row[u]
                    neighbors[j,k]=lookup[tuple(row)]
    t=time.monotonic();cost,y,attempts=run(pool,off,lens,qt.add,qt.FN,a.seed,a.restarts,a.steps,neighbors,a.proposal=='swap')
    residual=sum(qt.paf_classes(row) for row in y)-qt.cls_target
    out=dict(stage='Z43 compression only',order=3452,kind=a.kind,z=z,seed=a.seed,
             max_abs=a.max_abs,restarts=a.restarts,steps=a.steps,attempts=attempts,
             proposal=a.proposal,
             seconds=time.monotonic()-t,cost=cost,pool_sizes=lens.tolist(),
             y=y.tolist(),class_residual=residual.tolist(),spec=rec,Ngens=ng)
    assert np.array_equal(y.reshape(4,4,5).sum(2),z)
    if cost==0:out['independent_exact_check']=verify(y,rec,ng)
    out['status']='FEASIBLE' if cost==0 else 'NO_SOLUTION_FOUND'
    dest=Path(a.out);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('y','spec')}),flush=True)


if __name__=='__main__':main()
