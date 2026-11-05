"""Recheck saved lift candidates using explicit group translations, no search updates."""
import json
from pathlib import Path
import numpy as np


def main():
    y=np.array(json.load(open('data/h3452_stage1.json'))['y'],np.int64)
    epsilon=np.where((y-1)%14==0,1,-1)
    cls=[];labels=np.zeros(43,np.int64)
    for p in range(1,43):
        if labels[p]:continue
        c=sorted({p*pow(4,j,43)%43 for j in range(7)});cls.append(c);labels[c]=len(cls)
    # The ordering is independently reconstructed from the smallest representatives.
    prime=json.load(open('data/h3452_prime_projection.json'))
    v=np.array(prime['orbit_values'],np.int64);px=v[:,labels]
    pp=np.array([sum(int(row@np.roll(row,-g)) for row in px) for g in range(43)])
    assert np.array_equal(pp,np.array([863]+[0]*42))
    assert np.array_equal(v[:,0]*2,epsilon.sum(1))
    assert np.array_equal(px.sum(1)*2,y.sum(1))
    target=np.full((4,5,43),4,np.int64);target[0]=-12;target[0,0,0]=3440
    # Class representatives for <multiplier 4, inversion>, in first-occurrence order.
    representatives=[];seen=set()
    for q in range(4):
        for u in range(5):
            for p in range(43):
                if (q,u,p) in seen:continue
                representatives.append(q*215+u*43+p)
                for j in range(7):
                    seen.add((q,u,p*pow(4,j,43)%43))
                    seen.add((q,(-u)%5,-p*pow(4,j,43)%43))
    files=['exact_full_seed1','repair_cells10_seed1','repair_cells4_seed3','tabu_seed1',
           'tabu_seed10','tabu_seed100','tabu_seed200','double_fibre_seed1','double_fibre_seed20']
    results=[]
    for name in files:
        d=json.load(open(f'logs/h3452/{name}.json'))
        signs=np.array(d['orbit_signs'],np.int64).reshape(4,20,7)
        assert np.isin(signs,[-1,1]).all()
        xs=signs[:,:,labels].reshape(4,4,5,43)
        assert np.array_equal(xs.sum(3).reshape(4,20),y)
        assert np.array_equal(xs[:,:,:,0].reshape(4,20),epsilon)
        paf=np.zeros((4,5,43),np.int64)
        for q in range(4):
            for u in range(5):
                for p in range(43):
                    shifted=np.roll(np.roll(xs[:,np.arange(4)^q],-u,axis=2),-p,axis=3)
                    paf[q,u,p]=int(np.sum(xs*shifted))
        res=paf-target;cr=res.ravel()[representatives]
        assert np.array_equal(cr,d['class_residual'])
        assert float(np.sum(cr**2))/256==d['cost']
        assert int(abs(res).max())==d['max_abs_residual']
        assert int(np.count_nonzero(res))==d['nonzero_lags']
        if name.startswith('double_fibre'):
            assert np.array_equal(signs.sum(1),2*v)
            hcomp=xs.sum((1,2));check=np.array([sum(int(row@np.roll(row,-p)) for row in hcomp) for p in range(43)])
            assert np.array_equal(check,np.array([3452]+[0]*42))
        results.append(dict(file=name,cost=d['cost'],nonzero_lags=int(np.count_nonzero(res)),
                            nonzero_classes=int(np.count_nonzero(cr)),max_abs_residual=int(abs(res).max())))
    out=dict(independent_group_translation_check=True,all_candidates_preserve_original_compression=True,
             prime_projection_exact=True,double_projection_candidates_checked=True,results=results,
             core_found=any(r['cost']==0 for r in results))
    Path('logs/h3452/lift_validation.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out))


if __name__=='__main__':main()
