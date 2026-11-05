"""Necessary Z43 quotient gate for the fixed Z43-fibre problem.

Compress along H=Z2^2 x Z5, orthogonally to the saved compression. All values
are even; after division by 2 the four integer sequences must have PAF sum
863 delta on Z43. Their values at zero are fixed by the saved certificate.
FEASIBLE is only a quotient solution; UNKNOWN is inconclusive.
"""
import argparse,json,time
import numpy as np
from ortools.sat.python import cp_model


def main():
    p=argparse.ArgumentParser();p.add_argument('--secs',type=float,default=60)
    p.add_argument('--out',required=True);a=p.parse_args()
    y=np.array(json.load(open('data/h3452_stage1.json'))['y'],np.int64)
    epsilon=np.where((y-1)%14==0,1,-1)
    zero=epsilon.sum(1)//2;rows=y.sum(1)//2
    orbit=np.zeros(43,np.int64);classes=[]
    for x in range(1,43):
        if orbit[x]:continue
        c=sorted({x*pow(4,j,43)%43 for j in range(7)})
        classes.append(c);orbit[c]=len(classes)
    K=np.zeros((43,7,7),np.int64)
    for lag in range(43):
        for u in range(43):K[lag,orbit[u],orbit[(u+lag)%43]]+=1
    m=cp_model.CpModel();X=[];sq=[]
    for i in range(4):
        row=[m.new_constant(int(zero[i]))]+[m.new_int_var(-10,10,f'x{i}_{j}') for j in range(6)]
        m.add(int(zero[i])+7*sum(row[1:])==int(rows[i]));X.append(row)
    products={}
    for lag in [0,1,3,7]:
        terms=[]
        for i in range(4):
            for u in range(7):
                for v in range(u,7):
                    coeff=int(K[lag,u,v]+(K[lag,v,u] if v!=u else 0))
                    if not coeff:continue
                    key=i,u,v
                    if key not in products:
                        pp=m.new_int_var(0 if u==v else -100,100,f'p{i}_{u}_{v}')
                        m.add_multiplication_equality(pp,[X[i][u],X[i][v]]);products[key]=pp
                    terms.append(coeff*products[key])
        m.add(sum(terms)==(863 if lag==0 else 0))
    solver=cp_model.CpSolver();solver.parameters.max_time_in_seconds=a.secs;solver.parameters.num_workers=4
    status=solver.solve(m)
    out=dict(status=solver.status_name(status),seconds=solver.wall_time,scope='necessary H-compression onto Z43 only',
             fixed_zero=zero.tolist(),row_sums=rows.tolist(),cyclotomic_classes=classes)
    if status in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        xx=np.array([[solver.value(v) for v in row] for row in X],np.int64)
        expanded=xx[:,orbit]
        paf=np.array([sum(int(row@np.roll(row,-g)) for row in expanded) for g in range(43)])
        assert np.array_equal(paf,np.array([863]+[0]*42))
        assert np.array_equal(expanded.sum(1),rows)
        out.update(orbit_values=xx.tolist(),paf=paf.tolist(),all_lags_exact=True)
    open(a.out,'w').write(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':main()
