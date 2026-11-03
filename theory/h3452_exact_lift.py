"""Exact fixed-compression lift model and bounded large-neighborhood repairs.

Spin products are linearized with Boolean ANDs. All PAF equations are integer
equations; cell cardinalities preserve the input Z43 compression exactly.
UNKNOWN is inconclusive. INFEASIBLE is scoped to the selected free cells.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from ortools.sat.python import cp_model
from h3452_compress import structure,verify
from h3452_lift import coefficients,evaluate


def residual(st,eps):
    return sum(st.paf_classes(st.expand(i,e)) for i,e in enumerate(eps))-st.cls_target


def model_for(st,K,S,eps,cells,objective):
    m=cp_model.CpModel();eps=np.asarray(eps,np.int64);base=eps.copy()
    free=[[] for _ in range(4)]
    for i,c in cells:free[i].extend(range(7*c+1,7*c+7))
    for i in range(4):base[i,free[i]]=-1
    const=sum(evaluate(e,K) for e in base)-st.cls_target
    assert not (const%4).any()
    expressions=[int(v//4) for v in const]
    bits={};aux=[]
    for i in range(4):
        field=np.einsum('jab,b->ja',S,base[i])
        for a in free[i]:
            b=m.new_bool_var(f'b{i}_{a}');bits[i,a]=b
            m.add_hint(b,int((eps[i,a]+1)//2))
            lin=2*field[:,a]+4*K[:,a,a]
            assert not (lin%4).any()
            for j in np.flatnonzero(lin):expressions[j]+=int(lin[j]//4)*b
        for ia,a in enumerate(free[i]):
            for b in free[i][ia+1:]:
                coeff=S[:,a,b]
                if not coeff.any():continue
                p=m.new_bool_var(f'p{i}_{a}_{b}')
                m.add_multiplication_equality(p,[bits[i,a],bits[i,b]])
                m.add_hint(p,int(eps[i,a]==1 and eps[i,b]==1))
                aux.append((i,a,b,p))
                for j in np.flatnonzero(coeff):expressions[j]+=int(coeff[j])*p
    for i,c in cells:
        m.add(sum(bits[i,a] for a in range(7*c+1,7*c+7))==int(np.count_nonzero(eps[i,7*c+1:7*c+7]==1)))
    squares=[];old=residual(st,eps)
    assert not (old%4).any()
    budget=int(np.sum((old//4)**2))
    if objective:
        bound=int(np.sqrt(budget))
        for j,expr in enumerate(expressions):
            r=m.new_int_var(-bound,bound,f'r{j}');m.add(r==expr)
            ss=m.new_int_var(0,budget,f'sq{j}');m.add_multiplication_equality(ss,[r,r])
            m.add_hint(r,int(old[j]//4));m.add_hint(ss,int((old[j]//4)**2));squares.append(ss)
        m.add(sum(squares)<=budget);m.minimize(sum(squares))
    else:
        for expr in expressions:m.add(expr==0)
    # Verify the Boolean polynomial at several independently evaluated assignments.
    # Each expression is checked using the model's integer linear-expression proto.
    rng=np.random.default_rng(3452)
    for trial in range(4):
        test=eps.copy()
        for i,c in cells:
            test[i,7*c+1:7*c+7]=rng.permutation(test[i,7*c+1:7*c+7])
        values={v.index:int((test[i,a]+1)//2) for (i,a),v in bits.items()}
        values.update({p.index:int(test[i,a]==1 and test[i,b]==1) for i,a,b,p in aux})
        predicted=[]
        for expr in expressions:
            if isinstance(expr,(int,np.integer)):predicted.append(int(expr));continue
            flat=cp_model.FlatIntExpr(expr)
            predicted.append(flat.offset+sum(int(c)*values[v.index] for c,v in zip(flat.coeffs,flat.vars)))
        assert np.array_equal(4*np.array(predicted),residual(st,test))
    return m,bits,budget


def run(st,K,S,eps,cells,secs,workers,seed,objective=True):
    t=time.monotonic();m,bits,budget=model_for(st,K,S,eps,cells,objective)
    s=cp_model.CpSolver();s.parameters.max_time_in_seconds=secs
    s.parameters.num_workers=workers;s.parameters.random_seed=seed
    status=s.solve(m)
    out=dict(status=s.status_name(status),solver_seconds=s.wall_time,total_seconds=time.monotonic()-t,
             free_cells=[list(c) for c in cells],objective=objective,seed=seed,workers=workers,
             initial_squared_residual=budget*16,branches=s.num_branches,conflicts=s.num_conflicts,
             polynomial_checks=4)
    result=eps.copy()
    if status in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        for (i,a),b in bits.items():result[i,a]=2*s.value(b)-1
        rr=residual(st,result);out['squared_residual']=int(rr@rr)
        if objective:assert int(round(s.objective_value))*16==int(rr@rr)
        else:assert not rr.any()
    return result,out


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',default='logs/h3452/lift_seed1.json')
    p.add_argument('--out',required=True);p.add_argument('--secs',type=float,default=5)
    p.add_argument('--workers',type=int,default=1);p.add_argument('--seed',type=int,default=1)
    p.add_argument('--rounds',type=int,default=20);p.add_argument('--cells',type=int,default=12)
    p.add_argument('--full',action='store_true');a=p.parse_args()
    st,rec=structure('s7');cert=json.load(open('data/h3452_stage1.json'));y=np.array(cert['y'])
    verify(y,rec,[[0,0,0,1]])
    eps=np.array(json.load(open(a.input))['orbit_signs'],np.int64)
    K,S=coefficients(st);rng=np.random.default_rng(a.seed);runs=[]
    initial=residual(st,eps);allcells=[(i,c) for i in range(4) for c in range(20)]
    for it in range(1 if a.full else a.rounds):
        if a.full:cells=allcells
        else:cells=[allcells[k] for k in sorted(rng.choice(80,min(a.cells,80),replace=False))]
        eps,r=run(st,K,S,eps,cells,a.secs,a.workers,a.seed+it,not a.full)
        r['round']=it;runs.append(r)
        xs=np.array([st.expand(i,e) for i,e in enumerate(eps)])
        assert np.array_equal(xs.reshape(4,20,43).sum(2),y)
        full=sum(st.paf(x) for x in xs)-st.F
        out=dict(stage='full core lift',input=a.input,order=3452,orbit_signs=eps.tolist(),
                 class_residual=full[st.cls_rep].tolist(),cost=float(np.sum(full[st.cls_rep]**2))/256,
                 initial_cost=float(initial@initial)/256,max_abs_residual=int(abs(full).max()),
                 nonzero_lags=int(np.count_nonzero(full)),compression_preserved=True,
                 status='CORE_FOUND' if not full.any() else 'NO_CORE_FOUND',runs=runs)
        path=Path(a.out);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps(dict(round=it,status=r['status'],cost=out['cost'],seconds=r['solver_seconds'])),flush=True)
        if not full.any():break


if __name__=='__main__':main()
