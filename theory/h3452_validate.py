"""Validate the stage-one certificate and compatible border; does not claim H3452."""
import json
from pathlib import Path
import numpy as np
from h3452_compress import structure, verify, solve
import bt_border as B


def main():
    c=json.load(open('data/h3452_stage1.json'));st,rec=structure('s7')
    assert verify(c['y'],rec,c['Ngens'])
    y=np.array(c['y'],np.int64);z=y.reshape(4,4,5).sum(2)
    assert np.array_equal(z,c['z'])
    b=json.load(open('logs/h3452/border.json'))
    assert np.array_equal(b['z'],z)
    X,R,C=(np.array(b[k],np.int64) for k in ('X','Rh','Ch'))
    Q=B.QData(st);Mc=B.gs_array(z,Q.add,Q.neg)
    assert X.shape==(12,12) and R.shape==(12,16) and C.shape==(16,12)
    assert all(np.isin(v,[-1,1]).all() for v in (X,R,C))
    assert np.array_equal(X@X.T,12*np.eye(12,dtype=int))
    assert np.array_equal(R@R.T,16*np.eye(12,dtype=int))
    assert np.array_equal(C@C.T,np.kron(np.eye(4,dtype=int),Q.Gmat))
    assert not (X@C.T+R@Mc.T).any()
    bad=y.copy();bad[0,0]+=2
    rejected=False
    try:verify(bad,rec,c['Ngens'])
    except AssertionError:rejected=True
    assert rejected
    # Recheck the saved lift's full residual, independently of incremental updates.
    lift=json.load(open('logs/h3452/lift_seed1.json'))
    xs=np.array([st.expand(i,eps) for i,eps in enumerate(lift['orbit_signs'])])
    assert np.isin(xs,[-1,1]).all()
    assert np.array_equal(xs.reshape(4,20,43).sum(2),y)
    residual=sum(st.paf(x) for x in xs)-st.F
    assert np.array_equal(residual[st.cls_rep],lift['class_residual'])
    assert residual.any(), 'update report if a core has been found'
    out=dict(stage1_all_lags_exact=True,border_identities_exact=True,
             corrupted_certificate_rejected=True,lift_residual_reproduced=True,
             full_core_found=False,hadamard_matrix_found=False)
    Path('logs/h3452/validation.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out))


if __name__=='__main__':main()
