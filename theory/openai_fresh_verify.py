#!/usr/bin/env python3
"""Independent GS build and exact integer Gram check for logged b0 controls.

Imports none of the construction/search modules. Also rejects a corrupted
matrix. Input cores may come from any solver; only small known orders are
claimed in this session.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def check(sequences):
    x=np.asarray(sequences,np.int64)
    assert x.ndim==2 and x.shape[0]==4 and np.all(np.abs(x)==1)
    v=x.shape[1];n=4*v
    # Direct correlation check written independently of the search modules.
    for t in range(v):
        value=sum(int(x[b,j])*int(x[b,(j+t)%v]) for b in range(4) for j in range(v))
        assert value==(4*v if t==0 else 0)
    idx=(np.arange(v)[None,:]-np.arange(v)[:,None])%v
    a,b,c,d=[row[idx] for row in x]
    reverse=np.eye(v,dtype=np.int64)[::-1]
    br=b@reverse;cr=c@reverse;dr=d@reverse
    h=np.block([[a,br,cr,dr],[-br,a,d.T@reverse,-c.T@reverse],
                [-cr,-d.T@reverse,a,b.T@reverse],[-dr,c.T@reverse,-b.T@reverse,a]])
    assert np.array_equal(h@h.T,n*np.eye(n,dtype=np.int64))
    digest=hashlib.sha256(h.astype(np.int8).tobytes()).hexdigest()
    h[0,0]*=-1
    assert not np.array_equal(h@h.T,n*np.eye(n,dtype=np.int64))
    return dict(order=n,integer_Gram_passed=True,corrupted_matrix_rejected=True,sha256_int8=digest)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('inputs',nargs='+')
    ap.add_argument('--out',default='logs/openai_fresh/independent_matrices.json')
    args=ap.parse_args()
    results=[]
    for name in args.inputs:
        data=json.loads(Path(name).read_text())
        assert data['kind']=='b0'
        for rec in data.get('hits',[]):
            results.append(dict(source=name,trial=rec['trial'],**check(rec['sequences'])))
        for rec in data.get('trials',[]) if isinstance(data.get('trials'),list) else []:
            if rec.get('sequences'):
                results.append(dict(source=name,trial=rec['trial'],**check(rec['sequences'])))
            other=rec.get('outside_requested_sums')
            if other:
                results.append(dict(source=name,trial=rec['trial'],outside_requested_sums=True,
                                    **check(other['sequences'])))
    Path(args.out).parent.mkdir(parents=True,exist_ok=True)
    Path(args.out).write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(dict(verified=len(results),orders=sorted({r['order'] for r in results}),
                         all_integer_Gram_checks=True)))
