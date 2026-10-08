"""Exact decision of small high-symmetry slices (r <= ~24) of the b>=3 problem: does ANY core exist that is invariant
under the symmetry group S?  (This is where an algebraic / cyclotomic construction would have to live.)
1. all 2^r S-invariant blocks; per-block filters: row sum |a| in a feasible profile, PSD |X^(chi)|^2 <= T^(chi) for all
   characters (T^ = 4b on H-perp minus 1, else 4v);  blocks normalised to a >= 0 (negation preserves PAF).
2. for every profile (a1,a2,a3,a4): exact 4-list join on the lag-class PAF vector by 64-bit linear hashing of pair sums,
   candidate matches re-verified exactly.
usage: bge3_slice_exact.py v c Qmods rmax [rmin]      (env GMAX/PAIRMAX for group enumeration)"""
import sys, os, json, itertools, numpy as np
v=int(sys.argv[1]); c=int(sys.argv[2]); Qs=sys.argv[3]; rmax=int(sys.argv[4]); rmin=int(sys.argv[5]) if len(sys.argv)>5 else 1
os.environ.setdefault('GMAX','400'); os.environ.setdefault('PAIRMAX','40')
src=open('scripts/bge3_cost.py').read(); cut=src.index("rng=np.random.default_rng(0)\nres=[]")
sys.argv=['bge3_cost.py',str(v),str(c),Qs,str(rmin),str(rmax)]
exec(src[:cut])
nQ=len(Q)
chis=np.array(np.unravel_index(np.arange(m),mods)).T
That=np.where((chis[:,-1]==0)&(chis[:,:nQ]!=0).any(1),4*b,4*v).reshape(mods)
def decide(lab,grp,r):
    reps=lagclasses(grp); L=len(reps); sz=np.bincount(lab)
    profs=profiles(sz)
    tgt=np.array([-4*b if int(x) in Hset else 4 for x in reps],np.int64)
    if not profs: return dict(status='no row-sum profile')
    A=sorted({a for p in profs for a in p})
    keepY=[]; keepP=[]; keepA=[]
    CH=1<<15
    for s in range(0,2**r,CH):
        idx_=np.arange(s,min(2**r,s+CH),dtype=np.int64)
        Y=(1-2*((idx_[:,None]>>np.arange(r))&1)).astype(np.int8)
        Xs=Y[:,lab]; a=Xs.sum(1).astype(np.int64)
        sel=np.isin(a,A)
        if not sel.any(): continue
        Y=Y[sel]; Xs=Xs[sel]; a=a[sel]
        F=np.fft.fftn(Xs.reshape((-1,)+tuple(mods)).astype(np.float64),axes=tuple(range(1,len(mods)+1)))
        Pw=np.abs(F)**2
        ok=(Pw<=That[None]+1e-6).reshape(len(Y),-1).all(1)
        if not ok.any(): continue
        Pf=np.fft.ifftn(Pw[ok],axes=tuple(range(1,len(mods)+1))).real.reshape(int(ok.sum()),-1)
        keepY.append(Y[ok]); keepP.append(np.rint(Pf[:,reps]).astype(np.int64)); keepA.append(a[ok])
    if not keepY: return dict(status='no block passes PSD filter',profiles=len(profs))
    Y=np.concatenate(keepY); P=np.concatenate(keepP); a=np.concatenate(keepA)
    rng_=np.random.default_rng(7); wv=rng_.integers(1,2**62,size=L,dtype=np.int64)|1
    H=(P*wv).sum(1); Ht=(tgt*wv).sum()
    lists={x:np.where(a==x)[0] for x in A}
    sols=[]; sizes={int(x):int(len(lists[x])) for x in A}
    for p in profs:
        if any(len(lists[x])==0 for x in p): continue
        L1,L2,L3,L4=[lists[x] for x in p]
        # pair hashes left
        h12=(H[L1][:,None]+H[L2][None,:]).ravel(); o=np.argsort(h12); h12s=h12[o]
        for i0 in range(0,len(L3),max(1,2**22//max(1,len(L4)))):
            blk=L3[i0:i0+max(1,2**22//max(1,len(L4)))]
            h34=(H[blk][:,None]+H[L4][None,:]).ravel(); need=Ht-h34
            pos=np.searchsorted(h12s,need); pos=np.clip(pos,0,len(h12s)-1); hit=np.where(h12s[pos]==need)[0]
            for hh in hit:
                # enumerate all left pairs with that hash; verify exactly
                lo=np.searchsorted(h12s,need[hh],'left'); hi=np.searchsorted(h12s,need[hh],'right')
                i3=blk[hh//len(L4)]; i4=L4[hh%len(L4)]
                for t in o[lo:hi]:
                    i1=L1[t//len(L2)]; i2=L2[t%len(L2)]
                    if np.array_equal(P[i1]+P[i2]+P[i3]+P[i4],tgt):
                        sols.append([Y[i].tolist() for i in (i1,i2,i3,i4)])
                        if len(sols)>=5: break
                if len(sols)>=5: break
            if len(sols)>=5: break
        if len(sols)>=5: break
    return dict(status='decided',profiles=len(profs),list_sizes=sizes,solutions_found=len(sols),example=(sols[0] if sols else None))
for S,(d,r,lab,grp,abel) in sorted(seen.items(),key=lambda kv:kv[1][1]):
    if r>rmax or r<rmin: continue
    L=len(lagclasses(grp))
    res=decide(lab,grp,r)
    print(json.dumps(dict(v=v,c=c,Q=Q,h=h,d=d,r=r,L=L,abelian=abel,**res),default=int),flush=True)
