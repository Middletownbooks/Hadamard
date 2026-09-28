"""Three-sided exact lift (fast replacement for the per-class exact join of the two-sided route).
For a yc class: per block, all lifts matching xc (z mod t) and yc (z mod p2), PSD-filtered.  Each lift is further
compressed along the Z3 coordinate: wc(a,z) = sum_b f(a,b,z) on Z2 x Z_h.  The wc-level PAF equations
sum_k PAF(wc_k)(s) = T_w(s) := sum_{beta} T_full(alpha,beta,z) are implied by the full equations, so
  (1) exact 4-list join over the DISTINCT wc values of each block (small lists) -> all valid wc 4-tuples;
  (2) for each valid wc 4-tuple, exact full join over the lifts with those wc values (small sub-lists).
Every hit is verified on the full group (sum_k PAF_k = T_full).  Exact: no sampling, no bucketing heuristics.
lib: wc_lift(g, comps, foot, xc, ycs, u) -> (solutions [F arrays], stats)"""
import sys, time, numpy as np
sys.path.insert(0,'scripts')
from bge3_twosided import enum_lifts_yc

def _lag_reps(cells_shape,u):
    """orbit representatives of <u,-1> on the finite abelian group Z_n1 x Z_n2 (mixed radix), excluding nothing."""
    n1,n2=cells_shape; seen=set(); reps=[]
    for x in range(n1*n2):
        if x in seen: continue
        a,z=divmod(x,n2); orb=set(); y=(a,z)
        while True:
            for w in (y,((-y[0])%n1,(-y[1])%n2)): orb.add(w[0]*n2+w[1])
            y=((y[0]*u)%n1,(y[1]*u)%n2)
            if y==(a,z): break
        seen|=orb; reps.append(x)
    return np.array(reps)

def _paf2(Wv,n1,n2,reps):
    F=Wv.reshape(-1,n1,n2).astype(float)
    return np.rint(np.fft.ifft2(np.abs(np.fft.fft2(F))**2).real).reshape(len(Wv),-1)[:,reps].astype(np.int64)

def _pairs_keys(P1,P2,W):
    """all pairs (i,j) of rows of P1 x P2 with 64-bit linear hash of P1[i]+P2[j]; returns keys, I, J"""
    h1=(P1.astype(np.uint64)*W).sum(1); h2=(P2.astype(np.uint64)*W).sum(1)
    K=(h1[:,None]+h2[None,:]).ravel(); I=np.repeat(np.arange(len(P1)),len(P2)); J=np.tile(np.arange(len(P2)),len(P1))
    return K,I,J

def wc_lift(g,comps,foot,xc,ycs,u,log=print):
    t0=time.time(); co=g['co']; v=g['v']; lab=g['lab']; m=g['m']; ADD=g['ADD']; T_full=g['T_full']
    h=int(co[:,2].max())+1; mods=np.array([2,3,h]); chars=co[co[:,2]!=0]
    ph=np.exp(2j*np.pi*((co[:,None,:]*chars[None,:,:])/mods).sum(-1))
    wcell=(co[:,0]%2)*h+co[:,2]; reps=_lag_reps((2,h),u)
    Tw=np.zeros(2*h,np.int64); np.add.at(Tw,wcell,T_full)
    Fs=[]; Wd=[]; inv=[]
    for k in range(4):
        Y=enum_lifts_yc(g,comps,xc[k],ycs[k],foot); keep=[]
        for a in range(0,len(Y),50000):
            F=Y[a:a+50000][:,lab].astype(np.float32); keep.append(((np.abs(F@ph)**2)<=4*v+1e-3).all(1))
        Y=Y[np.concatenate(keep)] if keep else Y; F=Y[:,lab].astype(np.int8)
        W=np.zeros((len(F),2*h),np.int16); np.add.at(W.T,wcell,F.T)
        u_,iv=np.unique(W,axis=0,return_inverse=True); Fs.append(F); Wd.append(u_); inv.append(iv.ravel())
    stats=dict(psd_lifts=[int(len(f)) for f in Fs],distinct_wc=[int(len(w)) for w in Wd])
    P=[_paf2(Wd[k],2,h,reps) for k in range(4)]; P[0]=P[0]-Tw[reps]
    rng=np.random.default_rng(11); Wt=rng.integers(1,2**62,size=len(reps),dtype=np.int64).astype(np.uint64)
    KA,IA,JA=_pairs_keys(P[0],P[1],Wt); KC,IC,JC=_pairs_keys(-P[2],-P[3],Wt)
    o=np.argsort(KC); KC=KC[o]; IC=IC[o]; JC=JC[o]
    lo=np.searchsorted(KC,KA); hi=np.searchsorted(KC,KA,side='right'); hit=np.nonzero(hi>lo)[0]
    tuples=[]
    for q in hit:
        for r in range(lo[q],hi[q]):
            a,b,c,d=int(IA[q]),int(JA[q]),int(IC[r]),int(JC[r])
            if np.array_equal(P[0][a]+P[1][b],-(P[2][c]+P[3][d])): tuples.append((a,b,c,d))
    stats['wc_pairs']=[int(len(KA)),int(len(KC))]; stats['wc_tuples']=len(tuples)
    # full exact join inside each valid wc tuple.  PAF is constant on the lag classes of <u,-1>, so matching at the
    # class representatives is equivalent to matching at all lags; each lift gets a 64-bit linear hash once.
    from bge3_residue import class_reps
    creps=class_reps(g,u,h); Tr=T_full[creps]
    Wf=np.random.default_rng(13).integers(1,2**62,size=len(creps),dtype=np.int64).astype(np.uint64)
    PR=[]; HK=[]
    for k in range(4):
        R=np.empty((len(Fs[k]),len(creps)),np.int64)
        for a0 in range(0,len(Fs[k]),20000):
            x=Fs[k][a0:a0+20000].reshape(-1,2,3,h).astype(float)
            R[a0:a0+20000]=np.rint(np.fft.ifftn(np.abs(np.fft.fftn(x,axes=(1,2,3)))**2,axes=(1,2,3)).real).reshape(len(x),-1)[:,creps]
        PR.append(R); HK.append((R.astype(np.uint64)*Wf).sum(1))
    hT=np.uint64((Tr.astype(np.uint64)*Wf).sum())
    idx=[]
    for k in range(4):
        order=np.argsort(inv[k],kind='stable'); srt=inv[k][order]; idx.append((order,np.searchsorted(srt,np.arange(len(Wd[k])+1))))
    def members(k,w): order,b=idx[k]; return order[b[w]:b[w+1]]
    sols=[]; npairs=0; ncand=0
    for (a,b,c,d) in tuples:
        ia=members(0,a); ib=members(1,b); ic=members(2,c); id_=members(3,d)
        ka=(HK[0][ia][:,None]+HK[1][ib][None,:]-hT).ravel(); kc=(np.uint64(0)-(HK[2][ic][:,None]+HK[3][id_][None,:])).ravel()
        npairs+=len(ka)+len(kc)
        oc=np.argsort(kc); kcs=kc[oc]; l2=np.searchsorted(kcs,ka); h2=np.searchsorted(kcs,ka,side='right')
        for q in np.nonzero(h2>l2)[0]:
            for r in range(l2[q],h2[q]):
                j2=oc[r]; x0,x1=divmod(int(q),len(ib)); x2,x3=divmod(int(j2),len(id_))
                sel=(ia[x0],ib[x1],ic[x2],id_[x3]); ncand+=1
                if not np.array_equal(sum(PR[k][sel[k]] for k in range(4)),Tr): continue
                Fq=np.array([Fs[k][sel[k]] for k in range(4)]).astype(np.int64)
                Pq=sum(np.array([int(np.dot(Fq[k],Fq[k][ADD[:,gg]])) for gg in range(m)]) for k in range(4))
                if np.array_equal(Pq[1:],T_full[1:]): sols.append(Fq)
    stats['candidates']=ncand
    stats['subjoin_pairs']=int(npairs); stats['solutions']=len(sols); stats['secs']=round(time.time()-t0,1)
    return sols,stats
