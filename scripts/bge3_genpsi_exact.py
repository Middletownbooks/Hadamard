"""Generalized bordered architecture (report/bge3_answer.md section 6): Psi = 3 independent characters of Q = Z2^3
(b = 3, c = 8, v = 8h + 3). Target: sum_k PAF_k(q,z) = 4v delta - 4 sum_{psi in Psi} psi(q).
Exact search in the slice of blocks invariant under z -> -z: per-block PSD filter |X^(chi)|^2 <= T^(chi), then a
per-row-sum-profile 4-list join by pair hashing. Result [M]: no core for h = 3 or h = 5 (h = 1 is impossible by parity).
usage: bge3_genpsi_exact.py h"""
import sys, itertools, numpy as np, json
h=int(sys.argv[1]); c=8; b=3; v=c*h+b; m=c*h
Qco=np.array(list(itertools.product(range(2),repeat=3)))
# orbits of z -> -z on Z_h within each row
zorb=[]; seen=set()
for z in range(h):
    if z in seen: continue
    o=sorted({z,(-z)%h}); seen|=set(o); zorb.append(o)
lab=np.zeros((8,h),int); r=0
for q in range(8):
    for o in zorb:
        for z in o: lab[q,z]=r
        r+=1
lab=lab.ravel()                         # element (q,z) at index q*h+z
Psi=[np.array(p) for p in [(1,0,0),(0,1,0),(0,0,1)]]
psiQ=np.array([(-1)**((Qco*p).sum(1)%2) for p in Psi])            # (3,8)
chiQ=np.array(list(itertools.product(range(2),repeat=3)))
inPsi=np.array([any((x==p).all() for p in Psi) for x in chiQ])     # over Q-characters in same order as fft index
That=np.full((8,h),4*v,float); That[inPsi,0]=4*b                    # T^ indexed (Q-char flat via fftn over (2,2,2), omega)
That=That.reshape(2,2,2,h)
# target PAF: T(q,z) = 4v delta - 4 sum_psi psi(q)
Tg=-4*psiQ.sum(0)[:,None]*np.ones((1,h)); Tg=Tg.astype(np.int64); Tg[0,0]=4*m
Tg=Tg.reshape(2,2,2,h)
print('r per block',r,'blocks',2**r,flush=True)
keep=[]; Ps=[]; As=[]
CH=1<<14
for s in range(0,2**r,CH):
    idx_=np.arange(s,min(2**r,s+CH),dtype=np.int64)
    Y=(1-2*((idx_[:,None]>>np.arange(r))&1)).astype(np.int8)
    X=Y[:,lab].reshape(-1,2,2,2,h).astype(float)
    F=np.fft.fftn(X,axes=(1,2,3,4)); Pw=np.abs(F)**2
    ok=(Pw<=That[None]+1e-6).reshape(len(Y),-1).all(1)
    if not ok.any(): continue
    Pf=np.rint(np.fft.ifftn(Pw[ok],axes=(1,2,3,4)).real).astype(np.int64).reshape(int(ok.sum()),-1)
    keep.append(Y[ok]); Ps.append(Pf); As.append(X[ok].reshape(int(ok.sum()),-1).sum(1).astype(np.int64))
Y=np.concatenate(keep); P=np.concatenate(Ps); A=np.concatenate(As)
print('PSD-filtered blocks',len(Y),'row sums',sorted(set(A.tolist()))[:20],flush=True)
rng=np.random.default_rng(1); wv=rng.integers(1,2**62,size=P.shape[1],dtype=np.int64)|1
Hh=(P*wv).sum(1); Ht=(Tg.reshape(-1)*wv).sum()
# per-profile 4-list joins (row sums a1..a4, sum a^2 = 4v)
Avals=sorted(set(A.tolist()))
profs=[p for p in itertools.combinations_with_replacement([x for x in Avals if x>=0],4) if sum(x*x for x in p)==4*v]
print('profiles',profs,flush=True)
found=None
for p in profs:
    for perm in sorted(set(itertools.permutations(p))):
        L=[np.where(A==x)[0] for x in perm]
        if any(len(l)==0 for l in L): continue
        S12=(Hh[L[0]][:,None]+Hh[L[1]][None,:]).ravel(); o=np.argsort(S12); S12s=S12[o]
        for i0 in range(len(L[2])):
            need=Ht-(Hh[L[2][i0]]+Hh[L[3]])
            pos=np.searchsorted(S12s,need); pos=np.clip(pos,0,len(S12s)-1); hit=np.where(S12s[pos]==need)[0]
            for j in hit:
                t=o[pos[j]]; i1,i2=L[0][t//len(L[1])],L[1][t%len(L[1])]; i3,i4=L[2][i0],L[3][j]
                if np.array_equal(P[i1]+P[i2]+P[i3]+P[i4],Tg.reshape(-1)):
                    found=[Y[i][lab].tolist() for i in (i1,i2,i3,i4)]; break
            if found: break
        print('profile',perm,'sizes',[len(l) for l in L],'found',found is not None,flush=True)
        if found: break
    if found: break
print(json.dumps(dict(h=h,v=v,found=found is not None,core=found)))
