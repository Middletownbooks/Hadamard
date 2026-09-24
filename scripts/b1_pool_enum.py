"""Complete enumeration of U-invariant +-1 blocks with row sum a and PSD <= bound at every frequency class.
Meet-in-the-middle over orbit bits: F(w) = F_lo(w) + F_hi(w), row sum = rs_lo + rs_hi.
usage (lib): enum_pool(R, a, bound) -> sorted uint64 array of orbit bitmasks (bit q = orbit q is +1)."""
import numpy as np, numba as nb

def _tables(R,lo,hi):
    nb_=hi-lo; n=1<<nb_
    bits=((np.arange(n)[:,None]>>np.arange(nb_)[None,:])&1)*2-1                  # (n, nb_) signs
    rs=bits@R.size[lo:hi]
    F=bits@R.S[1:,lo:hi].T                                                       # skip omega=0
    return rs.astype(np.int64), np.ascontiguousarray(F.real), np.ascontiguousarray(F.imag)

@nb.njit(parallel=True,cache=True)
def _scan(rsH,FrH,FiH,rsLs,oL,FrL,FiL,starts,ends,a,bound,shift,write,cnt,out,off):
    NH=rsH.shape[0]; F=FrH.shape[1]
    for h in nb.prange(NH):
        t=a-rsH[h]
        # binary search for t in sorted rsLs
        lo=np.searchsorted(rsLs,t); hi=np.searchsorted(rsLs,t+1)
        c=0; n=off[h] if write else 0
        for p in range(lo,hi):
            l=oL[p]; ok=True
            for f in range(F):
                xr=FrH[h,f]+FrL[l,f]; xi=FiH[h,f]+FiL[l,f]
                if xr*xr+xi*xi>bound: ok=False; break
            if ok:
                if write: out[n]=(np.uint64(h)<<np.uint64(shift))|np.uint64(l); n+=1
                else: c+=1
        if not write: cnt[h]=c
    return 0

def enum_pool(R,a,bound,hkeep=None):
    """hkeep: optional sorted array of high-half indices to scan (random subsample of the full pool)."""
    h=R.r//2
    rsL,FrL,FiL=_tables(R,0,h); rsH,FrH,FiH=_tables(R,h,R.r)
    if hkeep is not None: rsH,FrH,FiH=rsH[hkeep],FrH[hkeep],FiH[hkeep]
    oL=np.argsort(rsL,kind='stable'); rsLs=rsL[oL]
    cnt=np.zeros(len(rsH),np.int64); dummy=np.zeros(1,np.uint64); off=np.zeros(len(rsH),np.int64)
    _scan(rsH,FrH,FiH,rsLs,oL,FrL,FiL,0,0,a,float(bound),h,False,cnt,dummy,off)
    off=np.concatenate([[0],np.cumsum(cnt)[:-1]]); out=np.empty(int(cnt.sum()),np.uint64)
    _scan(rsH,FrH,FiH,rsLs,oL,FrL,FiL,0,0,a,float(bound),h,True,cnt,out,off)
    if hkeep is not None:
        hi=np.asarray(hkeep,np.uint64)[(out>>np.uint64(h)).astype(np.int64)]
        out=(hi<<np.uint64(h))|(out&np.uint64((1<<h)-1))
    out.sort(); return out

def bits_to_Y(bits,r):
    return (((bits[:,None]>>np.arange(r,dtype=np.uint64)[None,:])&np.uint64(1)).astype(np.int64))*2-1

@nb.njit(parallel=True,cache=True)
def _scan_hist(rsH,FrH,FiH,rsLs,oL,FrL,FiL,a,bound,fsel,H,cnt):
    """like _scan, but for every block passing the full PSD bound record round(P) at frequencies fsel in H."""
    NH=rsH.shape[0]; F=FrH.shape[1]; nfs=fsel.shape[0]
    for h in nb.prange(NH):
        t=a-rsH[h]; lo=np.searchsorted(rsLs,t); hi=np.searchsorted(rsLs,t+1); c=0
        for p in range(lo,hi):
            l=oL[p]; ok=True
            for f in range(F):
                xr=FrH[h,f]+FrL[l,f]; xi=FiH[h,f]+FiL[l,f]
                if xr*xr+xi*xi>bound: ok=False; break
            if ok:
                c+=1
                for j in range(nfs):
                    f=fsel[j]; xr=FrH[h,f]+FrL[l,f]; xi=FiH[h,f]+FiL[l,f]
                    H[j,int(xr*xr+xi*xi+0.5)]=True
        cnt[h]=c
    return 0

def pool_value_sets(R,a,bound,fsel):
    """exact sets of achievable round(P(w)) at frequency-class indices fsel (into R.freqs[1:]) over the complete pool."""
    h=R.r//2
    rsL,FrL,FiL=_tables(R,0,h); rsH,FrH,FiH=_tables(R,h,R.r)
    oL=np.argsort(rsL,kind='stable'); rsLs=rsL[oL]
    Hh=np.zeros((len(fsel),int(bound)+2),np.bool_); cnt=np.zeros(len(rsH),np.int64)
    _scan_hist(rsH,FrH,FiH,rsLs,oL,FrL,FiL,a,float(bound),np.asarray(fsel,np.int64),Hh,cnt)
    return int(cnt.sum()), [np.where(Hh[j])[0] for j in range(len(fsel))]
