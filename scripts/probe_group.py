import numpy as np, gzip, itertools
def load(N): return np.loadtxt(gzip.open(f'data/hadamard_{N}.csv.gz'),delimiter=',',dtype=np.int8).astype(np.int64)
def tuples(m, maxlen=4):
    out=[]
    def rec(rem, cur):
        if rem==1:
            if 2<=len(cur)<=maxlen: out.append(tuple(cur))
            return
        for d in range(2,rem+1):
            if rem%d==0 and len(cur)<maxlen: rec(rem//d, cur+[d])
    rec(m,[]); return out
def shift_index(mods, e, sign=1):
    m=int(np.prod(mods)); idx=np.arange(m)
    coords=np.array(np.unravel_index(idx,mods)).T
    coords=(coords+sign*np.array(e))%np.array(mods)
    return np.ravel_multi_index(coords.T,mods)
def rank1(T): return np.mean(T*T[0:1,:]*T[:,0:1]*T[0,0]==1)
for N,b in [(1388,5),(1436,7),(1916,3)]:
    H=load(N); v=N//4; m=v-b; off=4*b
    M=H[off:off+m,off:off+m]; M2=H[off:off+m,off+m:off+2*m]
    res=[]
    for mods in tuples(m):
        sc=[];sc2=[]
        for t in range(len(mods)):
            e=[0]*len(mods); e[t]=1
            p=shift_index(mods,e); q=shift_index(mods,e,-1)
            sc.append(rank1(M*M[p][:,p]))           # circ: M[i+e,j+e]
            sc2.append(rank1(M2*M2[p][:,q]))        # back: M2[i+e,j-e]
        res.append((min(sc),min(sc2),mods))
    res.sort(key=lambda t:-t[0]); print(N,m,res[:5],flush=True)
