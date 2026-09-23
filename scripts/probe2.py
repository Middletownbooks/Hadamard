import numpy as np, gzip
def circ_scaled_score(Bk, back=False):
    m=Bk.shape[0]
    # Q[i,j] = Bk[i,j]*Bk[i+1,j+1] (circulant) or Bk[i,j]*Bk[i+1,j-1] (back-circulant)
    sh = -1 if back else 1
    Q = Bk*np.roll(np.roll(Bk,-1,0),-sh,1)
    # rank-1 test
    T = Q*Q[0:1,:]*Q[:,0:1]*Q[0,0]
    return np.mean(T==1)
for N in [668,716,892,1132,1244,1388,1436,1676,1772,1852,1916,1948,1964]:
    H=np.loadtxt(gzip.open(f'data/hadamard_{N}.csv.gz'),delimiter=',',dtype=np.int8).astype(np.int32)
    v=N//4; out=[]
    for b in range(0,9):
        m=v-b
        for off in sorted({0,b,4*b}):
          if off+m<=N:
            Bk=H[off:off+m, off:off+m]
            out.append((round(circ_scaled_score(Bk),3),round(circ_scaled_score(Bk,True),3),b,off))
    out.sort(reverse=True); print(N,v,out[:3],flush=True)
