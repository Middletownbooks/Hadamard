"""Extract GS circulant cores (b=0, b=1) from normalised full Hadamard matrices.
H = D_r * H0 * D_c with H0 a GS array; recover each block's defining sequence
up to global sign (and the (-1)^k twist for even m), then choose the combination
that satisfies the periodic complementarity equations."""
import numpy as np, gzip, itertools, sys
def load(N):
    return np.loadtxt(gzip.open(f'data/hadamard_{N}.csv.gz'),delimiter=',',dtype=np.int8).astype(np.int64)
def seq_from_block(Bk, back):
    m=Bk.shape[0]; sols=[]
    # circulant: Bk[i,j]=r_i s_j a_{j-i}; back: Bk[i,j]= r_i s_j a_{i+j}
    for r1 in (1,-1):
        a=np.zeros(m,int); a[0]=1
        # r_0=1 => s_j a_j = Bk[0,j] (circ) ; row1: circ: r1 s_j a_{j-1}=Bk[1,j] -> a_j a_{j-1} = r1 Bk[1,j]Bk[0,j]
        # back: s_j a_j = Bk[0,j]; row1: r1 s_j a_{j+1} = Bk[1,j] -> a_{j+1} a_j = r1 Bk[1,j] Bk[0,j]
        if not back:
            for j in range(1,m): a[j]=a[j-1]*r1*Bk[1,j]*Bk[0,j]
            if a[m-1]*r1*Bk[1,0]*Bk[0,0]!=a[0]: continue
        else:
            for j in range(0,m-1): a[j+1]=a[j]*r1*Bk[1,j]*Bk[0,j]
            if a[0]*r1*Bk[1,m-1]*Bk[0,m-1]!=a[m-1]: continue
        # full check
        idx=((np.arange(m)[None,:]+(1 if back else -1)*np.arange(m)[:,None])%m)
        T=Bk*a[idx]  # should be r_i s_j rank-1
        if np.all(T*T[0:1,:]*T[:,0:1]*T[0,0]==1): sols.append(a)
    return sols
paf=lambda x: np.array([np.dot(x,np.roll(x,-s)) for s in range(len(x))])
def extract(N,b,off):
    H=load(N); v=N//4; m=v-b
    blocks=[]
    for c in range(4):
        Bk=H[off:off+m, off+c*m:off+(c+1)*m]
        s=seq_from_block(Bk,False) or seq_from_block(Bk,True)
        kind='circ' if seq_from_block(Bk,False) else 'back'
        blocks.append((s,kind))
    target = 0 if b==0 else -4
    good=[]
    for combo in itertools.product(*[s for s,_ in blocks]):
        X=np.array(combo); P=sum(paf(x) for x in X)
        if np.all(P[1:]==target): good.append(X)
    return blocks,good
if __name__=='__main__':
    todo=[(668,1,4),(716,1,4),(892,0,0),(1132,0,0),(1244,0,0),(1676,1,4),(1772,1,4),(1852,0,0),(1948,0,0),(1964,0,0)]
    for N,b,off in todo:
        blocks,good=extract(N,b,off)
        v=N//4
        print(N,v,'b',b,'kinds',[k for _,k in blocks],'cands',[len(s) for s,_ in blocks],'valid combos',len(good),flush=True)
        if good:
            X=good[0]
            # fix sign so row sums are >=0 (sign irrelevant to PAF)
            X=X*np.where(X.sum(1)<0,-1,1)[:,None]
            print('  row sums',X.sum(1))
            with open(f'data/core_v{v}_b{b}.txt','w') as f:
                f.write(f'v={v} N={N} b={b} m={v-b}\n')
                for k,x in zip('ABCD',X): f.write(k+'='+''.join('+' if t>0 else '-' for t in x)+'\n')
