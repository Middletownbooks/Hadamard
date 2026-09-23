"""Extract group-developed GS cores (b>=2) from normalised full matrices.
Group G = Z_mods (mixed radix, last index fastest), core offset 4b.
Block 0 is group-developed f(h-g), blocks 1..3 are f(g+h)."""
import numpy as np, gzip, itertools
def load(N): return np.loadtxt(gzip.open(f'data/hadamard_{N}.csv.gz'),delimiter=',',dtype=np.int8).astype(np.int64)
class G:
    def __init__(s,mods):
        s.mods=np.array(mods); s.m=int(np.prod(mods))
        s.coords=np.array(np.unravel_index(np.arange(s.m),mods)).T
    def idx(s,c):
        c=np.asarray(c)%s.mods
        return np.ravel_multi_index(tuple(np.moveaxis(c,-1,0)),s.mods)
    def add(s,i,j): return s.idx(s.coords[i]+s.coords[j])
    def sub(s,i,j): return s.idx(s.coords[i]-s.coords[j])
def solve(M,g,back):
    m=g.m; gens=[g.idx(np.eye(len(g.mods),dtype=int)[t]) for t in range(len(g.mods))]
    sols=[]
    for signs in itertools.product((1,-1),repeat=len(gens)):
        f=np.zeros(m,int); f[0]=1; stack=[0]; ok=True
        while stack and ok:
            j=stack.pop()
            for e,r in zip(gens,signs):
                # back: f(j+e)=f(j) r M[e,j]M[0,j]; circ: f(j-e)=f(j) r M[e,j]M[0,j]
                k = g.idx(g.coords[j]+g.coords[e]) if back else g.idx(g.coords[j]-g.coords[e])
                val=f[j]*r*M[e,j]*M[0,j]
                if f[k]==0: f[k]=val; stack.append(k)
                elif f[k]!=val: ok=False;break
        if not ok or np.any(f==0): continue
        I=np.arange(m)
        P = np.array([[0]]) 
        F = f[g.idx(g.coords[I][:,None,:]+g.coords[I][None,:,:])] if back else f[g.idx(g.coords[I][None,:,:]-g.coords[I][:,None,:])]
        T=M*F
        if np.all(T*T[0:1,:]*T[:,0:1]*T[0,0]==1): sols.append(f)
    return sols
def paf_group(g,f):
    # PAF(s)=sum_i f(i) f(i+s)
    return np.array([np.dot(f, f[g.idx(g.coords+g.coords[s])]) for s in range(g.m)])
if __name__=='__main__':
    for N,b,mods in [(1388,5,(2,3,57)),(1436,7,(2,2,8,11)),(1916,3,(2,2,119))]:
        H=load(N); v=N//4; m=v-b; off=4*b; g=G(mods)
        cands=[solve(H[off:off+m,off+c*m:off+(c+1)*m],g,c>0) for c in range(4)]
        print(N,v,b,mods,[len(c) for c in cands],flush=True)
        pafs=[[paf_group(g,f) for f in c] for c in cands]
        good=[]
        for combo in itertools.product(*[range(len(c)) for c in cands]):
            T=sum(pafs[k][combo[k]] for k in range(4))
            vals=sorted(set(T[1:].tolist()))
            if len(vals)<=2: good.append((combo,vals,T))
        print('  combos two-valued:',len(good))
        for combo,vals,T in good[:2]:
            X=np.array([cands[k][combo[k]] for k in range(4)])
            Hs=[s for s in range(1,m) if T[s]==-4*b]
            print('  values',vals,'#(-4b)',len(Hs),'+1 =',len(Hs)+1,'= m/(b+1)?',m//(b+1),'row sums',X.sum(1))
        combo,vals,T=good[0]; X=np.array([cands[k][combo[k]] for k in range(4)])
        X=X*np.where(X.sum(1)<0,-1,1)[:,None]
        with open(f'data/core_v{v}_b{b}.txt','w') as fh:
            fh.write(f'v={v} N={N} b={b} m={m} group={"x".join(map(str,mods))} (mixed radix, last index fastest)\n')
            for k,x in zip('ABCD',X): fh.write(k+'='+''.join('+' if t>0 else '-' for t in x)+'\n')
