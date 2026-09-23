"""Flat-pool generator for the k-tree join: tabu walks on F=sum_c PAF(c)^2, emitting every
`stride`-th state with F<=F0 and |row sum| in the target set, as 64-bit packed orbit signs
(sign-normalised so the row sum is positive).  usage: kt_gen.py v d rho rowsums(csv) secs stride seed out.npy"""
import sys, time, numpy as np, numba as nb
sys.path.insert(0,'scripts')
from engine import Red, build, init_state, flip

@nb.njit(cache=True)
def walk(y,M,F0,nsteps,tenure,seed,stride,d,targets,out,nout):
    np.random.seed(seed)
    r=y.shape[1]; L=M.shape[0]
    G,P=init_state(y,M)
    F=0
    for c in range(L): F+=P[0,c]*P[0,c]
    rs=y[0,0]
    for q in range(1,r): rs+=d*y[0,q]
    tabu=np.zeros(r,np.int64); n=nout
    for step in range(nsteps):
        best=1<<62; bq=-1; nb_=0
        for q in range(r):
            if tabu[q]>step: continue
            yq=y[0,q]; dF=0
            for c in range(L):
                dd=-2*yq*G[0,c,q]+2*M[c,q,q]; dF+=2*P[0,c]*dd+dd*dd
            if dF<best or (dF==best and np.random.random()<1.0/(nb_+1)):
                if dF<best: nb_=1
                else: nb_+=1
                best=dF; bq=q
        yq=y[0,bq]
        for c in range(L):
            dd=-2*yq*G[0,c,bq]+2*M[c,bq,bq]; F+=2*P[0,c]*dd+dd*dd
        flip(y,G,P,M,0,bq); tabu[bq]=step+1+tenure+np.random.randint(tenure//2+1)
        rs+= -2*yq*(1 if bq==0 else d)
        if F<=F0 and step%stride==0:
            a=abs(rs); ok=False
            for t in targets:
                if a==t: ok=True
            if ok:
                sgn=1 if rs>0 else -1
                w=np.uint64(0)
                for q in range(r):
                    if y[0,q]*sgn>0: w|=np.uint64(1)<<np.uint64(q)
                out[n]=w; n+=1
                if n>=out.shape[0]: return n
    return n

if __name__=='__main__':
    v,d,rho=int(sys.argv[1]),int(sys.argv[2]),float(sys.argv[3])
    targets=np.array([int(x) for x in sys.argv[4].split(',')]); secs=float(sys.argv[5]); stride=int(sys.argv[6]); seed=int(sys.argv[7]); outf=sys.argv[8]
    R=Red(v,d); M=build(R); F0=int(rho**2*R.m*R.L); rng=np.random.default_rng(seed)
    import os; WALK=int(os.environ.get("WALK","300000"))
    buf=np.zeros(1<<24,np.uint64); n=0; t0=time.time(); it=0; ci=0; tot=0
    def flush(b):
        global ci, tot
        b=np.unique(b); ones=np.zeros(len(b),np.int64)
        for q in range(1,R.r): ones+=((b>>np.uint64(q))&np.uint64(1)).astype(np.int64)
        rs=np.where(b&np.uint64(1),1,-1)+d*(2*ones-(R.r-1))
        for a in targets:
            np.save(f"{outf}_rs{a}_c{ci:04d}.npy",b[rs==a])
        ci+=1; tot+=len(b)
    while time.time()-t0<secs:
        y=rng.choice([-1,1],size=(1,R.r)).astype(np.int64)
        n=walk(y,M,F0,WALK,max(3,R.r//10),seed*1_000_003+it,stride,d,targets,buf,n); it+=1
        if n>len(buf)-(1<<21):
            flush(buf[:n]); n=0
    flush(buf[:n])
    print(dict(chunk_entries=tot,walks=it,secs=round(time.time()-t0,1)),flush=True)
