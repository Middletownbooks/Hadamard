"""Common utilities: load cores, group arithmetic, PAF/PSD, multiplier stabilisers."""
import numpy as np, re, glob, os
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

class Core:
    def __init__(s, path):
        S={}; hdr=None
        for l in open(path):
            if hdr is None: hdr=l.strip(); continue
            if len(l)>2 and l[0] in 'ABCD' and l[1]=='=':
                S[l[0]]=np.array([1 if c=='+' else -1 for c in l[2:].strip()])
        s.X=np.array([S[k] for k in 'ABCD'])
        s.m=s.X.shape[1]
        kv=dict(re.findall(r'(\w+)=(\S+)',hdr))
        s.v=int(kv['v']); s.b=int(kv.get('b',0))
        s.mods=tuple(int(t) for t in kv['group'].split('x')) if 'group' in kv else (s.m,)
        s.name=os.path.basename(path)
    # group helpers: element index <-> coords (mixed radix, last fastest)
    @property
    def shape(s): return s.mods
    def arr(s,x): return x.reshape(s.mods)
    def psd(s,X=None):
        X=s.X if X is None else X
        F=np.fft.fftn(X.reshape((X.shape[0],)+s.mods).astype(float),axes=tuple(range(1,len(s.mods)+1)))
        return (np.abs(F)**2).reshape(X.shape[0],-1)
    def paf(s,X=None):
        P=s.psd(X)
        R=np.fft.ifftn(P.reshape((P.shape[0],)+s.mods),axes=tuple(range(1,len(s.mods)+1))).real
        return np.rint(R.reshape(P.shape[0],-1)).astype(np.int64)
    def target(s):
        """required sum of PAF at each nonzero group element"""
        m=s.m; t=np.zeros(m,dtype=np.int64); t[0]=4*m
        if s.b==0: return t
        if s.b==1: t[1:]=-4; return t
        # b>=2: -4b on H\{e}, +4 off H, H = subgroup of index b+1; identify H from the data
        R=s.paf().sum(0); return R
    def verify(s):
        R=s.paf().sum(0)
        if s.b==0: return bool(np.all(R[1:]==0) and R[0]==4*s.m)
        if s.b==1: return bool(np.all(R[1:]==-4))
        vals=set(R[1:].tolist()); nH=int(np.sum(R[1:]==-4*s.b))+1
        return vals=={-4*s.b,4} and nH*(s.b+1)==s.m
    def coords(s):
        return np.array(np.unravel_index(np.arange(s.m),s.mods)).T
    def idx(s,c):
        c=np.asarray(c)%np.array(s.mods)
        return np.ravel_multi_index(tuple(np.moveaxis(c,-1,0)),s.mods)
    def units(s):
        """integers u in [1,e) coprime to exponent e: scalar multipliers"""
        from math import gcd
        e=int(np.lcm.reduce(s.mods)); return [u for u in range(1,e) if gcd(u,e)==1]

def all_cores():
    return [Core(p) for p in sorted(glob.glob(os.path.join(ROOT,'data','core_v*_b*.txt')),key=lambda p:int(re.search(r'core_v(\d+)',p).group(1)))]

def translate_match(x, y, core):
    """all translations c with x(i+c)=y(i) (group translations); returns list of element indices"""
    # correlation via FFT
    X=np.fft.fftn(core.arr(x).astype(float)); Y=np.fft.fftn(core.arr(y).astype(float))
    C=np.fft.ifftn(np.conj(Y)*X).real.reshape(-1)   # C[c]=sum_i y(i) x(i+c)
    return list(np.where(np.rint(C)==core.m)[0])

def stabiliser(core, k):
    """multipliers u (scalar) such that block k is invariant under i -> u*i + c for some c"""
    x=core.X[k]; co=core.coords(); out=[]
    for u in core.units():
        perm=core.idx(co*u)       # i -> u*i
        y=x[perm]                 # y(i)=x(u i)
        cs=translate_match(x,y,core)  # x(i+c)=x(u i)
        if cs: out.append((u,cs))
    return out
