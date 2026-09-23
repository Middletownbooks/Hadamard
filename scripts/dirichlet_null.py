"""CDF of T = sum_k (s_k - 1/4)^2 for s ~ Dirichlet(a,a,a,a), a=1 (complex characters)
or a=1/2 (real characters).  Monte Carlo + analytic small-ball tail, cached to data/."""
import numpy as np, os
from hadlib import ROOT
_cache={}
def _build(a, n=20_000_000, seed=1):
    f=os.path.join(ROOT,'data',f'dirT_a{a}.npy')
    if os.path.exists(f): return np.load(f)
    rng=np.random.default_rng(seed); out=[]
    for _ in range(n//2_000_000):
        g=rng.gamma(a,size=(2_000_000,4)); s=g/g.sum(1,keepdims=True)
        out.append(((s-0.25)**2).sum(1))
    T=np.sort(np.concatenate(out)).astype(np.float32); np.save(f,T); return T
def cdf(t, a=1.0):
    if a not in _cache: _cache[a]=_build(a)
    T=_cache[a]; n=len(T)
    t=np.atleast_1d(np.asarray(t,float))
    k=np.searchsorted(T,t.astype(np.float32),side='right')
    # small-ball regime: P = C t^{3/2}; C fitted on the 2000th order statistic
    t0=float(T[1999]); C=2000/n/t0**1.5
    if a==1.0: C=4*np.pi   # exact for t<=1/12
    return np.where(k<2000, C*t**1.5, k/n)
if __name__=='__main__':
    for a in (1.0,0.5):
        T=_build(a); print('a',a,'mean T',T.mean(),'-> per-share var',T.mean()/4)
        t0=float(T[1999]); print('  MC C',2000/len(T)/t0**1.5,' analytic Dir(1) 4pi=',4*np.pi)
