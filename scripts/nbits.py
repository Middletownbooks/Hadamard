"""Gaussian local estimate of n = -log2 P(random U-invariant 4-tuple is a b=0 solution).
P = sum over ordered signed row-sum 4-tuples with sum a^2 = 4v of prod P(row sum a_k) * p(R=0 | a),
p(R=0|a) = 4^(L-1) * Gaussian density at 0 of (R_1..R_{L-1}) with mean/cov from Monte Carlo of
blocks conditioned on their row sums (class L-1 dropped: implied by the profile).
Optionally restrict each block to 'flat' blocks (PAF ratio <= rho) -> n' (conditional on flatness)."""
import sys, json, numpy as np, itertools
sys.path.insert(0,'scripts')
from redlib import Red
from math import comb, log2
def block_stats(R,a,nsamp,rng,rho=None,pool=None):
    Y=R.sample(a,nsamp,rng) if pool is None else pool
    P=R.paf(R.expand(Y))[:,R.lags]
    if rho is not None:
        F=(P**2).sum(1); P=P[F<=rho**2*R.m*R.L]
    return P.mean(0), np.cov(P.T), len(P)
def nbits(v,d,nsamp=20000,seed=0,profiles=None):
    R=Red(v,d); rng=np.random.default_rng(seed)
    rs=[a for a in R.rowsums()]
    tot_blocks=2**R.r
    pa={a:R.count_rowsum(a)/tot_blocks for a in rs}
    stats={}
    Ptot=0.0; per={}
    for prof in (profiles or R.profiles(v)):
        # ordered signed assignments: permutations of multiset x signs (sign doesn't change PAF)
        perms=set(itertools.permutations(prof)); nsign=np.prod([2 if a!=0 else 1 for a in prof])
        mu=0; C=0
        for a in prof:
            if a not in stats: stats[a]=block_stats(R,a,nsamp,rng)[:2]
            mu=mu+stats[a][0]; C=C+stats[a][1]
        L=R.L; mu1=mu[:L-1]; C1=C[:L-1,:L-1]
        sign,logdet=np.linalg.slogdet(2*np.pi*C1)
        logdens=-0.5*logdet-0.5*mu1@np.linalg.solve(C1,mu1)
        lp=(L-1)*np.log(4)+logdens
        pprof=np.prod([pa[a] for a in prof])
        w=len(perms)*nsign*pprof*np.exp(lp)
        per[str(prof)]=dict(log2_p_given=float(lp/np.log(2)),weight_log2=float(np.log2(w)))
        Ptot+=w
    return -np.log2(Ptot), per
if __name__=='__main__':
    for arg in sys.argv[1:]:
        v,d=map(int,arg.split(','))
        n,per=nbits(v,d)
        print(json.dumps(dict(v=v,d=d,n_bits=round(float(n),1))),flush=True)
