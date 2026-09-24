"""Cost model for b=1 with the half-shift split (m=2h, h prime, d=1).
x_k on Z_2h -> u_k(j)=(x(j)+x(j+h))/2, w_k=(x(j)-x(j+h))/2 on Z_h, ternary, complementary supports S_k / S_k^c.
Constraints split: sum_k PAF_{u_k}(s) = const_u and sum_k PAF_{w_k}(s) = const_w for s != 0 in Z_h
((h-1)/2 independent lag classes each).  Gaussian local estimate of bits per half from Monte Carlo of random
supports (|S_k| ~ h/2 as in the real cores) and random signs; lattice spacing measured from the samples.
Freedom: sign bits on supports (sum|S_k| ~ 2h-1 total) + support choice (<= 4h bits)."""
import numpy as np, sys, json
from math import log2, lgamma, log
def paf(x):
    F=np.fft.fft(x,axis=-1); return np.rint(np.fft.ifft(np.abs(F)**2,axis=-1).real).astype(np.int64)
def lnC(n,k): return (lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1))/log(2)
rng=np.random.default_rng(0)
for v in [167,179,515,719]:
    m=v-1; h=m//2
    sizes=[h//2]*4                             # typical |S_k|; real cores: 39-44 (h=83), 43-46 (h=89)
    N=20000
    S=np.zeros((N,h))
    for k in range(4):
        sup=np.argsort(rng.random((N,h)),axis=1)[:,:sizes[k]]
        u=np.zeros((N,h)); np.put_along_axis(u,sup,rng.choice([-1,1],size=(N,sizes[k])),axis=1)
        S+=paf(u)
    L=(h-1)//2; vals=S[:,1:L+1]
    sd=vals.std(); g=np.gcd.reduce(np.abs(np.diff(np.unique(vals))).astype(int)) if True else 1
    bits_class=log2(np.sqrt(2*np.pi)*sd/g)
    n_half=L*bits_class
    free_signs=sum(sizes)                      # sign bits of the u-half (w-half similar with complements)
    free_supports=sum(lnC(h,s) for s in sizes) # support choice (shared by both halves)
    print(json.dumps(dict(v=v,h=h,classes_per_half=L,class_sd=round(float(sd),1),lattice=int(g),bits_per_class=round(bits_class,2),
        n_per_half=round(n_half,0),n_total=round(2*n_half,0),sign_bits_per_half=free_signs,support_bits=round(free_supports,0))))
