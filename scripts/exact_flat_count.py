import sys,numpy as np,json; sys.path.insert(0,'scripts')
from redlib import Red
v,d=int(sys.argv[1]),int(sys.argv[2]); R=Red(v,d)
assert R.r<=24
cnt=[]; Fs=[]
for start in range(0,2**R.r,2**18):
    idx=np.arange(start,min(2**R.r,start+2**18))
    Y=1-2*((idx[:,None]>>np.arange(R.r))&1)
    P=R.paf(R.expand(Y))[:,R.lags]; Fs.append((P**2).sum(1))
F=np.concatenate(Fs); ratio=np.sqrt(F/(R.m*R.L))
for rho in [0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.0]:
    print(rho, 'log2 #<=', round(float(np.log2(max(1,(ratio<=rho).sum()))),2))
print('min ratio',ratio.min())
