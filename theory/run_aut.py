"""Run Traces automorphism computation on a stored Hadamard matrix; write generators."""
import numpy as np, gzip, subprocess, time, sys
n=int(sys.argv[1]); out=sys.argv[2]
H=np.loadtxt(gzip.open(f'/home/user/Hadamard/data/hadamard_{n}.csv.gz'),delimiter=',',dtype=np.int64)
assert H.shape==(n,n) and np.array_equal(H@H.T, n*np.eye(n,dtype=np.int64)); print(n,'verified',flush=True)
inp=(f"{n}\n"+"\n".join(" ".join(map(str,r)) for r in H)).encode()
t=time.time(); r=subprocess.run(['/home/user/Hadamard/theory/hadaut_traces',out],input=inp,capture_output=True)
print(r.stderr.decode(), round(time.time()-t,1),'s',flush=True)
