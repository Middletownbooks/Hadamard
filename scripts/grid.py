"""generic background grid runner: python3 grid.py jobs.jsonl out.jsonl nworkers
each job: {"kind": "s3sched"|"hinge"|"floor", ...params}"""
import sys, json, time, numpy as np, multiprocessing as mp
sys.path.insert(0,'scripts')
def work(job):
    from engine import Red, build, random_y, run_joint_sched, run_hinge, flat_tabu, fingerprint
    R=Red(job['v'],job['d']); M=build(R); rng=np.random.default_rng(job['seed']+1000*job['v'])
    t=time.time(); out=dict(job)
    if job['kind']=='s3sched':
        y=random_y(R,rng); st,ok=run_joint_sched(y,M,job['lam0'],job['tau'],job['lam_min'],job['maxsteps'],job['tenure'],job['seed'])
    elif job['kind']=='hinge':
        y=random_y(R,rng); F0=int(job['rho']**2*R.m*R.L); st,ok=run_hinge(y,M,job['lam'],F0,job['maxsteps'],job['tenure'],job['seed'])
    elif job['kind']=='floor':
        y=random_y(R,rng)[:1]; flat_tabu(y,M,job['maxsteps'],job['tenure'],job['seed'],0); st,ok=job['maxsteps'],False
        X=R.expand(y); pa=R.paf(X)[0,1:]; out.update(paf_std=float(pa.std()),ratio=float(pa.std()/np.sqrt(R.m)),rowsum=int(X.sum()))
    out['secs']=round(time.time()-t,2)
    if job['kind']!='floor':
        out.update(steps=int(st),solved=bool(ok))
        if ok:
            sv,ratio=fingerprint(R,y); out.update(sv=sv,ratio=ratio,blocks_ratio=[float(z) for z in R.paf(R.expand(y))[:,1:].std(1)/np.sqrt(R.m)],
                                                  Y=[''.join('+' if t>0 else '-' for t in yy) for yy in y])
    return out
if __name__=='__main__':
    jobs=[json.loads(l) for l in open(sys.argv[1])]
    with mp.Pool(int(sys.argv[3])) as pool, open(sys.argv[2],'a') as f:
        for res in pool.imap_unordered(work,jobs):
            f.write(json.dumps(res)+'\n'); f.flush()
