r"""Production driver for theory/q523_wagner.py: W worker processes (seeds base..base+W-1), each looping over middle
values until the deadline or until any worker has written a verified core (stop file).  Per-worker JSON-line logs.
Every S-join is checkpointed (logs/q523/prod/ckpt/, pair indices + D-check progress); with --resume the workers first
finish the D-checks of checkpoints left by an interrupted run on the same y.
usage: python3 theory/q523_prod.py F q s y.json --tag name [--workers 4] [--hours 2] [--logN 32] [--K 24] [--rows 64]
          [--beta 1] [--seed 100] [--planted] [--resume]"""
import sys, os, json, time, glob, hashlib, multiprocessing as mp
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def worker(args):
    F, q, s, yfile, tag, w, seed, hours, logN, K, rows, beta, planted, resume, nw = args
    from q523_sdsolve import SDInst, planted_target
    import q523_wagner as W
    Y = json.load(open(yfile)); Y = Y['y'] if isinstance(Y, dict) else Y
    inst = SDInst(F, q, s, Y)
    if planted:
        Pt, _ = planted_target(inst, np.random.default_rng(seed + 999)); inst = SDInst(F, q, s, Y, T=Pt)
    d = os.path.join(ROOT, 'logs', 'q523', 'prod'); os.makedirs(d, exist_ok=True)
    lf = open(os.path.join(d, f"{tag}_w{w}.jsonl"), 'a')
    def log(rec):
        rec = dict(rec, worker=w, t=round(time.time(), 1)); lf.write(json.dumps(rec) + "\n"); lf.flush()
    out = os.path.join(ROOT, 'data', f"core_v{F * q + 1}_b1_q523_{tag}_w{w}.txt")
    stop = os.path.join(d, f"{tag}.STOP")
    log(dict(stage="start", seed=seed, logN=logN, K=K, rows=rows, beta=beta, planted=planted))
    yh = hashlib.sha1(json.dumps(Y).encode()).hexdigest()[:8] + ("p" if planted else "")
    cdir = os.path.join(d, 'ckpt'); os.makedirs(cdir, exist_ok=True)
    if resume:
        files = sorted(glob.glob(os.path.join(cdir, f"*_y{yh}_c*.npy")))[w::nw]
        if files:
            from q523_dcheck_nb import DCheckNB
            dchk = DCheckNB(inst, seed)
            for f in files:
                if os.path.exists(stop): break
                if W.resume_ckpt(inst, f, dchk, out, log, stop, seed, beta): return dict(cores=1)
    found, st = W.run(inst, beta, rows, K, logN, 10 ** 6, seed, 20, out, log, None, time.time() + 3600 * hours, stop, predict=True,
                      ckpt=os.path.join(cdir, f"{tag}_w{w}_y{yh}"))
    log(dict(stage="done", **st))
    return st

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--tag": "run", "--workers": "4", "--hours": "2", "--logN": "32", "--K": "24", "--rows": "64",
                             "--beta": "1", "--seed": "100"}
    pos = []; flags = set(); i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        elif a[i].startswith('--'): flags.add(a[i]); i += 1
        else: pos.append(a[i]); i += 1
    F, q, s = map(int, pos[:3]); W_ = int(opt["--workers"])
    jobs = [(F, q, s, pos[3], opt["--tag"], w, int(opt["--seed"]) + w, float(opt["--hours"]), float(opt["--logN"]),
             int(opt["--K"]), int(opt["--rows"]), float(opt["--beta"]), '--planted' in flags, '--resume' in flags, W_)
            for w in range(W_)]
    with mp.get_context("fork").Pool(W_) as pool:
        res = pool.map(worker, jobs)
    print(json.dumps(res))
