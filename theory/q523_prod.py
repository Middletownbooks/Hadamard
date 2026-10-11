r"""Production driver for theory/q523_wagner.py: W worker processes (seeds base..base+W-1), each looping over middle
values until the deadline or until any worker has written a verified core (stop file).  Per-worker JSON-line logs.
Every S-join is checkpointed (logs/q523/prod/ckpt/: hits and next pass during the join, then pair indices and D-check
progress); with --resume the workers first finish the D-checks left by an interrupted run on the same y, and a worker
restarted with the same command line (tag, seed, lists) continues its interrupted S-join.
usage: python3 theory/q523_prod.py F q s y.json --tag name [--workers 4] [--hours 2] [--logN 32] [--K 24] [--rows 64]
          [--beta 1] [--seed 100] [--planted] [--resume] [--big] [--pass 23]
--big: lazy parents (theory/q523_wagner2.py) for large --rows (up to 256) and logN; --pass: log2 left pairs per pass"""
import sys, os, json, time, glob, hashlib, multiprocessing as mp
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def worker(args):
    F, q, s, yfile, tag, w, seed, hours, logN, K, rows, beta, planted, resume, nw, big, plog = args
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
    cpre = os.path.join(cdir, f"{tag}_w{w}_y{yh}")
    if resume:
        # this worker's own interrupted D-checks (its earlier process is gone), then orphaned D-check checkpoints of
        # other runs: a live worker rewrites its progress file every 2000 checks (seconds), so idle > 5 min means dead
        own = sorted(glob.glob(cpre + "_c*.npy"), key=lambda f: int(f[len(cpre) + 2:-4]))
        other = sorted(f for f in glob.glob(os.path.join(cdir, f"*_y{yh}_c*.npy"))
                       if not os.path.basename(f).startswith(f"{tag}_w")
                       and time.time() - os.path.getmtime(f[:-4] + '.json' if os.path.exists(f[:-4] + '.json') else f) > 300)[w::nw]
        if own or other:
            from q523_dcheck2 import DCheck2
            dchk = DCheck2(inst, seed)
            for f in own + other:
                if os.path.exists(stop): break
                if W.resume_ckpt(inst, f, dchk, out, log, stop, seed, beta): return dict(cores=1)
    start = 0
    if resume:
        # same command line: continue after the last centre whose S-join this worker logged, or inside an
        # interrupted S-join (its checkpoint)
        for l in open(os.path.join(d, f"{tag}_w{w}.jsonl")):
            try: rec = json.loads(l)
            except ValueError: continue
            if rec.get('stage') == 'S-join': start = max(start, int(rec['centre']) + 1)
        for j in glob.glob(cpre + "_c*_join.npz"): start = max(start, int(j[len(cpre) + 2:].split('_')[0]))
    found, st = W.run(inst, beta, rows, K, logN, 10 ** 6, seed, 20, out, log, None, time.time() + 3600 * hours, stop, predict=True,
                      ckpt=cpre, big_middle=big, pass_log2=plog, start_centre=start)
    log(dict(stage="done", **st))
    return st

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--tag": "run", "--workers": "4", "--hours": "2", "--logN": "32", "--K": "24", "--rows": "64",
                             "--beta": "1", "--seed": "100", "--pass": "23"}
    pos = []; flags = set(); i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        elif a[i].startswith('--'): flags.add(a[i]); i += 1
        else: pos.append(a[i]); i += 1
    F, q, s = map(int, pos[:3]); W_ = int(opt["--workers"])
    jobs = [(F, q, s, pos[3], opt["--tag"], w, int(opt["--seed"]) + w, float(opt["--hours"]), float(opt["--logN"]),
             int(opt["--K"]), int(opt["--rows"]), float(opt["--beta"]), '--planted' in flags, '--resume' in flags, W_,
             '--big' in flags, float(opt["--pass"])) for w in range(W_)]
    # independent processes (not a Pool): a worker killed by the OOM killer or a signal ends alone, the others go on
    ctx = mp.get_context("fork"); procs = [ctx.Process(target=worker, args=(j,)) for j in jobs]
    for p_ in procs: p_.start()
    for p_ in procs: p_.join()
    print(json.dumps([dict(worker=w, exitcode=p_.exitcode) for w, p_ in enumerate(procs)]))
