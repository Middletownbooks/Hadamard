r"""Production run of the 647 x353 cross level: exact S-join + D-check over joint solutions [M], restart-safe.

Worker w of W takes the joint solutions of a priority list (logs/p647/prod/queue.json: orbit representatives of
y_B #0 under its stabiliser, best predicted first) at positions i = w (mod W).  For each joint solution:
  1. per-layer lists, distinct sigma, multiplicities (theory/p647_cross.Block);
  2. an optional restriction of the sigma lists (--rule; see theory/p647_restrict.py for the exact evaluation);
  3. the exact S-join (theory/p647_join.join4p), energy class by energy class; after each class its S-solutions
     are D-checked at once (theory/p647_cross.dcheck) and one line is appended to logs/p647/prod/w{w}_els.jsonl
     (joint solution, EL, pairs, S-solutions, D-candidates, D-hits, seconds): a restarted worker skips finished
     classes of a joint solution in progress and finished joint solutions (logs/p647/prod/w{w}.jsonl);
  4. every D-hit is assembled on Z_646 and checked by the full cyclic PAF; a core is written at once to
     logs/p647/CORE_w{w}_j{j}_{n}.json (and printed); completion and the independent verification are separate
     steps (theory/p647_finish.py, scripts/verify_h2588.py).
A worker never adopts another worker's checkpoints.  Stop: touch logs/p647/prod/STOP (checked between classes).
usage: python3 theory/p647_prod.py --worker W --workers N [--rule none|score:t:beta] [--queue FILE]"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p647_common as PC
import p647_cross as P
from p647_join import join4p
from nd_cross import orbit_maps

PROD = os.path.join(PC.LOG, "prod")

def restrict(blocks, rule):
    """keep masks per block; 'none' or 'score:t:beta' (keep sigma with mult >= t 2^(beta (e - e_pw)), e_pw = the
    pair-weighted mean energy of the block)"""
    if rule in (None, "", "none"): return [np.ones(B.N, bool) for B in blocks]
    kind, t, beta = rule.split(":"); t = float(t); beta = float(beta)
    assert kind == "score"
    out = []
    for B in blocks:
        e = B.Ssig[:, 0].astype(np.float64); pm = float((e * B.mult).sum() / B.mult.sum())
        out.append(B.mult >= t * np.exp2(beta * (e - pm)))
    return out

def run_one(X, yb, sol, j, w, rule, maps, NH, oi, log):
    elp = os.path.join(PROD, f"w{w}_els.jsonl")
    done = {}
    if os.path.exists(elp):
        for line in open(elp):
            try:
                r = json.loads(line)
                if r['j'] == j and r.get('rule') == rule: done[r['EL']] = r
            except Exception: pass
    t0 = time.time()
    blocks = P.build(X, yb, X.J.ya_full(sol['layers']), maps, NH, oi)
    keep = restrict(blocks, rule)
    idx = [np.nonzero(k_)[0] for k_ in keep]
    Ssig = [B.Ssig[i] for B, i in zip(blocks, idx)]
    TE = 2 * (X.TS + X.TD)
    _, st0 = join4p(Ssig, X.TS, stats_only=True)
    todo = [EL for EL in range(int(X.TS[0]) + 1) if EL not in done]
    info = {}; found_cores = []
    fo = open(elp, "a")
    def el_cb(EL, d): info[EL] = d
    def done_cb(EL, sols):
        d = info.get(EL, {})
        if sols:
            S = np.array([[idx[k][t[k]] for k in range(4)] for t in sols], np.int64)
            hits, cand = P.dcheck(blocks, S, TE)
        else:
            hits = np.zeros((0, 5), np.int64); cand = np.zeros(0, np.int64)
        for h in hits:
            xs = P.assemble(X, blocks, h[1:])
            ok = P.is_core(xs)
            if ok:
                path = os.path.join(PC.LOG, f"CORE_w{w}_j{j}_{len(found_cores)}.json")
                json.dump(dict(j=j, EL=EL, yb=[list(map(int, y)) for y in yb], ya=X.J.ya_full(sol['layers']),
                               pairs=[int(v) for v in h[1:]], core=xs.tolist()), open(path, "w"))
                found_cores.append(path)
                log(json.dumps(dict(stage="CORE", worker=w, j=j, EL=EL, path=path)))
        fo.write(json.dumps(dict(j=j, rule=rule, EL=EL, left=d.get('left', 0), right=d.get('right', 0), S=len(sols),
                                 cand=int(cand.sum()), hits=int(len(hits)), secs=d.get('secs', 0.0))) + "\n"); fo.flush()
        if os.path.exists(os.path.join(PROD, "STOP")): raise KeyboardInterrupt("STOP file")
    try:
        join4p(Ssig, X.TS, order=st0['order'], els=todo, el_cb=el_cb, done_cb=done_cb)
    finally:
        fo.close()
    # summary over all classes of this joint solution
    S = 0; cand = 0; hits = 0; left = 0; right = 0; secs = 0.0
    for line in open(elp):
        r = json.loads(line)
        if r['j'] == j and r.get('rule') == rule:
            S += r['S']; cand += r['cand']; hits += r['hits']; left += r['left']; right += r['right']; secs += r['secs']
    return dict(j=j, rule=rule, kept=[int(k_.sum()) for k_ in keep], N=[B.N for B in blocks], work=left + right,
                work_log2=round(math.log2(left + right), 3) if left + right else None, S=S, cand=cand,
                cand_log2=round(math.log2(cand), 3) if cand else None, hits=hits, cores=len(found_cores),
                join_secs=round(secs, 1), wall_secs=round(time.time() - t0, 1))

def main(w, W, rule, queue, log):
    os.makedirs(PROD, exist_ok=True)
    X, ys = PC.load_setup(); rec, sols = PC.joint_solutions(X, ys, 0)
    reps, oi, maps = orbit_maps(X); NH = P.lag_tensor(X, oi)
    q = json.load(open(queue))['order']
    mine = [j for i, j in enumerate(q) if i % W == w]
    sp = os.path.join(PROD, f"w{w}.jsonl")
    fin = set()
    if os.path.exists(sp):
        for line in open(sp):
            try:
                r = json.loads(line)
                if r.get('rule') == rule: fin.add(r['j'])
            except Exception: pass
    log(json.dumps(dict(stage="start", worker=w, workers=W, rule=rule, queue=len(mine), finished=len(fin))))
    for j in mine:
        if j in fin: continue
        if os.path.exists(os.path.join(PROD, "STOP")): log(json.dumps(dict(stage="stopped", worker=w))); return
        try:
            r = run_one(X, ys[0], sols[j], j, w, rule, maps, NH, oi, log)
        except KeyboardInterrupt:
            log(json.dumps(dict(stage="stopped", worker=w, j=j))); return
        with open(sp, "a") as f: f.write(json.dumps(r) + "\n")
        log(json.dumps(dict(stage="joint_done", worker=w, **r)))

if __name__ == "__main__":
    a = sys.argv[1:]
    opt = {"--worker": "0", "--workers": "1", "--rule": "none", "--queue": os.path.join(PROD, "queue.json")}
    i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: i += 1
    main(int(opt["--worker"]), int(opt["--workers"]), opt["--rule"], opt["--queue"], log=lambda s: print(s, flush=True))
