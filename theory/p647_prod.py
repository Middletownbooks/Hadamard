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
Class selection (--els K): only the K left-energy classes nearest the pair-weighted energy centre are joined
(section 3 of the report: K = 2 keeps about 55% of the D-candidates for 20% of the work).
Unit mode (--units): the class-merged sub-bicliques of logs/p647/prod/units.json (all lefts of a unit complete all its
rights at stage 1A, so one join4m covers every (left, right) joint solution; 410 units cover the 3,488 orbits with
2.78x fewer emitted sides than one join per orbit); checkpoints u{w}_els.jsonl and u{w}.jsonl.
usage: python3 theory/p647_prod.py --worker W --workers N [--els K] [--kernel p|q] [--rule none|score:t:beta] [--queue FILE]
       python3 theory/p647_prod.py --units --worker W --workers N --els 2"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p647_common as PC
import p647_cross as P
from p647_join import join4p, join4q, SL, pair_counts, hash_w, res_w
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

def central_els(blocks, Ssig, order, K, E0=645):
    """the K left-energy classes (with work on both sides) nearest the centre ((e_a + e_b) + (E0 - e_c - e_d)) / 2
    of the pair-weighted mean energies (report section 3: they hold about 55% of the yield for 20% of the work at K = 2)"""
    a, b, c, d = order
    pm = [float((B.Ssig[:, 0] * B.mult).sum() / B.mult.sum()) for B in blocks]
    centre = ((pm[a] + pm[b]) + (E0 - pm[c] - pm[d])) / 2
    w = hash_w(20, 11); u = res_w(20, 61, 11)
    Ls = [SL(s_, E0, w, u, 61) for s_ in Ssig]
    cL = pair_counts(Ls[a], Ls[b], E0, 61).sum(1); cR = pair_counts(Ls[c], Ls[d], E0, 61).sum(1)[::-1]
    live = np.nonzero((cL > 0) & (cR > 0))[0]
    pick = sorted(live, key=lambda e: abs(e - centre))[:K]
    return sorted(int(e) for e in pick), round(centre, 1)

def run_one(X, yb, sol, j, w, rule, maps, NH, oi, log, K=0, kernel="p"):
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
    keep = restrict(blocks, rule.split("|")[0])
    idx = [np.nonzero(k_)[0] for k_ in keep]
    Ssig = [B.Ssig[i] for B, i in zip(blocks, idx)]
    TE = 2 * (X.TS + X.TD)
    jf = join4q if kernel == "q" else join4p
    _, st0 = jf(Ssig, X.TS, stats_only=True)
    if K > 0:
        els, centre = central_els(blocks, Ssig, st0['order'], K)
    else:
        els, centre = list(range(int(X.TS[0]) + 1)), None
    todo = [EL for EL in els if EL not in done]
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
        jf(Ssig, X.TS, order=st0['order'], els=todo, el_cb=el_cb, done_cb=done_cb)
    finally:
        fo.close()
    # summary over all classes of this joint solution
    S = 0; cand = 0; hits = 0; left = 0; right = 0; secs = 0.0
    for line in open(elp):
        r = json.loads(line)
        if r['j'] == j and r.get('rule') == rule:
            S += r['S']; cand += r['cand']; hits += r['hits']; left += r['left']; right += r['right']; secs += r['secs']
    return dict(j=j, rule=rule, els=els if K > 0 else "all", centre=centre, kept=[int(k_.sum()) for k_ in keep], N=[B.N for B in blocks],
                work_all=st0['left_total'] + st0['right_total'], work=left + right,
                work_log2=round(math.log2(left + right), 3) if left + right else None, S=S, cand=cand,
                cand_log2=round(math.log2(cand), 3) if cand else None, hits=hits, cores=len(found_cores),
                join_secs=round(secs, 1), wall_secs=round(time.time() - t0, 1))

def main(w, W, rule, queue, log, K=0, kernel="p"):
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
            r = run_one(X, ys[0], sols[j], j, w, rule, maps, NH, oi, log, K=K, kernel=kernel)
        except KeyboardInterrupt:
            log(json.dumps(dict(stage="stopped", worker=w, j=j))); return
        with open(sp, "a") as f: f.write(json.dumps(r) + "\n")
        log(json.dumps(dict(stage="joint_done", worker=w, **r)))

# ---------------------------------------------------------------- unit mode: class-merged sub-bicliques (join4m)
from p647_join import join4m
from p647_deeptail import nested_counts

def deep_counts(X, bl, S, TE, cross):
    """exact D-candidates hitting the nested cross-lag subsets (live deep-tail calibration, report section 2.3)"""
    cnt = np.zeros(len(cross), np.int64); smp = np.zeros((1, len(TE)), np.int64)
    B = bl
    nested_counts(np.ascontiguousarray(S), B[0].pairs, B[0].bd, B[0].n1, B[1].pairs, B[1].bd, B[1].n1, B[2].pairs, B[2].bd,
                  B[2].n1, B[3].pairs, B[3].bd, B[3].n1, B[0].P0, B[0].P1, B[1].P0, B[1].P1, B[2].P0, B[2].P1, B[3].P0,
                  B[3].P1, TE, cross, cnt, smp, 1 << 62, 0)
    return cnt

def block_for(X, yb, sol, k, maps, NH, cache):
    """the Block of block k for the y_A option of joint solution sol (cached by (k, option))"""
    key = (k, json.dumps(sol['layers'][k]))
    if key not in cache:
        ya = X.J.ya_full(sol['layers'])
        L0 = X.layer_labellings(0, yb[k], ya[k]); L1 = X.layer_labellings(1, yb[k], ya[k])
        cache[key] = P.Block(L0, L1, maps, NH)
    return cache[key]

def run_unit(X, yb, sols, u, uid, w, K, maps, NH, log, E0=645):
    JJ = u['joint']; nl, nr = len(JJ), len(JJ[0])
    elp = os.path.join(PROD, f"u{w}_els.jsonl")
    done = set()
    if os.path.exists(elp):
        for line in open(elp):
            try:
                r = json.loads(line)
                if r['unit'] == uid: done.add(r['EL'])
            except Exception: pass
    t0 = time.time(); cache = {}
    bl = {(i, j): [block_for(X, yb, sols[JJ[i][j]], k, maps, NH, cache) for k in range(4)] for i in range(nl) for j in range(nr)}
    lefts = [(bl[(i, 0)][0].Ssig, bl[(i, 0)][1].Ssig) for i in range(nl)]
    rights = [(bl[(0, j)][2].Ssig, bl[(0, j)][3].Ssig) for j in range(nr)]
    for i in range(nl):
        for j in range(nr):
            assert bl[(i, j)][0] is bl[(i, 0)][0] and bl[(i, j)][1] is bl[(i, 0)][1]
            assert bl[(i, j)][2] is bl[(0, j)][2] and bl[(i, j)][3] is bl[(0, j)][3]
    build_secs = time.time() - t0
    _, st0 = join4m(lefts, rights, X.TS, stats_only=True)
    pm = lambda B: float((B.Ssig[:, 0] * B.mult).sum() / B.mult.sum())
    centre = float(np.mean([((pm(bl[(i, j)][0]) + pm(bl[(i, j)][1])) + (E0 - pm(bl[(i, j)][2]) - pm(bl[(i, j)][3]))) / 2
                            for i in range(nl) for j in range(nr)]))
    live = sorted(st0['per_EL_work'])
    els = sorted(sorted(live, key=lambda e: abs(e - centre))[:K]) if K > 0 else live
    TE = 2 * (X.TS + X.TD)
    cross = np.array([t for t, h in enumerate(X.lags) if h % 17 and h % 19], np.int64)
    info = {}; cores = []
    fo = open(elp, "a")
    def el_cb(EL, d): info[EL] = d
    def done_cb(EL, found):
        d = info.get(EL, {})
        per = {}
        for (i, j, s0, s1, s2, s3) in found: per.setdefault((i, j), []).append((s0, s1, s2, s3))
        cand = {}; hits = {}; nested = np.zeros(len(cross), np.int64)
        for (i, j), S in per.items():
            h, c = P.dcheck(bl[(i, j)], np.array(S, np.int64), TE)
            cand[f"{i},{j}"] = int(c.sum()); hits[f"{i},{j}"] = int(len(h))
            nested += deep_counts(X, bl[(i, j)], np.array(S, np.int64), TE, cross)
            for hh in h:
                xs = P.assemble(X, bl[(i, j)], hh[1:])
                if P.is_core(xs):
                    e = JJ[i][j]
                    path = os.path.join(PC.LOG, f"CORE_u{w}_unit{uid}_j{e}_{len(cores)}.json")
                    json.dump(dict(j=e, unit=uid, EL=EL, yb=[list(map(int, y)) for y in yb], ya=X.J.ya_full(sols[e]['layers']),
                                   pairs=[int(v) for v in hh[1:]], core=xs.tolist()), open(path, "w"))
                    cores.append(path); log(json.dumps(dict(stage="CORE", worker=w, unit=uid, j=e, EL=EL, path=path)))
        fo.write(json.dumps(dict(unit=uid, EL=EL, left=d.get('left', 0), right=d.get('right', 0), S=len(found),
                                 cand=int(sum(cand.values())), hits=int(sum(hits.values())), secs=d.get('secs', 0.0),
                                 joint={f"{i},{j}": JJ[i][j] for i in range(nl) for j in range(nr)}, cand_ij=cand,
                                 nested=nested.tolist())) + "\n"); fo.flush()
        if os.path.exists(os.path.join(PROD, "STOP")): raise KeyboardInterrupt("STOP file")
    try:
        join4m(lefts, rights, X.TS, els=[e for e in els if e not in done], el_cb=el_cb, done_cb=done_cb)
    finally:
        fo.close()
    S = 0; cand = 0; hits = 0; work = 0; secs = 0.0
    for line in open(elp):
        r = json.loads(line)
        if r['unit'] == uid: S += r['S']; cand += r['cand']; hits += r['hits']; work += r['left'] + r['right']; secs += r['secs']
    return dict(unit=uid, shape=[nl, nr], orbits=len(u['orbits']), H=u['H'], els=els, centre=round(centre, 1),
                work_all=st0['left_total'] + st0['right_total'], work=work, work_log2=round(math.log2(work), 3) if work else None,
                S=S, cand=cand, cand_log2=round(math.log2(cand), 3) if cand else None, hits=hits, cores=len(cores),
                build_secs=round(build_secs, 1), join_secs=round(secs, 1), wall_secs=round(time.time() - t0, 1))

def main_units(w, W, K, log, units_path=None):
    os.makedirs(PROD, exist_ok=True)
    X, ys = PC.load_setup(); rec, sols = PC.joint_solutions(X, ys, 0)
    reps, oi, maps = orbit_maps(X); NH = P.lag_tensor(X, oi)
    U = json.load(open(units_path or os.path.join(PROD, "units.json")))['units']
    sp = os.path.join(PROD, f"u{w}.jsonl"); fin = set()
    if os.path.exists(sp):
        for line in open(sp):
            try: fin.add(json.loads(line)['unit'])
            except Exception: pass
    mine = [uid for uid in range(len(U)) if uid % W == w]
    log(json.dumps(dict(stage="start_units", worker=w, workers=W, K=K, units=len(mine), finished=len(fin))))
    for uid in mine:
        if uid in fin: continue
        if os.path.exists(os.path.join(PROD, "STOP")): log(json.dumps(dict(stage="stopped", worker=w))); return
        try:
            r = run_unit(X, ys[0], sols, U[uid], uid, w, K, maps, NH, log)
        except KeyboardInterrupt:
            log(json.dumps(dict(stage="stopped", worker=w, unit=uid))); return
        with open(sp, "a") as f: f.write(json.dumps(r) + "\n")
        log(json.dumps(dict(stage="unit_done", worker=w, **r)))

if __name__ == "__main__":
    a = sys.argv[1:]
    opt = {"--worker": "0", "--workers": "1", "--rule": "none", "--queue": os.path.join(PROD, "queue.json"), "--els": "0", "--kernel": "p"}
    i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: i += 1
    if "--units" in a:
        main_units(int(opt["--worker"]), int(opt["--workers"]), int(opt["--els"]), log=lambda s: print(s, flush=True))
        sys.exit()
    rule = opt["--rule"] if int(opt["--els"]) == 0 else f"{opt['--rule']}|els{opt['--els']}"
    main(int(opt["--worker"]), int(opt["--workers"]), rule, opt["--queue"], log=lambda s: print(s, flush=True),
         K=int(opt["--els"]), kernel=opt["--kernel"])
