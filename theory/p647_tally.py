r"""Status of the 647 x353 production run (logs/p647/prod/) and of the calibration runs (logs/p647/calib/).

Per worker it sums the per-class checkpoint lines (w*_els.jsonl): joint solutions finished and in progress, emitted
pairs, S-solutions, D-candidates (= pair-level S-solutions, prod_k mult(sigma_k) summed over S-solutions) and D-hits;
and over everything the expected number of cores  sum D-candidates * 2^-C(D|S)  with C(D|S) given (calibrated value
from report/paley_647_answer.md section 2; default 34.6 bits).  Cores found: logs/p647/CORE_*.json.
usage: python3 theory/p647_tally.py [--cdgs 34.6] [--json]"""
import sys, os, json, glob, math
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROD = os.path.join(ROOT, "logs", "p647", "prod"); CAL = os.path.join(ROOT, "logs", "p647", "calib")

def main(cdgs=34.6, as_json=False):
    out = dict(workers={}, total={})
    T = dict(joint_done=0, joint_partial=0, pairs=0, S=0, cand=0, hits=0, secs=0.0)
    for f in sorted(glob.glob(os.path.join(PROD, "w*_els.jsonl"))):
        w = os.path.basename(f).split("_")[0]
        per = {}
        for line in open(f):
            try: r = json.loads(line)
            except Exception: continue
            k = (r['j'], r.get('rule'))
            d = per.setdefault(k, dict(pairs=0, S=0, cand=0, hits=0, secs=0.0, els=0))
            d['pairs'] += r['left'] + r['right']; d['S'] += r['S']; d['cand'] += r['cand']; d['hits'] += r['hits']
            d['secs'] += r['secs']; d['els'] += 1
        done = set()
        sp = os.path.join(PROD, f"{w}.jsonl")
        if os.path.exists(sp):
            for line in open(sp):
                try: r = json.loads(line); done.add((r['j'], r.get('rule')))
                except Exception: pass
        W = dict(joint_done=len(done), joint_partial=len(set(per) - done), pairs=sum(d['pairs'] for d in per.values()),
                 S=sum(d['S'] for d in per.values()), cand=sum(d['cand'] for d in per.values()),
                 hits=sum(d['hits'] for d in per.values()), secs=round(sum(d['secs'] for d in per.values()), 1))
        out['workers'][w] = W
        for k in T: T[k] += W[k]
    # calibration runs (unrestricted, complete) count as checked joint solutions too
    C = dict(joint_done=0, S=0, cand=0, hits=0)
    cp = os.path.join(CAL, "calib.jsonl")
    if os.path.exists(cp):
        for line in open(cp):
            r = json.loads(line); C['joint_done'] += 1; C['S'] += r.get('S_solutions', 0); C['cand'] += r.get('S_pairs', 0); C['hits'] += r.get('D_hits', 0)
    out['calibration'] = C
    cand = T['cand'] + C['cand']
    out['total'] = dict(**T, core_hours=round(T['secs'] / 3600, 2), cand_all=cand,
                        cand_log2=round(math.log2(cand), 3) if cand else None, cdgs=cdgs,
                        expected_cores=round(cand * 2.0 ** -cdgs, 3),
                        cores_found=len(glob.glob(os.path.join(ROOT, "logs", "p647", "CORE_*.json"))))
    if as_json: print(json.dumps(out))
    else:
        for w, W in out['workers'].items(): print(w, json.dumps(W))
        print('calibration', json.dumps(C)); print('TOTAL', json.dumps(out['total']))
    return out

if __name__ == "__main__":
    a = sys.argv[1:]
    main(float(a[a.index("--cdgs") + 1]) if "--cdgs" in a else 34.6, "--json" in a)
