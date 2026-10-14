r"""Status of the 647 x353 production run (logs/p647/prod/) and of the calibration runs (logs/p647/calib/).

Per worker it sums the per-class checkpoint lines (w*_els.jsonl): joint solutions finished and in progress, emitted
pairs, S-solutions, D-candidates (= pair-level S-solutions, prod_k mult(sigma_k) summed over S-solutions) and D-hits;
and over everything the expected number of cores  sum D-candidates * 2^-C(D|S)  with C(D|S) given (calibrated value
from report/paley_647_answer.md section 2.3; default 34.8 bits), the D-candidates counted once per true symmetry class
and energy class (theory/p647_classes.py; members of a class are copies).  Cores found: logs/p647/CORE_*.json.
usage: python3 theory/p647_tally.py [--cdgs 34.6] [--json]"""
import sys, os, json, glob, math
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROD = os.path.join(ROOT, "logs", "p647", "prod"); CAL = os.path.join(ROOT, "logs", "p647", "calib")

def main(cdgs=34.8, as_json=False):
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
    # unit mode (class-merged sub-bicliques): per (unit, i, j) joint solution, D-candidates over the processed classes
    U = dict(units_done=0, units_partial=0, pairs=0, S=0, cand=0, hits=0, secs=0.0)
    per_joint = {}
    for f in sorted(glob.glob(os.path.join(PROD, "u*_els.jsonl"))):
        for line in open(f):
            try: r = json.loads(line)
            except Exception: continue
            U['pairs'] += r['left'] + r['right']; U['S'] += r['S']; U['cand'] += r['cand']; U['hits'] += r['hits']; U['secs'] += r['secs']
            for ij, c in r.get('cand_ij', {}).items():
                e = r['joint'][ij]; per_joint[e] = per_joint.get(e, 0) + c
    ud = set(); up = set()
    for f in sorted(glob.glob(os.path.join(PROD, "u*.jsonl"))):
        if f.endswith("_els.jsonl"): continue
        for line in open(f):
            try: ud.add(json.loads(line)['unit'])
            except Exception: pass
    for f in sorted(glob.glob(os.path.join(PROD, "u*_els.jsonl"))):
        for line in open(f):
            try: up.add(json.loads(line)['unit'])
            except Exception: pass
    U['units_done'] = len(ud); U['units_partial'] = len(up - ud); U['secs'] = round(U['secs'], 1)
    out['units'] = U
    # true symmetry classes (theory/p647_classes.py: the stabiliser of y_B with per-block reversal, 436 classes of 48
    # joint solutions); members of one class have identical counts, so per (class, EL) the largest record counts once
    CJ = json.load(open(os.path.join(ROOT, "logs", "p647", "joint_yb0_classes.json")))
    cof = CJ['class_of']
    per_cel = {}
    def putc(c, EL, v): per_cel[(c, EL)] = max(per_cel.get((c, EL), 0), v)
    # live deep-tail calibration: exact D-candidates hitting nested cross-lag subsets vs the Gaussian of the
    # calibration fit (report section 2.3; bits per nested subset of the 12 cross lags, ranks 1..11).  A unit record
    # sums its members, which are copies (members / classes of them per class), so it is divided by that factor.
    GB = [3.9, 7.74, 11.51, 15.07, 18.68, 21.73, 25.19, 27.65, 30.45, 32.9, 34.63, 34.63]
    nest = [0] * 12; ncand = 0
    for f in sorted(glob.glob(os.path.join(PROD, "u*_els.jsonl"))):
        for line in open(f):
            try: r = json.loads(line)
            except Exception: continue
            for ij, c in r.get('cand_ij', {}).items(): putc(cof[r['joint'][ij]], r['EL'], c)
            if 'nested' in r and r.get('cand_ij'):
                fac = len(r['cand_ij']) // len({cof[e] for e in r['joint'].values()})
                ncand += r['cand'] // fac; nest = [a + b // fac for a, b in zip(nest, r['nested'])]
    for f in sorted(glob.glob(os.path.join(PROD, "w*_els.jsonl"))) + sorted(glob.glob(os.path.join(PROD, "c*_els.jsonl"))):
        for line in open(f):
            try: r = json.loads(line)
            except Exception: continue
            putc(cof[r['j']], r['EL'], r['cand'])
            if 'nested' in r and os.path.basename(f).startswith("c"):
                ncand += r['cand']; nest = [a + b for a, b in zip(nest, r['nested'])]
    out['deep_tail'] = dict(candidates=ncand, candidates_log2=round(math.log2(ncand), 3) if ncand else None,
                            rows=[dict(rank=min(q + 1, 11), hits=nest[q], expected=round(ncand * 2.0 ** -GB[q], 2)) for q in range(3, 12)])
    # class mode
    CM = dict(classes_done=0, classes_partial=0, S=0, cand=0, hits=0, secs=0.0)
    cd = set(); cp_ = set()
    for f in sorted(glob.glob(os.path.join(PROD, "c*.jsonl"))):
        if f.endswith("_els.jsonl"): continue
        for line in open(f):
            try: cd.add(json.loads(line)['class'])
            except Exception: pass
    for f in sorted(glob.glob(os.path.join(PROD, "c*_els.jsonl"))):
        for line in open(f):
            try: r = json.loads(line)
            except Exception: continue
            CM['S'] += r['S']; CM['cand'] += r['cand']; CM['hits'] += r['hits']; CM['secs'] += r['secs']; cp_.add(cof[r['j']])
    CM['classes_done'] = len(cd); CM['classes_partial'] = len(cp_ - cd); CM['secs'] = round(CM['secs'], 1)
    out['class_mode'] = CM
    # calibration runs (unrestricted, complete: every energy class) count as checked classes too
    C = dict(joint_done=0, S=0, cand=0, hits=0)
    cp = os.path.join(CAL, "calib.jsonl")
    calc = {}
    if os.path.exists(cp):
        for line in open(cp):
            r = json.loads(line); C['joint_done'] += 1; C['S'] += r.get('S_solutions', 0); C['cand'] += r.get('S_pairs', 0); C['hits'] += r.get('D_hits', 0)
            calc[cof[r['j']]] = max(calc.get(cof[r['j']], 0), r.get('S_pairs', 0))
    out['calibration'] = C
    per_class = {}
    for (c, EL), v in per_cel.items(): per_class[c] = per_class.get(c, 0) + v
    for c, v in calc.items(): per_class[c] = max(per_class.get(c, 0), v)
    cand = sum(per_class.values())
    # the old dedupe (3,488 orbits of the stabiliser without per-block reversal) for comparison: it counts each class
    # once per member orbit touched, i.e. up to 8 times
    om = json.load(open(os.path.join(ROOT, "logs", "p647", "joint_yb0_orbits.json"))).get('orbit_of')
    best = {}
    for e, c in per_joint.items(): best[om[e]] = max(best.get(om[e], 0), c)
    old = sum(best.values())
    # forecast for y_B #0 at the class mode's K energy classes: per H group (>= 79 / < 79) the mean D-candidates of
    # the classes finished in class mode, times the classes of that group not yet finished there
    Hc = CJ['H']; fin_c = set(); grp = {True: [], False: []}
    for f in sorted(glob.glob(os.path.join(PROD, "c*.jsonl"))):
        if f.endswith("_els.jsonl"): continue
        for line in open(f):
            try: r = json.loads(line)
            except Exception: continue
            fin_c.add(r['class'])
            if not r.get('pre_done'): grp[Hc[r['class']] >= 79.0].append(r['cand'])
    left = {g: sum(1 for c in range(CJ['classes']) if (Hc[c] >= 79.0) == g and c not in fin_c and per_class.get(c, 0) == 0)
            for g in (True, False)}
    fc = {}
    for g, name in ((True, 'H79'), (False, 'H78')):
        mean = sum(grp[g]) / len(grp[g]) if grp[g] else None
        fc[name] = dict(finished_in_class_mode=len(grp[g]), mean_cand_log2=round(math.log2(mean), 3) if mean else None,
                        untouched_classes=left[g], expected_remaining=round(left[g] * mean * 2.0 ** -cdgs, 3) if mean else None)
    out['forecast_yb0'] = fc
    D = out.get('class_mode', {})
    secs = T['secs'] + U['secs'] + D.get('secs', 0.0)
    out['total'] = dict(classes_touched=len(per_class), classes=CJ['classes'], core_hours=round(secs / 3600, 2), cand_classes=cand,
                        cand_log2=round(math.log2(cand), 3) if cand else None, cdgs=cdgs,
                        expected_cores=round(cand * 2.0 ** -cdgs, 3),
                        unit_mode_old_count_expected=round(old * 2.0 ** -cdgs, 3),
                        cores_found=len(glob.glob(os.path.join(ROOT, "logs", "p647", "CORE_*.json"))))
    if as_json: print(json.dumps(out))
    else:
        for w, W in out['workers'].items(): print(w, json.dumps(W))
        print('units', json.dumps(U)); print('class_mode', json.dumps(out['class_mode'])); print('calibration', json.dumps(C))
        print('deep_tail', json.dumps(out['deep_tail'])); print('forecast_yb0', json.dumps(out['forecast_yb0']))
        print('TOTAL', json.dumps(out['total']))
    return out

if __name__ == "__main__":
    a = sys.argv[1:]
    main(float(a[a.index("--cdgs") + 1]) if "--cdgs" in a else 34.8, "--json" in a)
