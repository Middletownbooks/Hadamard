r"""Class-deduplicated tally of a full control run (theory/p647_ctrl_full.py) [M, with E alongside].

Reads logs/p647/control/full_<tag>.jsonl (one line per processed joint solution) and the class maps
logs/p647/control/classes_<tag>.json written by a --classes run (the full stabiliser of y_B, per-block reversal
included; theory/p647_classes.py).  Joint solutions of one class have isomorphic cross levels, so
  * it checks that every processed member of a class has the same exact counts (S-solutions, D-candidates, D-hits);
  * it counts every class once: distinct D-candidates, and the expected number of distinct cores
    sum_class D-candidates * 2^-C(D|S), C(D|S) from the class's own Gaussian fit where one was made (every 10th
    processed joint solution), otherwise the median of the fits.
usage: python3 theory/p647_ctrl_tally.py [m mu qA qB]   (default 494 83 13 19)"""
import sys, os, json, math, statistics
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CTRL = os.path.join(ROOT, "logs", "p647", "control")

def main(m=494, mu=83, qA=13, qB=19, as_json=False):
    tag = f"{m}_{mu}_{qA}_{qB}"
    recs = [json.loads(l) for l in open(os.path.join(CTRL, f"full_{tag}.jsonl"))]
    cm = json.load(open(os.path.join(CTRL, f"classes_{tag}.json"))) if os.path.exists(os.path.join(CTRL, f"classes_{tag}.json")) else {}
    tests = {}
    for r in recs:
        if 'joint_total' in r: tests[r['yb']] = r
    js = [r for r in recs if 'jsol' in r]
    fits = [r['CDgS'] for r in js if r.get('CDgS') is not None and math.isfinite(r['CDgS'])]
    cmed = statistics.median(fits) if fits else None
    out = dict(tag=tag, yB_tested=len(tests), yB_with_joint=sum(1 for t in tests.values() if t.get('joint_total')),
               yB_undecided=sorted(i for i, t in tests.items() if t.get('joint_total') is None),
               CDgS_median=round(cmed, 2) if cmed else None, per_yb={})
    T = dict(classes=0, classes_done=0, cand=0, expected=0.0, hits=0, cores=0, mismatches=0, members_checked=0)
    for ib, t in sorted(tests.items()):
        if not t.get('joint_total'): continue
        mp = cm.get(str(ib))
        rows = [r for r in js if r['yb'] == ib]
        if mp is None:
            out['per_yb'][ib] = dict(joint=t['joint_total'], old_orbits_done=len(rows), classes=None, note="no class map yet")
            continue
        cof = mp['class_of']; ncls = len(mp['reps'])
        per = {}
        for r in rows: per.setdefault(cof[r['jsol']], []).append(r)
        mism = 0
        for c, rr in per.items():
            k0 = (rr[0]['S'], rr[0]['Spairs'], rr[0]['D_hits'])
            mism += sum(1 for r in rr[1:] if (r['S'], r['Spairs'], r['D_hits']) != k0)
            T['members_checked'] += len(rr)
        cand = sum(rr[0]['Spairs'] for rr in per.values())
        exp_ = 0.0
        for rr in per.values():
            cf = next((r['CDgS'] for r in rr if r.get('CDgS') is not None and math.isfinite(r['CDgS'])), cmed)
            if cf is not None: exp_ += rr[0]['Spairs'] * 2.0 ** -cf
        hits = sum(rr[0]['D_hits'] for rr in per.values()); cores = sum(rr[0]['cores'] for rr in per.values())
        out['per_yb'][ib] = dict(joint=t['joint_total'], classes=ncls, classes_done=len(per), old_orbits_done=len(rows),
                                 distinct_cand=cand, cand_log2=round(math.log2(cand), 2) if cand else None,
                                 expected_distinct_cores=round(exp_, 3), D_hits=hits, cores=cores, member_mismatches=mism)
        T['classes'] += ncls; T['classes_done'] += len(per); T['cand'] += cand; T['expected'] += exp_
        T['hits'] += hits; T['cores'] += cores; T['mismatches'] += mism
    T['cand_log2'] = round(math.log2(T['cand']), 2) if T['cand'] else None
    T['expected'] = round(T['expected'], 3)
    T['poisson_p0'] = round(math.exp(-T['expected']), 3)
    out['total'] = T
    if as_json: print(json.dumps(out))
    else:
        for ib, v in out['per_yb'].items(): print(ib, json.dumps(v))
        print(json.dumps({k: v for k, v in out.items() if k != 'per_yb'}))
    return out

if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if not x.startswith("--")]
    main(*(int(x) for x in a[:4]) if len(a) >= 4 else (), as_json="--json" in sys.argv)
