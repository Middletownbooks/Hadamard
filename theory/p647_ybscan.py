r"""Joint stage of the other y_B orbits of 647 x353 (#1 .. #494), capped and restart-safe [M].

For every y_B orbit (theory/p647_common.load_setup; 495 orbits, #0 is done: 20,928 joint solutions, 436 classes) the
exact joint test nd_joint_exhaust.JointTest.test is run with the pair lists capped at --maxpairs per block (the
uncapped test of y_B #1 needed more than the 15 GB of this machine).  A y_B over the cap, or one whose run is killed
inside the test twice (an 'attempt' line is written first; one retry), is recorded as undecided.  For a y_B with joint solutions:
  * the joint solutions are saved to logs/p647/joint_yb{k}.json (the format of p647_common.joint_solutions);
  * their classes under the full stabiliser of that y_B, per-block reversal included (theory/p647_classes.py), go to
    logs/p647/joint_yb{k}_classes.json, with the H (= sum of log2 label counts) of every class representative.
A class's full D-candidate yield is about 2^(H - 50.3) (calibration, report section 2.2), so the H histogram ranks the
y_B orbits for production.
Output: one line per y_B in logs/p647/ybscan.jsonl.
The joint test runs in a child under an address-space cap (--mem GB, theory/p647_captest.py).
usage: python3 theory/p647_ybscan.py [--from K] [--to K] [--ks k1,k2,...] [--maxpairs N] [--mem GB]"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p647_common as PC
import p647_classes as PCL
from p647_captest import capped_test

def main(k0=1, k1=None, maxpairs=4 * 10 ** 6, log=print, mem_gb=6.0, ks=None):
    X, ys = PC.load_setup(); k1 = len(ys) - 1 if k1 is None else k1
    ks = list(range(k0, k1 + 1)) if ks is None else ks
    outp = os.path.join(PC.LOG, "ybscan.jsonl")
    done = set(); attempted = {}
    if os.path.exists(outp):
        for line in open(outp):
            try: r = json.loads(line)
            except Exception: continue
            if 'joint' in r: done.add(r['k'])
            if 'attempt' in r: attempted[r['k']] = attempted.get(r['k'], 0) + 1
    log(json.dumps(dict(stage="start", yB_orbits=len(ys), range=[k0, k1], done=len(done), maxpairs=maxpairs)))
    with open(outp, "a") as fo:
        for k in ks:
            if k in done: continue
            if attempted.get(k, 0) >= 2:                      # killed twice inside the test: give up
                fo.write(json.dumps(dict(k=k, joint=None, note="run killed inside the joint test (memory); undecided")) + "\n")
                fo.flush(); continue
            fo.write(json.dumps(dict(k=k, attempt=1)) + "\n"); fo.flush()
            t0 = time.time()
            status, rec, sols = capped_test(X, ys[k], maxsol=10 ** 7, maxpairs=maxpairs, mem_gb=mem_gb)
            out = dict(k=k, joint=rec.get('joint'), S_solutions=rec.get('S_solutions'), pairs=rec.get('pairs'),
                       layer_options=rec.get('layer_options'), note=rec.get('note'), pairs_needed=rec.get('pairs_needed'),
                       secs=round(time.time() - t0, 1))
            if sols:
                path = os.path.join(PC.LOG, f"joint_yb{k}.json")
                json.dump(dict(k=k, yb=[list(map(int, y)) for y in ys[k]], rec=rec, sols=sols), open(path + ".tmp", "w"))
                os.replace(path + ".tmp", path)
                stab = PCL.stabiliser_full(ys[k])
                cls, reps = PCL.classes(X, sols, stab)
                H = [sols[r]['mult_log2'] for r in reps]
                json.dump(dict(k=k, classes=len(reps), stabiliser_full_on_Z38=len(stab), class_of=cls, reps=reps, H=H),
                          open(os.path.join(PC.LOG, f"joint_yb{k}_classes.json"), "w"))
                hs = np.array(H)
                out.update(stabiliser=len(stab), classes=len(reps), H_max=float(hs.max()), H_median=float(np.median(hs)),
                           classes_H_ge_79=int((hs >= 79.0).sum()),
                           yield_sum_log2=round(float(np.log2(np.sum(2.0 ** (hs - 50.3)))), 2))
            fo.write(json.dumps(out) + "\n"); fo.flush()
            log(json.dumps(out))

if __name__ == "__main__":
    a = sys.argv[1:]
    g = lambda f, d: int(a[a.index(f) + 1]) if f in a else d
    main(g("--from", 1), g("--to", None), g("--maxpairs", 4 * 10 ** 6), log=lambda s: print(s, flush=True),
         mem_gb=float(a[a.index("--mem") + 1]) if "--mem" in a else 6.0,
         ks=[int(x) for x in a[a.index("--ks") + 1].split(",")] if "--ks" in a else None)
