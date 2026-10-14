r"""Shared helpers for the 647 x353 cross level (theory/p647_*.py).

load_setup()       -> (X, ys): the Cross object of nd_cross.py for 646 353 17 19 and the 495 y_B orbit representatives
                      (read from the restart-safe enumeration log logs/nd/joint/enum_646_353_17_19.jsonl; no recompute).
joint_solutions(k) -> the joint solutions (y_A per block as layer pairs at the Z_19 orbit reps) of y_B #k, cached in
                      logs/p647/joint_yb{k}.json after the first exact JointTest.test (140 s for k = 0)."""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nd_cross import Cross
from nd_joint_exhaust import enumerate_yB

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOG = os.path.join(ROOT, "logs", "p647")
M_, MU, QA, QB = 646, 353, 17, 19

def load_setup():
    X = Cross(M_, MU, QA, QB)
    ys = enumerate_yB(M_, MU, QA, QB, os.path.join(ROOT, "logs", "nd", "joint", f"enum_{M_}_{MU}_{QA}_{QB}.jsonl"),
                      1, log=lambda s: None)
    return X, ys

def joint_solutions(X, ys, k=0, log=print):
    path = os.path.join(LOG, f"joint_yb{k}.json")
    if os.path.exists(path):
        d = json.load(open(path)); return d['rec'], d['sols']
    t0 = time.time()
    rec, sols = X.J.test(ys[k])
    os.makedirs(LOG, exist_ok=True)
    tmp = path + ".tmp"
    json.dump(dict(k=k, yb=[list(map(int, y)) for y in ys[k]], rec=rec, sols=sols), open(tmp, "w"))
    os.replace(tmp, path)
    log(json.dumps(dict(stage="joint_cached", k=k, joint=len(sols), secs=round(time.time() - t0, 1))))
    return rec, sols
