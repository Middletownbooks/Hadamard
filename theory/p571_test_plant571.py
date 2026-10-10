r"""Planted test of the S-level join at full v = 571 size (x27 slice, a real stage-1 solution y): random S-class
sublists of 2^14 per block, a random planted quadruple defines the target, and theory/p571_hsolve.join_all must return
it.  Checks that the low S-solution counts of theory/p571_s3calib.py are not a join failure.
usage: python3 theory/p571_test_plant571.py [trials]"""
import sys, os, time, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p571_hsolve import Inst, Stage2Block, _sig_eis1, join_all

if __name__ == "__main__":
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    inst = Inst(570, [103]); root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    y = np.load(os.path.join(root, 'theory/logs/p571/y571_zl0.npz'))['Y'][0]
    rng = np.random.default_rng(7); sub = 1 << 14
    blocks = [Stage2Block(inst, y[k]) for k in range(4)]
    sigs = []
    for b in blocks:
        A, B, _ = b.classes(rng.choice(b.n_classes, size=sub, replace=False).astype(np.int64)); sigs.append(_sig_eis1(A, B, inst.N))
    ok = 0
    for t in range(trials):
        pl = [int(rng.integers(sub)) for _ in range(4)]
        T = sum(sigs[k][pl[k]] for k in range(4))
        t0 = time.time(); hits, _ = join_all(sigs, T, seed=11 + t, tbits=20)
        found = any(list(h) == pl for h in hits); ok += found
        print(json.dumps(dict(trial=t, hits=len(hits), planted_found=found, secs=round(time.time() - t0, 1))), flush=True)
    print(json.dumps(dict(summary=True, found=ok, trials=trials)), flush=True)
    assert ok == trials
