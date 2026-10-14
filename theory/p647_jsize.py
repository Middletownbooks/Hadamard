r"""Size of the joint test of every y_B orbit, from its cheap stages only [M].

Replicates the first stages of nd_joint_exhaust.JointTest.test (per-layer options, the exact energy filters, the
energy-filtered pair lists, distinct sigma per block) without the two 4-list joins, so it runs in seconds and bounded
memory.  The join sizes it reports (sigma_distinct per block, and the 2+2 pair counts of the S-join, prod over the two
sides) rank the y_B orbits, so that a capped run can take them smallest first.
Each test runs under the address-space cap of theory/p647_captest.py (a child per y_B).
Output: one line per y_B in logs/p647/<out>.jsonl (k, layer_options, pairs, sigma_distinct, join_log2, slack).
usage: python3 theory/p647_jsize.py m mu qA qB OUT [--mem GB]"""
import sys, os, json, math, time
import multiprocessing as mp, resource, traceback
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nd_cross import Cross
from nd_joint_exhaust import enumerate_yB

def sizes(J, yb, maxpairs=4 * 10 ** 7):
    TS, TD = J.TS, J.TD
    Ls = [J.block_options(yb[k]) for k in range(4)]
    rec = dict(layer_options=[(len(a[0]), len(b[0])) for a, b in Ls])
    if min(min(s_) for s_ in rec['layer_options']) == 0: return dict(rec, empty="layer options")
    En = [[(V * V * J.nw).sum(1) for V, c in L] for L in Ls]
    slack = 2 * int(TS[0] + TD[0]) - sum(int(e.min()) for E in En for e in E)
    rec['energy_slack'] = slack
    if slack < 0: return dict(rec, empty="energy")
    P = []
    for k, ((V0, c0), (V1, c1)) in enumerate(Ls):
        E0, E1 = En[k]; bound = int(E0.min()) + int(E1.min()) + slack
        o1 = np.argsort(E1, kind='stable'); E1s = E1[o1]
        cntj = np.searchsorted(E1s, bound - E0, 'right')
        if int(cntj.sum()) > maxpairs: return dict(rec, pairs_block=k, pairs_needed=int(cntj.sum()), note="pairs over the size cap")
        I0 = np.repeat(np.arange(len(V0)), cntj)
        I1 = o1[np.concatenate([np.arange(c) for c in cntj])] if cntj.sum() else np.zeros(0, np.int64)
        sg = (V0[I0] + V1[I1]) // 2; dl = (V0[I0] - V1[I1]) // 2
        P.append(dict(sg=sg, se=(sg * sg * J.nw).sum(1), de=(dl * dl * J.nw).sum(1)))
        del I0, I1, dl
    for it in range(6):
        ms = [int(p_['se'].min()) if len(p_['se']) else 10 ** 9 for p_ in P]
        md = [int(p_['de'].min()) if len(p_['de']) else 10 ** 9 for p_ in P]
        changed = False
        for k, p_ in enumerate(P):
            ok = (p_['se'] <= TS[0] - (sum(ms) - ms[k])) & (p_['de'] <= TD[0] - (sum(md) - md[k]))
            if not ok.all():
                changed = True
                for key in p_: p_[key] = p_[key][ok]
        if any(len(p_['se']) == 0 for p_ in P): return dict(rec, pairs=[len(p_['se']) for p_ in P], empty="energy pruning")
        if not changed: break
    rec['pairs'] = [int(len(p_['se'])) for p_ in P]
    nd = [int(len(np.unique(p_['sg'], axis=0))) for p_ in P]
    rec['sigma_distinct'] = nd
    rec['join_log2'] = round(math.log2(nd[0] * nd[1] + nd[2] * nd[3]), 2)
    return rec

def _child(conn, J, yb, mem):
    try:
        resource.setrlimit(resource.RLIMIT_AS, (mem, mem)); conn.send(('ok', sizes(J, yb)))
    except MemoryError: conn.send(('memory', None))
    except BaseException: conn.send(('error', traceback.format_exc()[-300:]))
    finally: conn.close()

def main(m, mu, qA, qB, out, mem_gb=5.0, log=print):
    X = Cross(m, mu, qA, qB)
    if (m, mu) == (646, 353):
        enum = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "logs", "nd", "joint", f"enum_{m}_{mu}_{qA}_{qB}.jsonl")
    else:
        enum = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "logs", "p647", "control", f"enum_{m}_{mu}_{qA}_{qB}.jsonl")
    ys = enumerate_yB(m, mu, qA, qB, enum, 1, log=lambda s: None)
    outp = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "logs", "p647", out + ".jsonl")
    done = set()
    if os.path.exists(outp):
        for line in open(outp):
            try: done.add(json.loads(line)['k'])
            except Exception: pass
    ctx = mp.get_context('fork')
    with open(outp, "a") as fo:
        for k in range(len(ys)):
            if k in done: continue
            t0 = time.time()
            rd, wr = ctx.Pipe(duplex=False)
            p = ctx.Process(target=_child, args=(wr, X.J, ys[k], int(mem_gb * 2 ** 30))); p.start(); wr.close()
            try: st, rec = rd.recv()
            except (EOFError, OSError): st, rec = 'died', None
            p.join()
            r = dict(k=k, status=st, secs=round(time.time() - t0, 2))
            if st == 'ok': r.update(rec)
            elif st == 'error': r['error'] = rec
            fo.write(json.dumps(r) + "\n"); fo.flush()
            log(json.dumps(r))

if __name__ == "__main__":
    a = sys.argv[1:]
    main(*(int(x) for x in a[:4]), a[4], float(a[a.index("--mem") + 1]) if "--mem" in a else 5.0, log=lambda s: print(s, flush=True))
