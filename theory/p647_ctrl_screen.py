r"""Exact stage-1 screen of the product-shape cyclic slices (theory/p647_control.py --list) [M].
For each slice and each of its two compressions (along Z_p onto Z_2q, along Z_q onto Z_2p) the exact compressed
stage 1 (theory/nd_s1_sd.py: alpha-congruence + S/D split, exhaustive) is asked for one solution.  Both compressions
are necessary conditions, so one empty compression proves the slice empty.
usage: python3 theory/p647_ctrl_screen.py VMAX [--secs 600] > logs/p647/control/screen.jsonl"""
import sys, os, json, time, signal
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nd_s1_sd import Problem
from p647_control import product_slices

class TO(Exception): pass
def handler(s, f): raise TO()

if __name__ == "__main__":
    vmax = int(sys.argv[1]); secs = int(sys.argv[sys.argv.index("--secs") + 1]) if "--secs" in sys.argv else 600
    signal.signal(signal.SIGALRM, handler)
    for (v, p, q, dp, dq, mu) in product_slices(vmax):
        m = v - 1; res = {}
        for K in (p, q):                      # compress along Z_K
            t0 = time.time(); signal.alarm(secs)
            try:
                P = Problem(m, mu, K)
                out, sols = P.solve(max_sol=1, log=lambda s: None)
                res[K] = dict(sol=len(sols), complete=out['complete'], quads=out['alpha_quads'], S=out['S_solutions'],
                              f0=P.f0, secs=round(time.time() - t0, 1))
            except TO:
                res[K] = dict(sol=None, note=f"timeout {secs}s")
            except Exception as ex:
                res[K] = dict(sol=None, note=str(ex)[:100])
            finally:
                signal.alarm(0)
        empty = any(r.get('sol') == 0 for r in res.values())
        print(json.dumps(dict(v=v, p=p, q=q, Mp=dp, Mq=dq, mu=mu, along_p=res[p], along_q=res[q],
                              verdict="EMPTY" if empty else ("alive" if all(r.get('sol') for r in res.values()) else "undecided"))), flush=True)
