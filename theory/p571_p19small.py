r"""Exact Paley(19) b = 1 pipeline (theory/paley_b1pipe.run) at the small p = 19 orders v = 19 n + 1, as controls for
the rigidity of stage 1 in the <61> = Paley(19) slice of 571 (report section 7).  Logs go to
theory/logs/paley_b1pipe_v{v}_p19.log (tagged, so no repo log is overwritten).
usage: python3 theory/p571_p19small.py v n [max_cores]"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paley_b1pipe as P

if __name__ == "__main__":
    v, n = int(sys.argv[1]), int(sys.argv[2]); mc = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    res = P.run(v, 19, n, mc, tag='_p19')
    print(json.dumps(dict(v=v, n=n, cores=len(res), results=[{k: r[k] for k in r if k != 'X'} for r in res])), flush=True)
