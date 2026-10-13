"""Run the exact lag-0 + rows + Z_2 screen (theory/gp_z2screen.py) on every generalised-Paley slice of the open orders
with a positive slice Gaussian (logs/msweep_record, logs/msweep34; pure slices only: order of the multiplier group =
(q - 1)/s).  Output: one JSON line per slice.  usage: python3 theory/gp_screen_candidates.py [timeout_secs]"""
import sys, os, json, signal
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gp_z2screen import screen

# (v, q, s, slice log2E) from the multiplier sweeps
L = [(523, 29, 4, 92.1), (907, 151, 10, 66.6), (955, 53, 4, 61.8), (537, 67, 6, 59.1), (739, 41, 4, 57.3),
     (917, 229, 12, 54.1), (823, 137, 8, 52.4), (971, 97, 8, 46.1), (907, 151, 6, 32.1), (947, 43, 3, 30.8),
     (971, 97, 6, 28.4), (599, 23, 2, 27.3), (523, 29, 2, 11.5), (789, 197, 4, 9.3), (971, 97, 4, 9.3),
     (917, 229, 4, 8.7), (659, 47, 2, 2.9), (947, 43, 2, 1.7), (739, 41, 2, 1.5)]

class TO(Exception): pass
def _alarm(*a): raise TO()

if __name__ == "__main__":
    tmax = int(sys.argv[1]) if len(sys.argv) > 1 else 900
    signal.signal(signal.SIGALRM, _alarm)
    for v, q, s, g in L:
        F = (v - 1) // q; signal.alarm(tmax)
        try:
            d = screen(F, q, s); signal.alarm(0)
            surv = [r['comp'] for r in d['detail'] if r['rows_and_Z2_ok']]
            print(json.dumps(dict(v=v, F=F, q=q, s=s, f=(q - 1) // s, log2E=g, comps=d['compositions'],
                                  survivors=len(surv), rows_only=sum(r['rows_ok'] for r in d['detail']),
                                  surv=surv[:6])), flush=True)
        except TO:
            print(json.dumps(dict(v=v, F=F, q=q, s=s, log2E=g, timeout=True)), flush=True)
