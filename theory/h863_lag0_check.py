r"""[P, M] H-level lag-0 test for the k4 hosts of the transfer table (theory/h863_transfer_run.py): a y exists only if
sum over the 4 * 4m cells of y^2 = 4|G| - 12(p - 1), with y = e + f S, e = +-1, S = sum of e_cls class signs
(S = e_cls mod 2, |S| <= e_cls).  Exact subset-sum over the multiset of possible y^2 - 1 values (unbounded counts,
at most 16m cells).  usage: python3 theory/h863_lag0_check.py"""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h863_transfer_run import HOSTS
for label, m, p, f in HOSTS:
    e = (p - 1) // f; G = 4 * m * p; cells = 16 * m
    vals = sorted({(1 + f * S) ** 2 - 1 for S in range(-e, e + 1, 2)} | {(-1 + f * S) ** 2 - 1 for S in range(-e, e + 1, 2)})
    tgt = 4 * G - 12 * (p - 1) - cells
    # reach[k] = set of sums using exactly k nonzero cells (k <= cells)
    reach = {0}
    ok = tgt == 0
    for k in range(cells):
        reach = {s + v for s in reach for v in vals if v > 0 and s + v <= tgt}
        if tgt in reach: ok = True; break
        if not reach: break
    print(json.dumps(dict(label=label, m=m, p=p, f=f, e=e, lag0_target_minus_cells=tgt, y2m1_values=vals[:6], solvable=ok)))
