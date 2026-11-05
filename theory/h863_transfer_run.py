r"""[E] Driver: the lift-chain numbers of theory/h863_analog.py at the analogues (519, 583) and the transfer orders of
report/prompt_863_lift.md section 2 (H = Z2^2 x Z5 or Z2^2 x Z7 with a large prime Z_p), for multipliers of odd order f
(so -1 is not in S and the classes pair under negation, as at the 863 S7 host).  Synthetic y (H-level lag-0 energy,
random positions), two reps each; read at the mean.  Output: logs/h863/analog_transfer.jsonl.
usage: python3 theory/h863_transfer_run.py [N]"""
import sys, os, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h863_analog as A


def mult_of_order(p, f):
    for a in range(2, p):
        x = a; k = 1
        while x != 1:
            x = x * a % p; k += 1
        if k == f: return a
    return None


HOSTS = [  # (label, m, p, f)
    ('analogue 519', 3, 43, 7), ('analogue 583', 5, 29, 7),
    ('863 S7', 5, 43, 7),
    ('1223', 5, 61, 15), ('1223', 5, 61, 5),
    ('1423', 5, 71, 7), ('1583', 5, 79, 13), ('1783', 5, 89, 11), ('2063', 5, 103, 17),
    ('1039', 7, 37, 9), ('1039', 7, 37, 3), ('1319', 7, 47, 23), ('1487', 7, 53, 13), ('1879', 7, 67, 11),
    ('2215', 7, 79, 13),
]

if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
    out = 'logs/h863/analog_transfer.jsonl'
    for label, m, p, f in HOSTS:
        a = mult_of_order(p, f)
        host = A.Host(m, p, a)
        assert host.f == f and not host.neg_in_S
        rng = np.random.default_rng(1000 + p * 7 + m)
        for rep in range(2):
            try:
                E, S, info = host.synthetic_y(rng)
            except ValueError as ex:
                print(json.dumps(dict(label=label, m=m, p=p, f=f, error=str(ex)))); break
            res = A.chain_numbers(host, E, S, N, rng)
            res.update(info); res.update(label=label, N=N, rep=rep, read_at='mean')
            s = json.dumps(res); print(s, flush=True)
            with open(out, 'a') as fo: fo.write(s + '\n')
