r"""Gate A1/A2: the Z_13 level (trivial character of K) of the plug-in slice on Z_n x F_13, decided over all stage-1
profiles [P, M].  Necessary conditions for a core, all exact [P, theory/h2396_plugin_system.py]:
  per block i: D_i defects (sum D_i = e_plus + e_minus for an allocation of the defect identity), defect sum
  e^_i = sum_k e_i(k) with |e^_i| <= D_i and e^_i = D_i (mod 2); row sum r_i = sum_k u_i(k) = a_i + 12 e^_i with
  a_i = sum_k x0_i(k) in [-n, n], a_i = n (mod 2); sum_i r_i^2 = 4 (n even: one r_i = +-2, the rest 0; n odd: all +-1);
  t^_i = sum_k t_i(k) over the n - D_i clean positions: |t^_i| <= n - D_i, t^_i = n - D_i (mod 2);
  with s^_i = r_i - 13 e^_i:   sum_i (s^_i^2 + 13 t^_i^2) = 52 n + 4   and   sum_i s^_i t^_i = 0.
The order-2 character (n even) gives the same equations with psi-sums but adds new unknowns; it is not used here.
The script lists, for every allocation and multiset of D_i, which (e^, r) profiles admit some t^.  If none does, the
slice has no core [P + M].  Also reports, for checking, the profiles of the sampled stage-1 solutions.
usage: python3 theory/h2396_zq_profiles.py n [--show]"""
import sys, itertools, json
from functools import lru_cache


def allocations(n, q=13):
    R = 4 * (n - 1)
    return [(a, (R - (q + 1) * a) // (q - 3)) for a in range(R // (q + 1) + 1) if (R - (q + 1) * a) % (q - 3) == 0]


def that_ok(n, D, eh, r):
    """exists t^ with the Z_13-level equations for the profile (D_i, e^_i, r_i)"""
    sh = [ri - 13 * ei for ri, ei in zip(r, eh)]
    rest = 52 * n + 4 - sum(x * x for x in sh)
    if rest < 0 or rest % 13: return None
    R2 = rest // 13
    rng = [[t for t in range(-(n - d), n - d + 1) if (t - (n - d)) % 2 == 0 and t * t <= R2] for d in D]
    for t0 in rng[0]:
        for t1 in rng[1]:
            q1 = R2 - t0 * t0 - t1 * t1
            if q1 < 0: continue
            for t2 in rng[2]:
                q2 = q1 - t2 * t2
                if q2 < 0: continue
                for t3 in rng[3]:
                    if t3 * t3 == q2 and sh[0] * t0 + sh[1] * t1 + sh[2] * t2 + sh[3] * t3 == 0:
                        return (t0, t1, t2, t3)
    return None


def profiles(n):
    out = []
    rowsets = []
    if n % 2 == 0:
        for big in range(4):
            for sg in (2, -2): rowsets.append(tuple(sg if b == big else 0 for b in range(4)))
    else:
        rowsets = list(itertools.product((1, -1), repeat=4))
    for (ep, em) in allocations(n):
        E = ep + em
        for D in itertools.product(range(E + 1), repeat=4):
            if sum(D) != E: continue
            if max(D) > n: continue
            for eh in itertools.product(*[range(-d, d + 1, 2) for d in D]):
                for r in rowsets:
                    if any(abs(ri - 12 * ei) > n or (ri - 12 * ei - n) % 2 for ri, ei in zip(r, eh)): continue
                    w = that_ok(n, D, eh, r)
                    out.append(dict(alloc=(ep, em), D=D, eh=eh, r=r, ok=w is not None, witness=w))
    return out


if __name__ == "__main__":
    n = int(sys.argv[1])
    P = profiles(n)
    ok = [p for p in P if p['ok']]
    by_alloc = {}
    for p in P:
        a = str(p['alloc']); by_alloc.setdefault(a, [0, 0]); by_alloc[a][0] += 1; by_alloc[a][1] += p['ok']
    odd_counts = sorted({sum(d % 2 for d in p['D']) for p in ok})
    print(json.dumps(dict(n=n, profiles=len(P), alive=len(ok), by_alloc=by_alloc, odd_blocks_in_alive=odd_counts)))
    if '--show' in sys.argv:
        for p in ok[:40]: print(json.dumps(p))
