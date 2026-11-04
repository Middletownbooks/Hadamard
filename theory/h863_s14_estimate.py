r"""[E, M sizes] The S14 sibling at v = 863 (census key fbf4f510b2566552; the S7 host plus tau = (q1, q2, u, t) ->
(q2, q1, -u, t); r = 77 orbits per block).  Set-up of both hosts by OpenAI Codex (theory/h3452_compress.py, kind s14);
construction Kumar / EKLO.  No y is known for this host (Codex: two Z43-fibres reported INFEASIBLE by CP-SAT, not
independently confirmed), so a synthetic tau-invariant y with the H-level lag-0 energy is used and the Gaussian is read
at the mean.  Reports, given y: log2 fibre, C(all | y) and rank, log2 E[cores | y], and the 523-style Z2 split along
the tau-fixed involution z = (1, 1): C_S, h_D (exact, per pair of tau-orbits), D_def = C(all) - C_S - h_D.
usage: python3 theory/h863_s14_estimate.py [N] [reps]"""
import sys, os, json, math, itertools as it
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h863_lib as L
from h863_analog import gauss_bits_mean
from h863_split_deficits import row_pair_entropy

SW = [0, 2, 1, 3]                                       # swap of the two Z2 bits, q = 2 q1 + q2
TAU = [SW[k // 5] * 5 + (-(k % 5)) % 5 for k in range(20)]
ORB = sorted({tuple(sorted({k, TAU[k]})) for k in range(20)})
assert len(ORB) == 11
ZM = 3                                                  # z = (1, 1): q -> q ^ 3 commutes with the swap


def summap(members, key):
    rows = {}
    for q in range(4):
        for u in range(5):
            for t in range(43):
                k = key(q, u, t)
                if k in rows: continue
                v = np.zeros(72, np.int64)
                for g in members(q, u, t): v[L.LAGCLS[g]] += 1
                rows[k] = v
    return np.unique(np.array(list(rows.values())), axis=0)


M_S = summap(lambda q, u, t: [(q, u, t), (q ^ ZM, u, t)], lambda q, u, t: (min(q, q ^ ZM), u, t))


def synthetic_y(rng):
    """per block and tau-orbit: (e, S) with S in {0, +-2}; cells with S != 0 have e S = -2 (y^2 - 1 = 168), and their
    number (counting both cells of a paired orbit) is 17, so sum y^2 = 80 + 17 * 168 = 2936"""
    while True:
        E = np.ones((4, 20), np.int64); S = np.zeros((4, 20), np.int64)
        for i in range(4):
            for o in ORB:
                e = int(rng.choice([-1, 1]))
                for k in o: E[i, k] = e
        cnt = 0
        for idx in rng.permutation(4 * 11):
            i, oi = divmod(int(idx), 11); o = ORB[oi]
            if cnt + len(o) > 17: continue
            s = 2 * int(rng.choice([-1, 1]))
            for k in o: S[i, k] = s; E[i, k] = -np.sign(s)
            cnt += len(o)
            if cnt == 17: return E, S


def sample(E, S, N, rng):
    c = np.empty((N, 4, 20, 6), np.int8)
    for i in range(4):
        for o in ORB:
            n = (int(S[i, o[0]]) + 6) // 2
            sub = L._SUBS[n]
            r = 2 * sub[rng.integers(len(sub), size=N)] - 1
            for k in o: c[:, i, k] = r
    return c


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 12000
    reps = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    rng = np.random.default_rng(14)
    out = []
    for rep in range(reps):
        E, S = synthetic_y(rng)
        x = L.assemble(np.broadcast_to(E, (N, 4, 20)), sample(E, S, N, rng))
        # tau-invariance of the samples
        xs = x[0].reshape(4, 20, 43)
        assert all((xs[:, k] == xs[:, TAU[k]]).all() for k in range(20))
        V = L.class_paf(x).sum(1)
        ga = gauss_bits_mean(V); gs = gauss_bits_mean(V @ M_S.T)
        lf = sum(math.log2(math.comb(6, (int(S[i, o[0]]) + 6) // 2)) for i in range(4) for o in ORB)
        # h_D: unordered pairs of distinct orbits {O, O + z}; an orbit with O + z = O carries no D freedom
        def shift(o): return tuple(sorted({((k // 5) ^ ZM) * 5 + k % 5 for k in o}))
        seen = set(); hD = 0.0; stable = 0
        for i in range(4):
            for o in ORB:
                o2 = shift(o)
                assert o2 in ORB
                if o2 == o: stable += 1; continue
                key = (i, min(o, o2), max(o, o2))
                if key in seen: continue
                seen.add(key)
                hD += row_pair_entropy(int(S[i, o[0]]), int(S[i, o2[0]]))
        d = dict(host='863 S14 fbf4f510b2566552', rep=rep, N=N, read_at='mean', log2_fibre_y=round(lf, 2),
                 C_all=ga['bits'], rank_all=ga['rank'], log2E_cores_per_y=round(lf - ga['bits'], 2),
                 Z2_z='(1,1)', C_S=gs['bits'], rank_S=gs['rank'], h_D=round(hD, 2), stable_orbits=stable,
                 C_D_given_S=round(ga['bits'] - gs['bits'], 2), D_def=round(ga['bits'] - gs['bits'] - hD, 2))
        print(json.dumps(d), flush=True); out.append(d)
    with open('logs/h863/s14_estimate.jsonl', 'a') as fo:
        for d in out: fo.write(json.dumps(d) + '\n')


if __name__ == '__main__':
    main()
