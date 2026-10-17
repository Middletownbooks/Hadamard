r"""Validation checks of the s1z tools [M].

  canon    exact totals of theory/s1z_lib.ZOrb against the earlier canonical stage-1 counts of the H' = 1 controls
           (v = 35, q = 5: 1,728; v = 47, q = 7, s = 2: 36; logs/ab/paley_c6.jsonl), by enumerating every ordered
           compression quadruple and canonicalising the solutions (sign per block, block order)
  sd       the S/D exact counter (SDCount) against the direct counter on every multiset of 347's slice
  invar    exact counts are equal on random members of a multiset's symmetry orbit (347's slice)
  screens  theory/s1z_qzq.Screen: the compression identity at every level on the published 347 core; planted lifts on
           719-shaped analogs (hashed count = independent exact count, planted y never killed); the v = 47 control
           (count at the y of a known core = 5,184, the cores found there by the exact stage-2 run)
usage: python3 theory/s1z_validate.py canon|sd|invar|screens"""
import sys, os, json, math, itertools, random
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_s1sd as SD, ab_sweep as W
from s1z_lib import ZOrb, SDCount
from s1z_qzq import Screen

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
C6 = lambda v, h, q, s, f, H1, U: dict(v=v, c=6, b=5, h=h, m=6 * h, Q=[2, 3], H1=H1, q=q, s=s, f=f,
                                        mods=[2, 3] + H1 + [q], nQ=2, U=U)


def canon():
    for cfgd, want in ((C6(35, 5, 5, 2, 2, [], []), 1728), (C6(47, 7, 7, 2, 3, [], []), 36), (C6(47, 7, 7, 3, 2, [], []), 1728)):
        cfg = SD.cfg_from_json(json.dumps(cfgd)); Z = ZOrb(cfg)
        C = Z.compressions(); cls, idx = Z.classes(C); qd = Z.class_quads(cls)
        byc = {}
        for z in C: byc.setdefault(idx[z], []).append(z)
        canon_ = set(); ordered = 0
        for q in qd:
            for perm in set(itertools.permutations(q['m'])):
                for zq in itertools.product(*[byc[i] for i in perm]):
                    n, _, Ys = Z.exact_count(list(zq)); ordered += n
                    for Y in Ys:
                        canon_.add(tuple(sorted(max(tuple(int(x) for x in Y[k]), tuple(-int(x) for x in Y[k])) for k in range(4))))
        tot = sum(q['ordered'] * Z.n0(cls, q['m']) for q in qd)        # H' = 1: the lag-0 count is the stage-1 count
        print(json.dumps(dict(v=cfgd['v'], q=cfgd['q'], s=cfgd['s'], ordered=ordered, class_weighted=tot, canonical=len(canon_),
                              expected_canonical=want, ok=len(canon_) == want and ordered == tot)))


def sd():
    cfg = SD.cfg_from_json(open(os.path.join(ROOT, 'logs/ab/cfg_347_s6.json')).read())
    Z = ZOrb(cfg); C = Z.compressions(); cls, idx = Z.classes(C); qd = Z.class_quads(cls); S = SDCount(Z); bad = 0
    for t, q in enumerate(qd):
        zq = [cls[i]['rep'] for i in q['m']]
        a = S.count(zq)[0]; b = Z.exact_count(zq)[0]; bad += a != b
        print(json.dumps(dict(i=t, m=list(q['m']), sd=a, direct=b)))
    print(json.dumps(dict(multisets=len(qd), mismatches=bad)))


def invar():
    cfg = SD.cfg_from_json(open(os.path.join(ROOT, 'logs/ab/cfg_347_s6.json')).read())
    Z = ZOrb(cfg); C = Z.compressions(); cls, idx = Z.classes(C); qd = Z.class_quads(cls)
    byc = {}
    for z in C: byc.setdefault(idx[z], []).append(z)
    random.seed(5)
    for t in (25, 26, 4, 12):
        m = qd[t]['m']; rep = Z.exact_count([cls[i]['rep'] for i in m])[0]; vals = []
        perms = list(set(itertools.permutations(m)))
        for _ in range(6):
            perm = random.choice(perms); vals.append(Z.exact_count([random.choice(byc[i]) for i in perm])[0])
        print(json.dumps(dict(multiset=t, rep=rep, members=vals, ok=all(v == rep for v in vals))))


def screens():
    from ab_verify_published import load_core, recoordinate
    cfg = SD.cfg_from_json(open(os.path.join(ROOT, 'logs/ab/cfg_347_s6.json')).read()); sl = W.build(cfg)
    hdr, X0 = load_core(os.path.join(ROOT, 'data/core_v347_b5.txt'))
    X = recoordinate(X0, (2, 3, 57), (2, 3, 3, 19), [(0, 2), (1, 3), (2, 3), (2, 19)])
    assert sl.is_invariant(X) and (sl.G.paf(X).sum(0) == sl.T).all()
    for keep in ([], [0], [1], [0, 1]):
        S = Screen(cfg, keep=keep, sl=sl); G = sl.G
        proj = S.PZ.idx(np.concatenate([G.co[:, keep], G.co[:, -1:]], 1)) if keep else G.co[:, -1]
        Xb = np.zeros((4, S.PZ.n), np.int64)
        for k in range(4): np.add.at(Xb[k], proj, X[k])
        print(json.dumps(dict(check='347 core identity', keep=keep, ok=bool((S.PZ.paf(Xb).sum(0) == S.T_full).all()))))
    for (n, q) in ((3, 7), (3, 5), (1, 13)):
        cfgd = C6(6 * n * q + 5, n * q, q, 2, (q - 1) // 2, [n] if n > 1 else [], [["neg", [1]]])
        cfg = SD.cfg_from_json(json.dumps(cfgd)); sl = W.build(cfg)
        for seed in range(3):
            rng = np.random.default_rng(seed)
            Xp = rng.choice(np.array([-1, 1]), size=(4, sl.r))[:, sl.orb]; TG = sl.G.paf(Xp).sum(0)
            Y = sl.compress(Xp); yb = [Y[k][sl.breps] for k in range(4)]
            for keep in ([], [0], [1], [0, 1]):
                S = Screen(cfg, keep=keep, sl=sl); S.set_target(TG)
                c1 = S.count(yb); c2 = S.count_exact_small(yb)
                print(json.dumps(dict(check='planted', v=cfgd['v'], seed=seed, keep=keep, hashed=c1, exact=c2, ok=abs(c1 - c2) < 0.5 and c2 >= 1)))
    cfg = SD.cfg_from_json(json.dumps(C6(47, 7, 7, 2, 3, [], []))); sl = W.build(cfg)
    X = np.load(os.path.join(ROOT, 'logs/ab/mono6/cores/core_v47_q7_s2_99c1dc75_31561075d0a6d80e_1.npy'))
    assert (sl.G.paf(X).sum(0) == sl.T).all()
    Y = sl.compress(X); yb = [Y[k][sl.breps] for k in range(4)]
    for keep in ([], [0], [1], [0, 1]):
        print(json.dumps(dict(check='v = 47 control', keep=keep, count=Screen(cfg, keep=keep, sl=sl).count(yb))))


if __name__ == '__main__':
    dict(canon=canon, sd=sd, invar=invar, screens=screens)[sys.argv[1]]()
