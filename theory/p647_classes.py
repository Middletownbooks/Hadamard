r"""True symmetry classes of the joint solutions of y_B #0 (647 x353): the stabiliser of y_B including per-block
reversal [P, M].

The group of theory/p647_common.group_elements, (g x)_k(z) = s_k x_{pi(k)}(u z + t_k), applies the multiplier u to all
four blocks at once.  Reversal of a single block, x_k(z) -> x_k(-z), is a symmetry too [P]:
  * PAF(x(-.))(h) = PAF(x)(-h) = PAF(x)(h), so sum_k PAF(x_k) is unchanged and cores map to cores;
  * -1 commutes with 353, so x353-invariance is kept;
  * on the Z_34 compression (y_B) reversal is the identity, because -1 = 13^2 lies in <353 mod 17> = <13>; so every
    per-block reversal fixes every y_B, and only y_A (the Z_38 compression) moves, by z -> -z in that block.
So the stabiliser of y_B in the full group is the old one times the 8 per-block reversal patterns (eps_0 = 1, since a
global -1 is already a multiplier), and joint solutions in one orbit of it have isomorphic cross levels (the per-layer
lists, distinct sigma, S-signatures and D-contributions correspond under z -> -z, so the same counts of S-solutions,
D-candidates and cores).

This script computes the orbits of the 20,928 joint solutions under the full stabiliser, checks the isomorphism on
real data (sorted S-signature and D-contribution lists of a joint solution and of its image agree), and maps the
production units (logs/p647/prod/units.json) and the old 3,488 orbits to the classes.
Output: logs/p647/joint_yb0_classes.json (class_of per joint solution, representatives, class sizes, unit -> classes).
usage: python3 theory/p647_classes.py [--check N]"""
import sys, os, json, itertools, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p647_common as PC

def group_full(m=PC.M_):
    """(u, eps, t, s, pi): (g x)_k(z) = s_k x_{pi(k)}(eps_k u z + t_k), eps_0 = 1"""
    for (u, t, s, pi) in PC.group_elements(m):
        for e in itertools.product((1, -1), repeat=3):
            yield u, (1,) + e, t, s, pi

def act_full(y, g, F):
    u, e, t, s, pi = g; y = np.asarray(y); b = np.arange(F)
    return np.array([s[k] * y[pi[k]][(e[k] * u * b + t[k]) % F] for k in range(4)])

def stabiliser_full(yb, FB=34, FA=38, m=PC.M_):
    """stabiliser of y_B (on Z_FB) in the full group, as distinct actions on Z_FA"""
    yb = np.asarray(yb); seenB = set(); keysB = set()
    for g in group_full(m):
        kb = (g[0] % FB, g[1], tuple(x % FB for x in g[2]), g[3], g[4])
        if kb in seenB: continue
        seenB.add(kb)
        if np.array_equal(act_full(yb, g, FB), yb): keysB.add(kb)
    acts = {}
    for g in group_full(m):
        if (g[0] % FB, g[1], tuple(x % FB for x in g[2]), g[3], g[4]) in keysB:
            acts.setdefault((g[0] % FA, g[1], tuple(x % FA for x in g[2]), g[3], g[4]), g)
    return list(acts.values())

def classes(X, sols, stab, FA=38):
    keys = []
    for s_ in sols:
        ya = np.array(X.J.ya_full(s_['layers']))
        keys.append(min(tuple(act_full(ya, g, FA).ravel().tolist()) for g in stab))
    rep = {}
    for i, kk in enumerate(keys): rep.setdefault(kk, i)
    cid = {kk: c for c, kk in enumerate(sorted(rep, key=lambda kk: rep[kk]))}
    return [cid[kk] for kk in keys], [rep[kk] for kk in sorted(rep, key=lambda kk: rep[kk])]

def image_index(X, sols, j, g, FA=38):
    """index of the joint solution g(sols[j]) (None if it is not among sols)"""
    if not hasattr(image_index, "tab"):
        image_index.tab = {tuple(np.array(X.J.ya_full(s_['layers'])).ravel().tolist()): i for i, s_ in enumerate(sols)}
    ya = np.array(X.J.ya_full(sols[j]['layers']))
    return image_index.tab.get(tuple(act_full(ya, g, FA).ravel().tolist()))

def check_iso(X, ys, sols, j, k, log=print):
    """reverse block k of joint solution j: the image joint solution's block k has the reversed layer lists, the same
    sorted S-signatures and the same sorted D-contribution multiset (the other blocks are identical)"""
    import p647_cross as PX
    from nd_cross import orbit_maps
    reps, oi, maps = orbit_maps(X); NH = PX.lag_tensor(X, oi)
    e = [1, 1, 1, 1]; e[k] = -1
    g = (1, tuple(e), (0, 0, 0, 0), (1, 1, 1, 1), (0, 1, 2, 3))
    j2 = image_index(X, sols, j, g)
    if j2 is None: return dict(j=j, k=k, image=None)
    B1 = PX.build(X, ys[0], X.J.ya_full(sols[j]['layers']), maps, NH, oi)
    B2 = PX.build(X, ys[0], X.J.ya_full(sols[j2]['layers']), maps, NH, oi)
    ro = np.array([oi[(-int(r)) % X.mp] for r in reps])
    f = lambda A: sorted(map(bytes, np.ascontiguousarray(A).astype(np.int8)))
    lists_rev = f(B1[k].A0[:, ro]) == f(B2[k].A0) and f(B1[k].A1[:, ro]) == f(B2[k].A1)
    others = all(np.array_equal(B1[q].A0, B2[q].A0) and np.array_equal(B1[q].A1, B2[q].A1) for q in range(4) if q != k)
    ss = lambda B: np.array(sorted(map(tuple, B.Ssig.tolist())))
    sig_same = np.array_equal(ss(B1[k]), ss(B2[k]))
    srt = lambda B: np.sort(np.concatenate([B.mult]))
    mult_same = np.array_equal(srt(B1[k]), srt(B2[k]))
    return dict(j=j, k=k, image=int(j2), lists_reversed=bool(lists_rev), other_blocks_identical=bool(others),
                S_signatures_same=bool(sig_same), sigma_multiplicities_same=bool(mult_same))

def main(ncheck=3, log=print):
    X, ys = PC.load_setup(); rec, sols = PC.joint_solutions(X, ys, 0)
    t0 = time.time()
    old = PC.stabiliser_on(ys[0], PC.M_, 34, 38)
    stab = stabiliser_full(ys[0])
    log(json.dumps(dict(stage="stabiliser", old_on_Z38=len(old), full_on_Z38=len(stab), secs=round(time.time() - t0, 1))))
    cls, reps = classes(X, sols, stab)
    sizes = np.bincount(cls)
    old_orb = json.load(open(os.path.join(PC.LOG, "joint_yb0_orbits.json")))['orbit_of']
    orb_per_class = {}
    for i, c in enumerate(cls): orb_per_class.setdefault(c, set()).add(old_orb[i])
    U = json.load(open(os.path.join(PC.LOG, "prod", "units.json")))['units']
    unit_classes = [sorted({cls[e] for row in u['joint'] for e in row}) for u in U]
    H = [sols[r]['mult_log2'] for r in reps]
    out = dict(joint=len(sols), stabiliser_full_on_Z38=len(stab), stabiliser_old_on_Z38=len(old), classes=len(reps),
               class_size_hist={int(s): int((sizes == s).sum()) for s in np.unique(sizes)},
               old_orbits_per_class_hist={int(n): sum(1 for v in orb_per_class.values() if len(v) == n)
                                          for n in sorted({len(v) for v in orb_per_class.values()})},
               units=len(U), classes_per_unit_hist={int(n): sum(1 for uc in unit_classes if len(uc) == n)
                                                    for n in sorted({len(uc) for uc in unit_classes})},
               classes_covered_by_units=len({c for uc in unit_classes for c in uc}),
               H_hist={str(h): int(sum(1 for x in H if x == h)) for h in sorted(set(H))})
    log(json.dumps(out))
    json.dump(dict(out, class_of=cls, reps=reps, H=H, unit_classes=unit_classes), open(os.path.join(PC.LOG, "joint_yb0_classes.json"), "w"))
    rng = np.random.default_rng(1)
    for t in range(ncheck):
        r = check_iso(X, ys, sols, int(rng.integers(len(sols))), int(rng.integers(1, 4)))
        log(json.dumps(dict(stage="iso_check", **r)))
    return out

if __name__ == "__main__":
    a = sys.argv[1:]
    main(int(a[a.index("--check") + 1]) if "--check" in a else 3, log=lambda s: print(s, flush=True))
