r"""Random stage-1 solutions of a GP slice (b >= 3) with a prescribed Q-compression quadruple [M for every solution
returned: verified exactly; the search itself is a random sample, not exhaustive].

Per block k, the alpha vectors (one sign per G'-orbit) whose Q-compression is congruent to z_k (mod 2 f0) are the
product over the cosets of Q of the per-coset sign patterns with the right residue.  Random alpha-triples
(alpha_0, alpha_1, alpha_2) are completed by a hash lookup of alpha_3 with
sum_k PAF(alpha_k) = T1 (mod 4 f0) at every lag (theory/ab_s1sd.py, alpha-congruence), and every alpha-quadruple
whose rows are admissible goes through the exact S/D joins of ab_s1sd.S1SD (block, dlist, ejoin4), restricted to the
compressions z_0..z_3.  Output: one JSONL line per solution (y, z, alpha-quadruple).
usage: python3 theory/ab_s1_sample.py '<cfg json>' '<z json 4 x c>' OUT.jsonl [--seed 1] [--secs 600] [--want 10]
                                      [--memgb 3.2]"""
import sys, os, json, time, itertools, resource
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_s1sd as SD
from ab_s1sd import ejoin4


def alpha_lists(P, Z):
    sl = P.sl; qc = P.qfull; f0 = P.f0; md = 2 * f0
    nO = len(sl.breps); w = sl.bsz
    # orbit -> its coset profile (an orbit may meet several cosets when U moves Q)
    prof = np.zeros((nO, Z.shape[1]), np.int64)
    for g in range(len(sl.borb)): prof[sl.borb[g], qc[g]] += 1
    # group orbits by support (cosets met); per group enumerate sign patterns and their coset-sum vectors
    groups = defaultdict(list)
    for j in range(nO): groups[tuple(np.nonzero(prof[j])[0].tolist())].append(j)
    out = []
    for k in range(4):
        parts = []
        for supp, orbs in groups.items():
            pats = []
            for sg in itertools.product((1, -1), repeat=len(orbs)):
                zs = sum(s_ * prof[j] for s_, j in zip(sg, orbs))
                if all((int(zs[q]) - int(Z[k, q])) % md == 0 for q in supp): pats.append(sg)
            if not pats: return None
            parts.append((orbs, pats))
        tot = int(np.prod([len(p) for _, p in parts]))
        A = np.zeros((tot, nO), np.int64)
        idx = np.indices([len(p) for _, p in parts]).reshape(len(parts), -1).T
        for pi, (orbs, pats) in enumerate(parts):
            pa = np.array(pats, np.int64)[idx[:, pi]]
            A[:, orbs] = pa
        out.append(A)
    return out


def keys_of(P, A, wv):
    full = A[:, P.sl.borb]
    mod = 4 * P.f0
    Pf = np.rint(P.Gp.paf(full)).astype(np.int64)[:, P.reps] % mod
    return (Pf @ wv), Pf, full.sum(1)


def main():
    a = sys.argv[1:]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    memgb = float(opt('--memgb', '3.2'))
    resource.setrlimit(resource.RLIMIT_AS, (int(memgb * 2 ** 30),) * 2)
    cfg = SD.cfg_from_json(a[0]); Z = np.array(json.loads(a[1]), np.int64); outp = a[2]
    seed = int(opt('--seed', '1')); secs = float(opt('--secs', '600')); want = int(opt('--want', '10'))
    rng = np.random.default_rng(seed)
    P = SD.S1SD(cfg, zallow=Z.tolist())
    rq, Ainfo = SD.row_quads_of(cfg, sl=P.sl)
    mod, rmod = 4 * P.f0, 2 * P.f0
    okrows = set()
    for q in set(tuple(sorted(x)) for x in rq):
        for perm in set(itertools.permutations(q)):
            for sg in itertools.product((1, -1), repeat=4):
                okrows.add(tuple((s_ * r) % rmod for s_, r in zip(sg, perm)))
    AL = alpha_lists(P, Z)
    if AL is None: print(json.dumps(dict(stage='alpha', empty=True))); return
    wv = rng.integers(1, 2 ** 61, size=len(P.reps), dtype=np.int64)
    K = []; PF = []; R = []
    for A in AL:
        h, Pf, r = keys_of(P, A, wv); K.append(h); PF.append(Pf); R.append(r)
    T = P.T1[P.reps] % mod
    look = defaultdict(list)
    for i, hh in enumerate(K[3].tolist()): look[hh].append(i)
    print(json.dumps(dict(stage='alpha', sizes=[len(x) for x in AL], keys3=len(look))), flush=True)
    t0 = time.time(); tries = 0; quads = 0; nsol = 0; seen = set()
    TS0, TD0 = int(P.TS[0]), int(P.TD[0]); TS = P.TS[P.preps]; TD = P.TD[P.preps]
    fo = open(outp, 'a')
    while time.time() - t0 < secs and nsol < want:
        B = 4096
        i0 = rng.integers(len(AL[0]), size=B); i1 = rng.integers(len(AL[1]), size=B); i2 = rng.integers(len(AL[2]), size=B)
        need = (T[None, :] - PF[0][i0] - PF[1][i1] - PF[2][i2]) % mod
        hn = need @ wv
        tries += B
        for t in range(B):
            for i3 in look.get(int(hn[t]), ()):
                if not ((PF[0][i0[t]] + PF[1][i1[t]] + PF[2][i2[t]] + PF[3][i3]) % mod == T).all(): continue
                rows = tuple(int(R[k][ii]) % rmod for k, ii in enumerate((i0[t], i1[t], i2[t], i3)))
                if rows not in okrows: continue
                qd = (int(i0[t]), int(i1[t]), int(i2[t]), int(i3))
                if qd in seen: continue
                seen.add(qd); quads += 1
                blocks = [P.block(AL[k][ii], TS0, TD0) for k, ii in enumerate(qd)]
                for (ia, ib, ic, id_) in ejoin4(*[b_['Ssig'] for b_ in blocks], TS):
                    us = [blocks[0]['U'][ia], blocks[1]['U'][ib], blocks[2]['U'][ic], blocks[3]['U'][id_]]
                    Ds = [P.dlist(b_, u) for b_, u in zip(blocks, us)]
                    for (da, db, dc, dd) in ejoin4(*[x[1] for x in Ds], TD):
                        ws = [Ds[0][0][da], Ds[1][0][db], Ds[2][0][dc], Ds[3][0][dd]]
                        Y = [P.y_of(B_, u, w_).tolist() for B_, u, w_ in zip(blocks, us, ws)]
                        if not P.verify(Y): continue
                        nsol += 1
                        fo.write(json.dumps(dict(y=Y, z=Z.tolist(), alpha_quad=list(qd))) + "\n"); fo.flush()
        if tries % (B * 64) == 0:
            print(json.dumps(dict(stage='progress', tries=tries, alpha_quads=quads, solutions=nsol, secs=round(time.time() - t0, 1))), flush=True)
    print(json.dumps(dict(stage='end', tries=tries, alpha_quads=quads, solutions=nsol, secs=round(time.time() - t0, 1))), flush=True)


if __name__ == "__main__":
    main()
