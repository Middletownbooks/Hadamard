r"""Exact type-level compression screens of stage 2 over ALL stage-1 candidates of a configuration [M].

For every class multiset of compression quadruples with a non-zero lag-0 count (theory/s1z_lib.py), every block type
compatible with its compression (per coset component, multisets of letters per orbit class with the component's coset
sums; theory/s1z_qzq.py for why the type is all that matters), every type quadruple with sum_k E_k = T1(0); each type
quadruple is screened by the exact compression conditions of stage 2 (theory/s1z_qzq.Screen) at the levels
Qbar x Z_q for Qbar = Q compressed to the kept axes: [] (the Z_q level), each single axis, all axes (Q x Z_q), in that
order; a zero count at any level kills the type quadruple, hence every stage-1 solution with those types [M].
If every type quadruple of every multiset dies, the configuration has no core [M].
Output JSONL: one line per multiset (type quadruples, deaths per level, survivors), the survivors themselves in
OUT.survivors.jsonl (canonical letters per orbit for each block, and the counts per level).
usage: python3 theory/s1z_typescreen.py CFG.json OUT.jsonl [--levels '[[],[0],[1],[0,1]]'] [--max-tq 2000000] [--part i/W]
(restart-safe: multisets already in OUT are skipped; with --part, OUT and the survivors file get a .pI suffix)"""
import sys, os, json, math, itertools, time
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_s1sd as SD
from s1z_lib import ZOrb
from s1z_qzq import Screen


def comp_classes(Z, scr, i):
    """orbit classes of component i: list of (orbits, profile within the component, weight, letters)"""
    C = Z.comps[i]; cl = defaultdict(list)
    for j in Z.corbs[i]:
        key = (tuple(int(x) for x in Z.prof[j, C]), int(Z.w[j]), scr.ocls[j], tuple(Z.letters[j]))
        cl[key].append(j)
    return [(orbs, key[0], key[1], key[3]) for key, orbs in sorted(cl.items())]


def comp_types(Z, classes, zc, cap):
    """letter assignments (canonical: sorted letters per class) with component compression zc, energy <= cap"""
    out = []
    def rec(t, acc_z, acc_e, cur):
        if t == len(classes):
            if tuple(acc_z) == tuple(zc): out.append((dict(cur), acc_e))
            return
        orbs, prof, w, lets = classes[t]
        for ms in itertools.combinations_with_replacement(lets, len(orbs)):
            s = sum(ms); e = acc_e + w * sum(y * y for y in ms)
            if e > cap: continue
            nz = [a + s * p for a, p in zip(acc_z, prof)]
            for j, y in zip(orbs, ms): cur[j] = y
            rec(t + 1, nz, e, cur)
            for j in orbs: cur.pop(j, None)
    rec(0, [0] * len(zc), 0, {})
    return out


def main():
    a = sys.argv[1:]
    cfg = SD.cfg_from_json(open(a[0]).read() if a[0].endswith('.json') else a[0]); outp = a[1]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    nQ = cfg['nQ']
    levels = json.loads(opt('--levels', json.dumps([[]] + [[i] for i in range(nQ)] + [list(range(nQ))])))
    max_tq = int(opt('--max-tq', '2000000'))
    part, nW = (int(x) for x in opt('--part', '0/1').split('/'))
    base = outp
    if nW > 1: outp = outp.replace('.jsonl', f'.p{part}.jsonl')
    import glob
    done = set()
    for f in {base, outp} | set(glob.glob(base.replace('.jsonl', '.p*.jsonl'))):
        if os.path.exists(f) and not f.endswith('.survivors.jsonl'):
            for l in open(f):
                r = json.loads(l)
                if r.get('stage') == 'ms': done.add(r['i'])
    t0 = time.time()
    Z = ZOrb(cfg)
    Cz = Z.compressions(); cls, idx = Z.classes(Cz); qd = Z.class_quads(cls)
    cache = {}
    for q in qd: q['n0'] = Z.n0(cls, q['m'], cache)
    qd = [q for q in qd if q['n0'] > 0]
    scr = [Screen(cfg, keep=k, sl=Z.sl) for k in levels]
    ccls = [comp_classes(Z, scr[0], i) for i in range(len(Z.comps))]
    ctcache = {}
    def block_types(z):
        z = np.asarray(z)
        lists = []
        for i, C in enumerate(Z.comps):
            zc = tuple(int(z[q]) for q in C)
            key = (i, zc)
            if key not in ctcache: ctcache[key] = comp_types(Z, ccls[i], zc, Z.Ecap)
            lists.append(ctcache[key])
        byE = defaultdict(list)
        for combo in itertools.product(*lists):
            E = sum(e for _, e in combo)
            if E > Z.Ecap: continue
            yb = np.zeros(Z.L1, np.int64)
            for d, _ in combo:
                for j, y in d.items(): yb[j] = y
            byE[E].append(yb)
        return byE
    fo = open(outp, 'a'); fs = open(outp.replace('.jsonl', '') + '.survivors.jsonl', 'a')
    head = dict(stage='setup', multisets=len(qd), levels=levels, L1=Z.L1, T0=Z.T0, part=f'{part}/{nW}', secs=round(time.time() - t0, 1))
    print(json.dumps(head), flush=True); fo.write(json.dumps(head) + "\n")
    tot_tq = 0; tot_surv = 0
    for mi, q in enumerate(qd):
        if mi % nW != part or mi in done: continue
        zq = [cls[i]['rep'] for i in q['m']]
        B = [block_types(z) for z in zq]
        dead = [0] * len(levels); surv = 0; ntq = 0
        E01 = defaultdict(list)
        for e0 in B[0]:
            for e1 in B[1]: E01[e0 + e1].append((e0, e1))
        for e2 in B[2]:
            for e3 in B[3]:
                for (e0, e1) in E01.get(Z.T0 - e2 - e3, ()):
                    for Y4 in itertools.product(B[0][e0], B[1][e1], B[2][e2], B[3][e3]):
                        ntq += 1
                        if ntq > max_tq: raise RuntimeError('too many type quadruples')
                        cnts = []
                        for li, S in enumerate(scr):
                            c = S.count(Y4); cnts.append(c)
                            if c == 0: dead[li] += 1; break
                        else:
                            surv += 1
                            fs.write(json.dumps(dict(ms=mi, m=list(q['m']), z=[list(z) for z in zq],
                                                     y=[yb.tolist() for yb in Y4], counts=cnts,
                                                     log2_lifts_qzq=round(math.log2(cnts[-1]), 2))) + "\n")
        tot_tq += ntq; tot_surv += surv
        rec = dict(stage='ms', i=mi, m=list(q['m']), ordered=q['ordered'], log2N0=round(math.log2(q['n0']), 3),
                   type_quads=ntq, dead_per_level=dead, survivors=surv, secs=round(time.time() - t0, 1))
        print(json.dumps(rec), flush=True); fo.write(json.dumps(rec) + "\n"); fo.flush(); fs.flush()
    summ = dict(stage='summary', type_quads=tot_tq, survivors=tot_surv, secs=round(time.time() - t0, 1))
    print(json.dumps(summ), flush=True); fo.write(json.dumps(summ) + "\n"); fo.close(); fs.close()


if __name__ == '__main__':
    main()
