r"""H(4892) slice, step 2 driver: exact stage 2 for every stage-1 solution of data/h4892_stage1.json [M].

For each of the 384 stage-1 classes (up to block order and block sign, which map stage-2 fibres bijectively onto
each other [P]):
  1. per block, the complete class-sign list (h4892_stage2.block_list), its 15-class PAF vectors over G, and its
     S-level s = A P (A = h4892_slevel.smap(L), the compression along L; L = G' by default, so s is the PAF of the
     compression to Z13 on its two nonzero lag classes);
  2. pairing {a, b} | {c, d} minimising |L_a||L_b| + |L_c||L_d|; pivot = the largest block, restricted to the least
     pattern of each <z -> 2z> orbit (h4892_stage2 docstring) -- complete up to that symmetry [P];
  3. S-level 2+2 join over the distinct s values (exact, dictionary);
  4. for every matched pair of S-buckets: all left pairs and all right pairs of those buckets, a 64-bit linear hash
     (mod 2^64, exact on sums) of the 15-vector, sort + searchsorted, every hash hit re-checked on the exact vectors;
  5. each hit (a core) is expanded to G and checked by the definition on all 1221 nonzero lags.
Planted controls (--planted N): a random stage-1 class, random patterns from the four full lists, the target
replaced by their PAF sum; recovery means some <z -> 2z> image of the planted quadruple is returned.
Logs one JSON line per y to --log.
usage: python3 theory/h4892_stage2_run.py [--planted N] [--L G'|Z47|Z2] [--log FILE] [--only i,j,...]"""
import sys, os, json, time, math
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h4892_stage1 as S1, h4892_stage2 as T, h4892_slevel as SL

RNG = np.random.default_rng(20261106)
HV = RNG.integers(-2 ** 62, 2 ** 62, size=T.NL, dtype=np.int64)


def v6(y): return [y[o[0]] for o in S1.ORB]


class Block:
    def __init__(self, v, A):
        self.v = v; self.U = T.block_list(v); self.P = T.paf_vec(self.U)
        self.H = (self.P * HV).sum(1)                         # wraps mod 2^64: linear and exact on sums
        self.s = self.P @ A.T
        self.minmask = T.orbit_min_mask(self.U)
    def buckets(self, mask=None):
        d = defaultdict(list)
        idx = np.arange(len(self.U)) if mask is None else np.nonzero(mask)[0]
        for i in idx: d[tuple(self.s[i])].append(i)
        return {k: np.array(v, np.int64) for k, v in d.items()}


def solve(blocks, Tvec, A, pivot=True):
    """blocks: four Block; Tvec: 15-vector target; returns (cores as index quadruples, stats)"""
    nb = [len(b.U) for b in blocks]
    pairings = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
    (a, b), (c, d) = min(pairings, key=lambda p: nb[p[0][0]] * nb[p[0][1]] + nb[p[1][0]] * nb[p[1][1]])
    piv = int(np.argmax(nb)) if pivot else -1
    Bk = [blk.buckets(blk.minmask if k == piv else None) for k, blk in enumerate(blocks)]
    TS = tuple(int(x) for x in (A @ Tvec))
    def pairbk(x, y):
        out = defaultdict(list)
        for kx in Bk[x]:
            for ky in Bk[y]: out[tuple(p + q for p, q in zip(kx, ky))].append((kx, ky))
        return out
    LB, RB = pairbk(a, b), pairbk(c, d)
    tH = np.int64((np.array(Tvec, np.int64) * HV).sum())
    def mk(x, y, lst):
        I, J = [], []
        for kx, ky in lst:
            g, h = np.meshgrid(Bk[x][kx], Bk[y][ky], indexing='ij'); I.append(g.ravel()); J.append(h.ravel())
        return np.concatenate(I), np.concatenate(J)
    nL = nR = nmatch = 0; hits = []; nhash = 0
    for ks, lstL in LB.items():
        need = tuple(t - x for t, x in zip(TS, ks))
        if need not in RB: continue
        nmatch += 1
        I, J = mk(c, d, RB[need])
        rH = blocks[c].H[I] + blocks[d].H[J]; o = np.argsort(rH, kind='stable'); rs = rH[o]; nR += len(rs)
        A_, B_ = mk(a, b, lstL)
        for c0 in range(0, len(A_), 1 << 23):
            ia, ib = A_[c0:c0 + (1 << 23)], B_[c0:c0 + (1 << 23)]
            nd = tH - (blocks[a].H[ia] + blocks[b].H[ib]); nL += len(nd)
            pos = np.searchsorted(rs, nd); pos[pos >= len(rs)] = len(rs) - 1
            for h in np.nonzero(rs[pos] == nd)[0]:
                lo = np.searchsorted(rs, nd[h], 'left'); hi = np.searchsorted(rs, nd[h], 'right')
                for qq in o[lo:hi]:
                    nhash += 1
                    q4 = {a: ia[h], b: ib[h], c: I[qq], d: J[qq]}
                    tot = sum(blocks[k].P[q4[k]] for k in range(4))
                    if (tot == np.array(Tvec)).all(): hits.append(tuple(int(q4[k]) for k in range(4)))
    st = dict(sizes=nb, pairing=[[a, b], [c, d]], pivot=piv, pivot_size=int(len(np.nonzero(blocks[piv].minmask)[0])) if piv >= 0 else None,
              s_distinct=[len(x) for x in Bk], s_buckets_left=len(LB), s_buckets_right=len(RB), s_matched=nmatch,
              left_pairs=int(nL), right_pairs=int(nR), hash_hits=nhash, cores=len(hits))
    return hits, st


def check_core(blocks, q4):
    """definition check on G (all 1221 nonzero lags), independent of the class bookkeeping"""
    X = [T.expand(blocks[k].U[q4[k]]) for k in range(4)]
    tot = np.zeros(T.NG, np.int64)
    for x in X:
        for h in range(T.NG):
            tot[h] += sum(int(x[g]) * int(x[T.addG(g, h)]) for g in range(T.NG))
    return bool(tot[0] == 4888 and (tot[1:] == -4).all()), X


def main():
    a = sys.argv[1:]
    npl = int(a[a.index('--planted') + 1]) if '--planted' in a else 0
    Ln = a[a.index('--L') + 1] if '--L' in a else "G'"
    logp = a[a.index('--log') + 1] if '--log' in a else None
    only = [int(x) for x in a[a.index('--only') + 1].split(',')] if '--only' in a else None
    A = SL.smap(Ln)
    D = json.load(open('data/h4892_stage1.json'))
    Ys = D['canonical']
    cache = {}
    def blk(y):
        v = tuple(v6(y))
        if v not in cache: cache[v] = Block(list(v), A)
        return cache[v]
    Tvec = np.full(T.NL, -4, np.int64)
    logf = open(logp, 'a') if logp else None
    # planted controls first
    for t in range(npl):
        # the pivot restriction is complete for targets invariant under z -> 2z (the real target, constant -4);
        # a planted target is not invariant, so the planted pivot pattern is drawn among the pivot's orbit minima
        # (even plants) or the pivot is switched off (odd plants); recovery = the exact planted quadruple returned
        i = int(RNG.integers(len(Ys))); bl = [blk(y) for y in Ys[i]]
        usepiv = (t % 2 == 0); piv = int(np.argmax([len(b.U) for b in bl]))
        q4 = []
        for k, b in enumerate(bl):
            pool = np.nonzero(b.minmask)[0] if (usepiv and k == piv) else np.arange(len(b.U))
            q4.append(int(pool[int(RNG.integers(len(pool)))]))
        Tp = sum(bl[k].P[q4[k]] for k in range(4))
        t0 = time.time(); hits, st = solve(bl, Tp, A, pivot=usepiv)
        rec = dict(planted=t, y=i, q4=q4, pivot_on=usepiv, recovered=tuple(q4) in set(hits),
                   secs=round(time.time() - t0, 2), **st)
        print(json.dumps(rec), flush=True)
        if logf: logf.write(json.dumps(rec) + "\n"); logf.flush()
    if npl and '--planted-only' in a: return
    total = dict(y=0, cores=0, left_pairs=0, right_pairs=0, s_matched=0)
    for i, Y in enumerate(Ys):
        if only is not None and i not in only: continue
        bl = [blk(y) for y in Y]
        t0 = time.time(); hits, st = solve(bl, Tvec, A)
        cores = []
        for h in hits:
            ok, X = check_core(bl, h); cores.append(dict(q4=list(h), verified=ok,
                                                       seqs=[x.tolist() for x in X]))
        rec = dict(y=i, L=Ln, secs=round(time.time() - t0, 2), **st, core_list=cores)
        print(json.dumps({k: v for k, v in rec.items() if k != 'core_list'}), flush=True)
        if logf: logf.write(json.dumps(rec) + "\n"); logf.flush()
        total['y'] += 1; total['cores'] += len(cores)
        for k in ('left_pairs', 'right_pairs', 's_matched'): total[k] += st[k]
    print('TOTAL', json.dumps(total), flush=True)
    if logf: logf.write(json.dumps(dict(total=total, L=Ln)) + "\n")


if __name__ == "__main__":
    main()
