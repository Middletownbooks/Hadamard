r"""Exact compression screens of stage 2 for whole compression quadruples (all block types at once) [P, M].

theory/s1z_qzq.py: for every core, sum_k PAF_{Qbar x Z_q}(xbar_k) = T' at each compression level Qbar x Z_q, and the
lifts of a stage-1 block contribute orbit by orbit.  Here the types are not enumerated: for a compression z_k, a DP over
the G'-orbits of block k (state = partial coset sums, partial stage-1 energy, partial xbar at the level's orbit
representatives; one transition per letter and fibre pattern) gives the set B_k of pairs (E, xbar) over ALL letter
assignments with compression z_k and ALL their lifts.  A core with this compression quadruple needs (E_k, xbar_k) in B_k
with sum_k E_k = T1(0) (stage-1 lag 0) and sum_k PAF(xbar_k) = T' (the level): a 2+2 join on (E, 64-bit hash of the
PAF signature) decides it; hash collisions can only create false matches, so "no match" is exact [M].  A dead
quadruple has no core, whatever its stage-1 solutions are.
usage: python3 theory/s1z_zqscreen.py CFG.json OUT.jsonl [--levels '[[],[0],[1]]'] [--max-states 3e7] [--part i/W]"""
import sys, os, json, math, time, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_s1sd as SD, ab_gauss as AGs
from s1z_lib import ZOrb
from s1z_qzq import Screen


class AggScreen:
    def __init__(self, Z, keep, max_states=3 * 10 ** 7):
        self.Z = Z; self.S = Screen(Z.cfg, keep=keep, sl=Z.sl); S = self.S
        self.max_states = max_states
        self.nr = len(S.oreps)
        # per orbit: list of (letter, z-contribution over Q cosets, E contribution, xbar contribution (unique, all patterns))
        self.trans = []
        for j in range(Z.L1):
            opts = []
            for y in Z.letters[j]:
                pats = S.patterns(j, y)
                for xb in pats: opts.append((int(y), tuple(int(v) for v in xb)))
            self.trans.append(opts)
        # orbit order: by component, so finished components can be checked and dropped
        self.order = [j for C in Z.corbs for j in C]
        self.comp_end = {}
        pos = 0
        for i, C in enumerate(Z.corbs):
            pos += len(C); self.comp_end[pos - 1] = i

    def block_set(self, z):
        """distinct (E, signature-hash) pairs of block lifts with compression z (all letter assignments, all lifts)"""
        Z = self.Z; S = self.S; c = Z.c; nr = self.nr
        z = np.asarray(z, np.int64)
        # state columns: [z over the current component's cosets..., E, xbar(nr)]
        st = np.zeros((1, 1 + nr), np.int64)            # E, xbar
        zc = np.zeros((1, c), np.int64)
        for t, j in enumerate(self.order):
            opts = self.trans[j]; prof = Z.prof[j]; w = int(Z.w[j])
            parts_E = []; parts_z = []
            for (y, xb) in opts:
                dE = w * y * y
                ns = st.copy(); ns[:, 0] += dE; ns[:, 1:] += np.array(xb, np.int64)
                nz = zc + prof[None, :] * y
                ok = ns[:, 0] <= Z.Ecap
                parts_E.append(ns[ok]); parts_z.append(nz[ok])
            st = np.concatenate(parts_E); zc = np.concatenate(parts_z)
            if t in self.comp_end:
                Cq = Z.comps[self.comp_end[t]]
                ok = (zc[:, Cq] == z[Cq][None, :]).all(1)
                st = st[ok]; zc = zc[ok]; zc[:, Cq] = 0
            # dedup
            if len(st) > 1:
                comb = np.concatenate([zc, st], 1)
                comb = np.unique(comb, axis=0)
                zc = comb[:, :c]; st = comb[:, c:]
            if len(st) > self.max_states: raise RuntimeError(f'too many states ({len(st)})')
            if not len(st): break
        if not len(st): return np.zeros(0, np.int64), np.zeros(0, np.int64), 0
        E = st[:, 0]; X = st[:, 1:]
        Pf = S.PZ.paf(X[:, S.orb])[:, S.lreps]
        H = Pf @ S.wv
        keys = np.unique(np.stack([E, Pf[:, 0], H], 1), axis=0)
        return keys[:, 0], keys[:, 2], len(st)

    def alive(self, zq, cache=None):
        """True if some (E_k, xbar_k) in B_k meets sum E_k = T1(0) and the level's PAF equation"""
        Z = self.Z; S = self.S
        cache = {} if cache is None else cache
        B = []
        for z in zq:
            key = tuple(z)
            if key not in cache: cache[key] = self.block_set(z)
            B.append(cache[key])
        if min(len(b[0]) for b in B) == 0: return False, 0
        T0 = Z.T0
        def pairs(a, b):
            Ea, Ha, _ = a; Eb, Hb, _ = b
            E = (Ea[:, None] + Eb[None, :]).ravel(); H = (Ha[:, None] + Hb[None, :]).ravel()
            ok = E <= T0
            return E[ok], H[ok]
        nL = len(B[0][0]) * len(B[1][0]); nR = len(B[2][0]) * len(B[3][0])
        if nL > 4 * 10 ** 8 or nR > 4 * 10 ** 8: raise RuntimeError(f'join too large ({nL}, {nR})')
        EL, HL = pairs(B[0], B[1]); ER, HR = pairs(B[2], B[3])
        kl = HL * np.int64(1000003) + EL
        kr = (np.int64(S.th) - HR) * np.int64(1000003) + (np.int64(T0) - ER)
        kls = np.sort(kl)
        pos = np.searchsorted(kls, kr)
        pos = np.minimum(pos, len(kls) - 1)
        hits = int((kls[pos] == kr).sum()) if len(kls) else 0
        return hits > 0, hits


def main():
    a = sys.argv[1:]
    cfg = SD.cfg_from_json(open(a[0]).read() if a[0].endswith('.json') else a[0]); outp = a[1]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    nQ = cfg['nQ']
    levels = json.loads(opt('--levels', json.dumps([[]] + [[i] for i in range(nQ)])))
    max_states = int(float(opt('--max-states', '3e7')))
    part, nW = (int(x) for x in opt('--part', '0/1').split('/'))
    if nW > 1: outp = outp.replace('.jsonl', f'.p{part}.jsonl')
    done = set()
    if os.path.exists(outp):
        for l in open(outp):
            r = json.loads(l)
            if r.get('stage') == 'ms': done.add(r['i'])
    t0 = time.time()
    Z = ZOrb(cfg)
    Cz = Z.compressions(); cls, idx = Z.classes(Cz); qd = Z.class_quads(cls)
    cache = {}
    for q in qd: q['n0'] = Z.n0(cls, q['m'], cache)
    qd = [q for q in qd if q['n0'] > 0]
    scr = [AggScreen(Z, k, max_states) for k in levels]
    caches = [dict() for _ in levels]
    fo = open(outp, 'a')
    head = dict(stage='setup', multisets=len(qd), levels=levels, part=f'{part}/{nW}', secs=round(time.time() - t0, 1))
    fo.write(json.dumps(head) + "\n"); print(json.dumps(head), flush=True)
    for mi, q in enumerate(qd):
        if mi % nW != part or mi in done: continue
        zq = [cls[i]['rep'] for i in q['m']]
        rec = dict(stage='ms', i=mi, m=list(q['m']), ordered=q['ordered'], log2N0=round(math.log2(q['n0']), 3), dead_at=None, hits=[])
        try:
            for li, A in enumerate(scr):
                ok, h = A.alive(zq, caches[li]); rec['hits'].append(h)
                if not ok: rec['dead_at'] = levels[li]; break
        except (RuntimeError, MemoryError) as e:
            rec['error'] = repr(e)[:150]
        rec['secs'] = round(time.time() - t0, 1)
        fo.write(json.dumps(rec) + "\n"); fo.flush(); print(json.dumps(rec), flush=True)
    fo.close()


if __name__ == '__main__':
    main()
