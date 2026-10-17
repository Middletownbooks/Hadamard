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
        self._ccache = {}
        self._gkey = []
        for i, C in enumerate(Z.corbs):
            touched = set()
            for j in C:
                for (y, xb) in self.trans[j]: touched |= {k for k, v in enumerate(xb) if v}
            self._gkey.append(tuple(sorted(touched)))

    def comp_set(self, i, zc_target):
        """(E, xbar) pairs over the letter assignments and lifts of component i's orbits with compression zc_target
        (rows of an int64 array, deduplicated)"""
        Z = self.Z; nr = self.nr; C = Z.comps[i]
        key = (i, tuple(zc_target))
        if key in self._ccache: return self._ccache[key]
        st = np.zeros((1, len(C) + 1 + nr), np.int64)           # zc (component cosets), E, xbar
        orbs = Z.corbs[i]
        # remaining-range pruning on the component compression
        lo = np.zeros((len(orbs) + 1, len(C)), np.int64); hi = np.zeros((len(orbs) + 1, len(C)), np.int64)
        for t in range(len(orbs) - 1, -1, -1):
            pj = Z.prof[orbs[t], C]; L = Z.letters[orbs[t]]
            lo[t] = lo[t + 1] + np.minimum(pj * min(L), pj * max(L)); hi[t] = hi[t + 1] + np.maximum(pj * min(L), pj * max(L))
        tgt = np.asarray(zc_target, np.int64)
        for t, j in enumerate(orbs):
            pj = Z.prof[j, C]; w = int(Z.w[j]); parts = []
            for (y, xb) in self.trans[j]:
                d = np.concatenate([pj * y, [w * y * y], np.array(xb, np.int64)])
                ns = st + d[None, :]
                rem = tgt[None, :] - ns[:, :len(C)]
                ok = (ns[:, len(C)] <= Z.Ecap) & (rem >= lo[t + 1]).all(1) & (rem <= hi[t + 1]).all(1)
                parts.append(ns[ok])
            st = np.unique(np.concatenate(parts), axis=0)
            if len(st) > self.max_states: raise RuntimeError(f'too many states in a component ({len(st)})')
        st = st[(st[:, :len(C)] == tgt[None, :]).all(1)][:, len(C):]
        st = np.unique(st, axis=0)
        self._ccache[key] = st
        return st

    def block_set(self, z, max_product=6 * 10 ** 8):
        """distinct (E, signature-hash) pairs of block lifts with compression z (all letter assignments, all lifts):
        per coset component an isolated DP; components with the same image in Qbar are combined first (deduplicated
        Minkowski sums: their xbar coordinates coincide), then the groups (disjoint coordinates) are multiplied out in
        chunks and hashed directly"""
        Z = self.Z; S = self.S; nr = self.nr
        z = np.asarray(z, np.int64)
        sets = [self.comp_set(i, tuple(int(z[q]) for q in C)) for i, C in enumerate(Z.comps)]
        if any(len(x) == 0 for x in sets): return np.zeros(0, np.int64), np.zeros(0, np.int64), 0
        # group by the xbar coordinates a component can touch
        groups = {}
        for i, x in enumerate(sets):
            nz = tuple(np.nonzero(np.abs(x[:, 1:]).sum(0))[0].tolist()) if len(x) else ()
            groups.setdefault(self._gkey[i], []).append(i)
        G = []
        for gk, comps in sorted(groups.items()):
            acc = np.zeros((1, 1 + nr), np.int64)
            for i in comps:
                x = sets[i]; out = []
                CH = max(1, 2_000_000 // max(1, len(x)))
                for c0 in range(0, len(acc), CH):
                    s_ = (acc[c0:c0 + CH][:, None, :] + x[None, :, :]).reshape(-1, 1 + nr)
                    out.append(np.unique(s_[s_[:, 0] <= Z.Ecap], axis=0))
                acc = np.unique(np.concatenate(out), axis=0)
                if len(acc) > self.max_states: raise RuntimeError(f'too many group states ({len(acc)})')
            G.append(acc)
        tot = math.prod(len(g) for g in G)
        if tot > max_product: raise RuntimeError(f'group product too large ({tot})')
        # multiply the groups out in chunks: E and xbar add (coordinates of different groups are disjoint)
        keys = []
        def rec(gi, accE, accX):
            if gi == len(G):
                ok = accE <= Z.Ecap
                E = accE[ok]; X = accX[ok]
                if len(E):
                    Pf = S.PZ.paf(X[:, S.orb])[:, S.lreps]; H = Pf @ S.wv
                    keys.append(np.unique(np.stack([E, H], 1), axis=0))
                return
            g = G[gi]
            CH = max(1, 1_000_000 // max(1, len(g)))
            for c0 in range(0, len(accE), CH):
                e = (accE[c0:c0 + CH][:, None] + g[None, :, 0]).ravel()
                X = (accX[c0:c0 + CH][:, None, :] + g[None, :, 1:]).reshape(-1, nr)
                rec(gi + 1, e, X)
        rec(0, np.zeros(1, np.int64), np.zeros((1, nr), np.int64))
        if not keys: return np.zeros(0, np.int64), np.zeros(0, np.int64), 0
        out = np.unique(np.concatenate(keys), axis=0)
        return out[:, 0], out[:, 1], tot

    def alive(self, zq, cache=None, early=True, max_bucket=6 * 10 ** 7):
        """True if some (E_k, xbar_k) in B_k meets sum E_k = T1(0) and the level's PAF equation.  Energy-bucketed 2+2
        join: for every split EL + ER = T1(0), the left pairs with energy EL against the right pairs with energy ER;
        with early=True it stops at the first match (hits then counts matches found up to that point)."""
        Z = self.Z; S = self.S
        cache = {} if cache is None else cache
        B = []
        for z in zq:
            key = tuple(z)
            if key not in cache: cache[key] = self.block_set(z)
            B.append(cache[key])
        if min(len(b[0]) for b in B) == 0: return False, 0
        T0 = Z.T0
        def groups(b):
            E, H, _ = b; o = np.argsort(E, kind='stable'); Es = E[o]
            u, st = np.unique(Es, return_index=True); en = np.r_[st[1:], len(Es)]
            return {int(a): H[o][i:j] for a, i, j in zip(u, st, en)}
        G = [groups(b) for b in B]
        def pair_index(g0, g1):
            idx = {}
            for e0 in g0:
                for e1 in g1:
                    if e0 + e1 <= T0: idx.setdefault(e0 + e1, []).append((e0, e1))
            return idx
        L = pair_index(G[0], G[1]); R = pair_index(G[2], G[3])
        th = np.int64(S.th); hits = 0
        for EL, lst in L.items():
            rl = R.get(T0 - EL)
            if not rl: continue
            nl = sum(len(G[0][a]) * len(G[1][b]) for a, b in lst); nr = sum(len(G[2][a]) * len(G[3][b]) for a, b in rl)
            if min(nl, nr) > max_bucket: raise RuntimeError(f'energy bucket too large ({nl}, {nr})')
            # build the smaller side as the sorted table, stream the other in chunks
            if nl <= nr:
                tab = np.sort(np.concatenate([(G[0][a][:, None] + G[1][b][None, :]).ravel() for a, b in lst]))
                for a, b in rl:
                    h2 = G[2][a]
                    for c0 in range(0, len(h2), max(1, 4_000_000 // max(1, len(G[3][b])))):
                        need = (th - (h2[c0:c0 + 4_000_000 // max(1, len(G[3][b]))][:, None] + G[3][b][None, :])).ravel()
                        pos = np.minimum(np.searchsorted(tab, need), len(tab) - 1)
                        hits += int((tab[pos] == need).sum())
                        if early and hits: return True, hits
            else:
                tab = np.sort(np.concatenate([(th - (G[2][a][:, None] + G[3][b][None, :])).ravel() for a, b in rl]))
                for a, b in lst:
                    h0 = G[0][a]
                    for c0 in range(0, len(h0), max(1, 4_000_000 // max(1, len(G[1][b])))):
                        need = (h0[c0:c0 + 4_000_000 // max(1, len(G[1][b]))][:, None] + G[1][b][None, :]).ravel()
                        pos = np.minimum(np.searchsorted(tab, need), len(tab) - 1)
                        hits += int((tab[pos] == need).sum())
                        if early and hits: return True, hits
        return hits > 0, hits


def main():
    a = sys.argv[1:]
    cfg = SD.cfg_from_json(open(a[0]).read() if a[0].endswith('.json') else a[0]); outp = a[1]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    nQ = cfg['nQ']
    levels = json.loads(opt('--levels', json.dumps([[]] + [[i] for i in range(nQ)])))
    max_states = int(float(opt('--max-states', '3e7')))
    part, nW = (int(x) for x in opt('--part', '0/1').split('/'))
    run_config(cfg, outp, levels, max_states, part, nW)


def run_config(cfg, outp, levels, max_states=3 * 10 ** 7, part=0, nW=1):
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


class CharAggScreen:
    """Z_q level for one character chi of Q (real, or real on U-invariant compressions) [P, M]:
    xbar_k^chi(t) = sum_{g' in G'} chi(q(g')) x_k(g', t) on Z_q satisfies sum_k PAF_{Z_q}(xbar_k^chi) = T^chi,
    T^chi(tau) = sum_{g in G, z(g) = tau} chi(q(g)) T(g) (report/paley6_answer.md 3.1).  Orbit j adds
    (sum_q chi(q) n_{j,q}) times its (per-point) fibre pattern, so per block the set of (E, xbar^chi) is a
    Minkowski sum over the orbits; aggregated over all letter assignments with the block's compression."""
    def __init__(self, Z, chi_vals, max_states=3 * 10 ** 7):
        self.Z = Z; self.S = Screen(Z.cfg, keep=[], sl=Z.sl); S = self.S
        self.max_states = max_states
        chi = np.asarray(chi_vals, complex)
        wts = [complex((chi * Z.prof[j]).sum()) for j in range(Z.L1)]
        if any(abs(w.imag) > 1e-9 for w in wts): raise ValueError('character not real on these orbits')
        self.wt = [int(round(w.real)) for w in wts]
        # target: T^chi on Z_q
        G = Z.sl.G; nQ = Z.cfg['nQ']
        qg = Z.Q.idx(G.co[:, :nQ]); chig = chi[qg].real
        Tq = np.zeros(S.PZ.n, float); np.add.at(Tq, G.co[:, -1], chig * np.asarray(Z.sl.T, float))
        self.T = np.rint(Tq).astype(np.int64)[S.lreps]
        self.th = int(self.T.astype(np.int64) @ S.wv)
        self.nr = len(S.oreps)
        self.trans = []
        for j in range(Z.L1):
            nj = int(Z.prof[j].sum()); opts = []
            for y in Z.letters[j]:
                for xb in S.patterns(j, y):
                    per_point = np.array(xb, np.int64) // nj
                    assert (per_point * nj == np.array(xb)).all()
                    opts.append((int(y), tuple(int(v) for v in self.wt[j] * per_point)))
            self.trans.append(opts)
        self._ccache = {}

    comp_set = AggScreen.comp_set

    def block_set(self, z):
        Z = self.Z; S = self.S; nr = self.nr
        z = np.asarray(z, np.int64)
        sets = [self.comp_set(i, tuple(int(z[q]) for q in C)) for i, C in enumerate(Z.comps)]
        if any(len(x) == 0 for x in sets): return np.zeros(0, np.int64), np.zeros(0, np.int64), 0
        acc = np.zeros((1, 1 + nr), np.int64)
        for x in sets:
            out = []
            CH = max(1, 2_000_000 // max(1, len(x)))
            for c0 in range(0, len(acc), CH):
                s_ = (acc[c0:c0 + CH][:, None, :] + x[None, :, :]).reshape(-1, 1 + nr)
                s_ = s_[s_[:, 0] <= Z.Ecap]
                out.append(np.unique(s_, axis=0))
            acc = np.unique(np.concatenate(out), axis=0)
            if len(acc) > self.max_states: raise RuntimeError(f'too many block states ({len(acc)})')
        E = acc[:, 0]; X = acc[:, 1:]
        Pf = S.PZ.paf(X[:, S.orb])[:, S.lreps]
        H = Pf @ S.wv
        out = np.unique(np.stack([E, H], 1), axis=0)
        return out[:, 0], out[:, 1], len(acc)

    def alive(self, zq, cache=None, early=True):
        save = (self.S.th,)
        self.S.th = self.th
        try:
            return AggScreen.alive(self, zq, cache, early)
        finally:
            self.S.th = save[0]


def real_characters(Z):
    """the characters of Q that are real on the U-invariant compressions (one per complex-conjugate pair)"""
    Q = Z.Q; mods = Q.mods; out = []; seen = set()
    for ks in itertools.product(*[range(m) for m in mods]):
        vals = np.exp(2j * np.pi * (Q.co * np.array(ks)[None, :] / np.array(mods)[None, :]).sum(1))
        wts = [complex((vals * Z.prof[j]).sum()) for j in range(Z.L1)]
        if any(abs(w.imag) > 1e-9 for w in wts): continue
        key = tuple(int(round(w.real)) for w in wts)
        if key in seen: continue
        seen.add(key); out.append((tuple(ks), vals))
    return out
