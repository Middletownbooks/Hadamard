r"""Q-level compression patterns of the Z_q all-character bordered family [M] (Kumar, search/NOTES.md section 18,
https://github.com/vatskishan03/hadamard, Apache-2.0; this code is new).

For n = q p with |K| = p - 1 even, the compressions y_i of a core to Q = G/K are even, y_i = 2 z_i, and the
core condition gives  sum_i PAF_{z_i} = q delta  on Q  (so sum_i |z_i|^2 = q and sum_i (sum z_i)^2 = q).
A pattern is a quadruple (z_1..z_4) of integer vectors on Q with that property.  PAF_z is invariant under
translating, negating and reversing z, so we first list the classes of single blocks by (weight, PAF), then the
multisets of four (weight, PAF) classes with weights summing to q and PAFs summing to q delta, then all quadruples
of block vectors realising them, reduced up to: per-block translation, negation, reversal; permutation of the
blocks; a common automorphism of Q.  Q is given by its cyclic factors (q = 9: Z_9 or Z_3^2).
usage: python3 theory/zq_patterns.py q [mods] [--out FILE]
"""
import sys, os, json, itertools, math, argparse
import numpy as np


class QGroup:
    def __init__(self, mods):
        self.mods = tuple(mods)
        self.n = int(np.prod(mods))
        self.co = np.array(list(itertools.product(*[range(m) for m in mods])), np.int64).reshape(self.n, len(mods))
        self.idx = lambda c: int(np.ravel_multi_index(tuple(np.mod(c, self.mods)), self.mods))
        self.add = np.array([[self.idx(self.co[a] + self.co[b]) for b in range(self.n)] for a in range(self.n)])
        self.neg = np.array([self.idx(-self.co[a]) for a in range(self.n)])
        # automorphisms (as permutations): images of the generators
        auts = []
        gens = [self.idx([int(i == j) for i in range(len(mods))]) for j in range(len(mods))]
        M = np.array(mods)
        for imgs in itertools.product(range(self.n), repeat=len(mods)):
            # e_j -> imgs[j] is well defined iff mods[j] * imgs[j] = 0
            if any(((mods[j] * self.co[imgs[j]]) % M).any() for j in range(len(mods))):
                continue
            perm = np.array([self.idx(sum(int(self.co[g][j]) * self.co[imgs[j]] for j in range(len(mods))))
                             for g in range(self.n)])
            if len(set(perm.tolist())) == self.n:
                auts.append(perm)
        self.auts = auts

    def order(self, g):
        c = self.co[g]; k = 1
        while (k * c % np.array(self.mods)).any():
            k += 1
        return k

    def paf(self, z):
        z = np.asarray(z, np.int64)
        return np.array([int((z * z[self.add[:, d]]).sum()) for d in range(self.n)])


def blocks_by_weight(Qg, wmax, amax=3):
    """all integer vectors z on Q with 0 < |z|^2 <= wmax and z[0] > 0 (translation + negation normalised) -> dict
    weight -> {paf tuple: list of z}"""
    n = Qg.n
    out = {}
    vals = list(range(-amax, amax + 1))

    def rec(pos, rem, cur):
        if pos == n:
            w = wmax - rem
            if w == 0 or cur[0] <= 0:
                return
            z = np.array(cur)
            p = tuple(Qg.paf(z)[1:].tolist())
            out.setdefault(w, {}).setdefault(p, []).append(tuple(cur))
            return
        for v in vals:
            if v * v <= rem:
                if pos == 0 and v <= 0:
                    continue
                cur.append(v); rec(pos + 1, rem - v * v, cur); cur.pop()
    rec(0, wmax, [])
    return out


def canon_quad(Qg, quad):
    """canonical form of a quadruple of block vectors under per-block translation/negation/reversal, block
    permutation and a common automorphism"""
    best = None
    n = Qg.n
    for a in Qg.auts:
        blks = []
        for z in quad:
            z = np.asarray(z)
            za = np.zeros(n, np.int64); za[a] = z            # z∘a^{-1}
            reps = []
            for zz in (za, za[Qg.neg]):
                for s in range(n):
                    zt = zz[Qg.add[:, s]]
                    for sg in (1, -1):
                        reps.append(tuple((sg * zt).tolist()))
            blks.append(max(reps))
        key = tuple(sorted(blks, reverse=True))
        if best is None or key > best:
            best = key
    return best


def patterns(q, mods=None, classes_only=False):
    mods = mods or (q,)
    Qg = QGroup(mods)
    B = blocks_by_weight(Qg, q)
    zero = tuple([0] * (Qg.n - 1))
    classes = [(0, zero)] + [(w, p) for w in sorted(B) for p in B[w]]
    target = np.zeros(Qg.n - 1, np.int64)
    found = []
    # multisets of 4 classes, weights descending, sum q, paf sum 0
    byw = {}
    for i, (w, p) in enumerate(classes):
        byw.setdefault(w, []).append(i)
    P = np.array([c[1] for c in classes], np.int64)
    W = np.array([c[0] for c in classes])
    # meet in the middle: weights w1 >= w2 >= w3 >= w4, pairs (1,2) against pairs (3,4) on the PAF sum
    base = 4 * q + 1
    wv = base ** np.arange(Qg.n - 1, dtype=np.int64)
    enc = lambda V: ((V + 2 * q) * wv).sum(-1)
    def pairs(wa, wb):
        A = np.array(byw.get(wa, []), np.int64); B = np.array(byw.get(wb, []), np.int64)
        if not len(A) or not len(B):
            return np.zeros((0, 2), np.int64), np.zeros(0, np.int64)
        I = np.repeat(A, len(B)); J = np.tile(B, len(A))
        keep = J >= I if wa == wb else np.ones(len(I), bool)
        I, J = I[keep], J[keep]
        S = P[I] + P[J]
        ok = np.abs(S).max(1, initial=0) <= q - wa - wb if S.shape[1] else np.ones(len(I), bool)
        return np.stack([I[ok], J[ok]], 1), enc(S[ok])
    for w1 in range(q, -1, -1):
        for w2 in range(min(w1, q - w1), -1, -1):
            p12, k12 = pairs(w1, w2)
            if not len(k12):
                continue
            for w3 in range(min(w2, q - w1 - w2), -1, -1):
                w4 = q - w1 - w2 - w3
                if w4 > w3 or w4 < 0:
                    continue
                p34, k34 = pairs(w3, w4)
                if not len(k34):
                    continue
                o = np.argsort(k34, kind='stable'); k34s = k34[o]; p34s = p34[o]
                need = enc(-(P[p12[:, 0]] + P[p12[:, 1]]))
                lo = np.searchsorted(k34s, need, 'left'); hi = np.searchsorted(k34s, need, 'right')
                for t in np.flatnonzero(hi > lo):
                    i1, i2 = map(int, p12[t])
                    for u in range(lo[t], hi[t]):
                        i3, i4 = map(int, p34s[u])
                        if w3 == w2 and i3 < i2:
                            continue
                        found.append((i1, i2, i3, i4))
    # expand to block-vector quadruples and canonicalise
    pats = {}
    if classes_only:                     # one representative quadruple per multiset of (weight, PAF) block classes
        out = []
        for cl in found:
            quad = [tuple([0] * Qg.n) if classes[i][0] == 0 else B[classes[i][0]][classes[i][1]][0] for i in cl]
            z = np.array(quad)
            tot = sum(Qg.paf(zz) for zz in z)
            assert tot[0] == q and not tot[1:].any()
            out.append(dict(z=z.tolist(), weights=[int((zz * zz).sum()) for zz in z], rows=[int(zz.sum()) for zz in z],
                            support=[int((zz != 0).sum()) for zz in z],
                            n_vectors=[1 if classes[i][0] == 0 else len(B[classes[i][0]][classes[i][1]]) for i in cl]))
        return Qg, out
    for cl in found:
        lists = []
        for i in cl:
            w, p = classes[i]
            lists.append([tuple([0] * Qg.n)] if w == 0 else B[w][p])
        for quad in itertools.product(*lists):
            key = canon_quad(Qg, quad)
            pats.setdefault(key, quad)
    out = []
    for key in sorted(pats):
        z = np.array(key)
        tot = sum(Qg.paf(zz) for zz in z)
        assert tot[0] == q and not tot[1:].any()
        out.append(dict(z=z.tolist(), weights=[int((zz * zz).sum()) for zz in z], rows=[int(zz.sum()) for zz in z],
                        support=[int((zz != 0).sum()) for zz in z]))
    return Qg, out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('q', type=int); ap.add_argument('mods', nargs='?', default=None)
    ap.add_argument('--out', default=None)
    ap.add_argument('--classes-only', action='store_true')
    a = ap.parse_args()
    mods = tuple(map(int, a.mods.split(','))) if a.mods else (a.q,)
    Qg, pats = patterns(a.q, mods, classes_only=a.classes_only)
    print(json.dumps(dict(q=a.q, mods=mods, n_auts=len(Qg.auts), n_patterns=len(pats), classes_only=a.classes_only)))
    for p in pats:
        print(json.dumps(p))
    if a.out:
        json.dump(dict(q=a.q, mods=list(mods), patterns=pats), open(a.out, 'w'))
