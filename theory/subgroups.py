"""Subgroups of Aut(G) = G1 x U for the AbGroup/AutGroup encoding (abgroup.py).

G1 = product of the GL(k, p) parts (small, possibly non-abelian), U = product of unit groups of the
cyclic components (abelian).  Subgroups of G1 x U are enumerated with Goursat's lemma:
    S = {(a, b) : theta(a A0) = b B0},  A0 <| A1 <= G1,  B0 <= B1 <= U,  theta: A1/A0 ~ B1/B0.
Because U is abelian, A1/A0 must be abelian.  A1 is taken up to conjugacy in G1 (conjugate symmetry
groups give equivalent slices: conjugation by an automorphism of G maps solutions to solutions).
Returns subgroups as generator lists of AutGroup elements (tuples of component-element indices).
"""
import itertools

import numpy as np


def _closure(gens, mul, e):
    S = {e}
    fr = [e]
    while fr:
        nf = []
        for s in fr:
            for g in gens:
                t = mul(s, g)
                if t not in S:
                    S.add(t)
                    nf.append(t)
        fr = nf
    return frozenset(S)


# ---------------------------------------------------------------- the abelian part U
class Abelian:
    """additive group Z_{n1} x ... x Z_{nk}; elements are tuples"""

    def __init__(self, ns):
        self.ns = tuple(ns)
        self.e = tuple(0 for _ in ns)

    def mul(self, x, y):
        return tuple((a + b) % n for a, b, n in zip(x, y, self.ns))

    def neg(self, x):
        return tuple((-a) % n for a, n in zip(x, self.ns))

    def elements(self):
        return itertools.product(*[range(n) for n in self.ns])

    def order(self, x):
        o = 1
        for a, n in zip(x, self.ns):
            from math import gcd
            o = np.lcm(o, n // gcd(a, n))
        return int(o)

    def subgroups(self):
        """all subgroups, via Sylow decomposition (each Sylow part enumerated by closure of <=3 elements)"""
        from sympy import factorint
        N = int(np.prod(self.ns))
        primes = sorted(factorint(N)) if N > 1 else []
        sylow_subs = []
        for p in primes:
            # elements of p-power order
            els = [x for x in self.elements() if _is_ppower(self.order(x), p)]
            subs = set()
            subs.add(frozenset([self.e]))
            for a in els:
                subs.add(_closure([a], self.mul, self.e))
            changed = True
            cur = set(subs)
            while changed:
                changed = False
                new = set()
                lst = list(cur)
                for A in lst:
                    for a in els:
                        if a in A:
                            continue
                        B = _closure(list(A) + [a], self.mul, self.e) if len(A) < 64 else _closure(_gens_of(A, self) + [a], self.mul, self.e)
                        if B not in cur and B not in new:
                            new.add(B)
                if new:
                    cur |= new
                    changed = True
            sylow_subs.append(list(cur))
        out = []
        for combo in itertools.product(*sylow_subs) if sylow_subs else [()]:
            gens = []
            for H in combo:
                gens += _gens_of(H, self)
            out.append(_closure(gens, self.mul, self.e) if gens else frozenset([self.e]))
        return list(set(out))


def _is_ppower(o, p):
    while o % p == 0:
        o //= p
    return o == 1


def _gens_of(H, grp):
    """a small generating set of the subgroup H (greedy)"""
    gens = []
    cur = frozenset([grp.e])
    for x in sorted(H, key=lambda t: -grp.order(t) if hasattr(grp, 'order') else 0):
        if x not in cur:
            gens.append(x)
            cur = _closure(gens, grp.mul, grp.e)
            if len(cur) == len(H):
                break
    return gens


# ---------------------------------------------------------------- small (possibly non-abelian) G1
class Finite:
    def __init__(self, table, e):
        self.T = table
        self.e = e
        self.N = table.shape[0]
        self.inv = np.zeros(self.N, dtype=int)
        for a in range(self.N):
            for b in range(self.N):
                if table[a, b] == e:
                    self.inv[a] = b
                    break

    def mul(self, a, b):
        return int(self.T[a, b])

    def subgroups(self):
        subs = {frozenset([self.e])}
        cyc = {_closure([a], self.mul, self.e) for a in range(self.N)}
        subs |= cyc
        cl = list(cyc)
        for A in cl:
            for B in cl:
                subs.add(_closure(list(A) + list(B), self.mul, self.e))
        # joins with a third cyclic group (covers 3-generated subgroups of small groups)
        cur = list(subs)
        for A in cur:
            for B in cl:
                if not B <= A:
                    subs.add(_closure(list(A) + list(B), self.mul, self.e))
        return list(subs)

    def conj(self, H, g):
        return frozenset(self.mul(self.mul(g, h), int(self.inv[g])) for h in H)

    def subgroups_up_to_conj(self):
        seen = set()
        reps = []
        for H in sorted(self.subgroups(), key=len):
            if H in seen:
                continue
            cls = {self.conj(H, g) for g in range(self.N)}
            seen |= cls
            reps.append(H)
        return reps

    def is_normal_in(self, K, H):
        return all(self.conj(K, h) == K for h in H)

    def commutator_subgroup(self, H):
        cs = [self.mul(self.mul(a, b), self.mul(int(self.inv[a]), int(self.inv[b]))) for a in H for b in H]
        return _closure(cs, self.mul, self.e)


def product_table(parts):
    """multiplication table of the direct product of component groups (list of (els, mult, ident))"""
    sizes = [len(p[0]) for p in parts]
    N = int(np.prod(sizes)) if sizes else 1
    idx = list(itertools.product(*[range(s) for s in sizes]))
    pos = {t: i for i, t in enumerate(idx)}
    T = np.zeros((N, N), dtype=np.int32)
    for i, a in enumerate(idx):
        for j, b in enumerate(idx):
            T[i, j] = pos[tuple(int(parts[k][1][a[k], b[k]]) for k in range(len(parts)))]
    e = pos[tuple(p[2] for p in parts)]
    return T, e, idx


def goursat_subgroups(AG, max_q=None):
    """subgroups of Aut(G) for an AutGroup AG; returns list of generator lists (AutGroup tuples)."""
    comps = AG.G.comps
    gl_idx = [i for i, c in enumerate(comps) if c[0] == 'E']
    cy_idx = [i for i, c in enumerate(comps) if c[0] == 'C']
    # G1 = product of GL parts
    T1, e1, idx1 = product_table([AG.parts[i] for i in gl_idx])
    G1 = Finite(T1, e1)
    # U = product of unit groups, written additively via discrete logs
    # each unit group is written as a product of cyclic factors: (gens, orders); for odd prime powers
    # one generator, for 2^e (e >= 3) the pair (-1, 5)
    unit_info = []
    ns = []
    for i in cy_idx:
        els, mult, ident = AG.parts[i]
        q = comps[i][1]
        nU = len(els)
        if q % 2 == 0 and q >= 8:
            gens, ords = [q - 1, 5], [2, nU // 2]
        else:
            gens, ords = None, None
            for u in els:
                o = 1
                x = u
                while x != 1:
                    x = x * u % q
                    o += 1
                if o == nU:
                    gens, ords = [u], [nU]
                    break
            assert gens is not None, f'unit group mod {q}: no generator found'
        unit_info.append((gens, ords, q, els))
        ns += ords
    U = Abelian(ns)

    def to_aut_tuple(a1, b):
        """a1 index in G1, b tuple of exponents (concatenated over components) -> AutGroup element tuple"""
        x = list(AG.ident())
        for k, i in enumerate(gl_idx):
            x[i] = idx1[a1][k]
        pos = 0
        for k, i in enumerate(cy_idx):
            gens, ords, q, els = unit_info[k]
            u = 1
            for g, e in zip(gens, b[pos:pos + len(gens)]):
                u = u * pow(g, e, q) % q
            pos += len(gens)
            x[i] = els.index(u)
        return tuple(x)

    Usubs = U.subgroups()
    out = []
    A1s = G1.subgroups_up_to_conj()
    all_g1 = G1.subgroups()
    for A1 in A1s:
        D = G1.commutator_subgroup(A1)
        # normal subgroups A0 of A1 containing D
        A0s = [K for K in all_g1 if K <= A1 and D <= K and G1.is_normal_in(K, A1)]
        for A0 in A0s:
            # quotient Q = A1/A0: cosets
            cos = []
            seen = set()
            for a in sorted(A1):
                if a in seen:
                    continue
                c = frozenset(G1.mul(a, k) for k in A0)
                seen |= c
                cos.append(c)
            qn = len(cos)
            if max_q and qn > max_q:
                continue
            coset_of = {a: i for i, c in enumerate(cos) for a in c}
            qmul = lambda i, j: coset_of[G1.mul(min(cos[i]), min(cos[j]))]
            qe = coset_of[e1]
            # abstract abelian structure of Q: pick generators greedily
            qgens = []
            cur = frozenset([qe])
            for i in sorted(range(qn), key=lambda i: -len(_closure([i], qmul, qe))):
                if i not in cur:
                    qgens.append(i)
                    cur = _closure(qgens, qmul, qe)
                    if len(cur) == qn:
                        break
            for B1 in Usubs:
                if len(B1) % qn:
                    continue
                for B0 in Usubs:
                    if not B0 <= B1 or len(B1) != qn * len(B0):
                        continue
                    # quotient B1/B0 ; enumerate homomorphic images of qgens (theta must be an isomorphism)
                    bcos = []
                    bseen = set()
                    for b in sorted(B1):
                        if b in bseen:
                            continue
                        c = frozenset(U.mul(b, k) for k in B0)
                        bseen |= c
                        bcos.append(c)
                    bcoset_of = {b: i for i, c in enumerate(bcos) for b in c}
                    bmul = lambda i, j: bcoset_of[U.mul(min(bcos[i]), min(bcos[j]))]
                    be = bcoset_of[U.e]
                    for imgs in itertools.product(range(len(bcos)), repeat=len(qgens)):
                        # theta defined on words in qgens; check well-defined & bijective via closure of pairs
                        pairs = _closure([(qg, im) for qg, im in zip(qgens, imgs)],
                                         lambda x, y: (qmul(x[0], y[0]), bmul(x[1], y[1])), (qe, be))
                        if len(pairs) != qn:
                            continue  # not a function or not injective
                        if len({p[1] for p in pairs}) != qn:
                            continue
                        # generators of S
                        gens = []
                        for k in _gens_of_set(A0, G1):
                            gens.append(to_aut_tuple(k, U.e))
                        for b in _gens_abelian(B0, U):
                            gens.append(to_aut_tuple(e1, b))
                        for qg, im in zip(qgens, imgs):
                            gens.append(to_aut_tuple(min(cos[qg]), min(bcos[im])))
                        out.append(gens)
    return out


def _gens_of_set(H, G1):
    gens = []
    cur = frozenset([G1.e])
    for x in sorted(H):
        if x not in cur:
            gens.append(x)
            cur = _closure(gens, G1.mul, G1.e)
            if len(cur) == len(H):
                break
    return gens


def _gens_abelian(H, U):
    gens = []
    cur = frozenset([U.e])
    for x in sorted(H, key=lambda t: -U.order(t)):
        if x not in cur:
            gens.append(x)
            cur = _closure(gens, U.mul, U.e)
            if len(cur) == len(H):
                break
    return gens
