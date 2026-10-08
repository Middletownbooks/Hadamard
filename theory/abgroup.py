"""Finite abelian groups G = Z_{m1} x ... x Z_{mk}, automorphisms, symmetry subgroups, orbits and classes.

Elements are indexed 0..n-1 in mixed radix (last coordinate fastest).  Automorphisms are stored as
index permutations.  Components:
  ('E', p, k)  elementary abelian Z_p^k   (automorphisms GL(k, p), given as k x k matrices mod p)
  ('C', q)     cyclic Z_q (q prime power) (automorphisms = units mod q)
A symmetry group S is given by generators; orbits and lag/frequency classes are connected components,
so S never has to be enumerated.

Characters are indexed by G itself through <x, xi> = sum_j x_j xi_j / m_j (mod 1).  For S-invariant
X one has hat X(alpha^T xi) = hat X(xi), so frequency classes are orbits of {alpha^T} and -1.
"""
import itertools
from math import gcd

import numpy as np


class AbGroup:
    def __init__(self, comps):
        self.comps = comps
        mods = []
        for c in comps:
            if c[0] == 'E':
                mods += [c[1]] * c[2]
            else:
                mods.append(c[1])
        self.mods = tuple(mods)
        self.n = int(np.prod(mods))
        self.coords = np.array(np.unravel_index(np.arange(self.n), self.mods)).T  # n x k
        self.k = len(mods)

    def idx(self, coords):
        c = np.asarray(coords) % np.array(self.mods)
        return np.ravel_multi_index(tuple(np.moveaxis(c, -1, 0)), self.mods)

    def neg_perm(self):
        return self.idx(-self.coords)

    def add_perm(self, g):
        return self.idx(self.coords + self.coords[g])

    def __repr__(self):
        return 'x'.join(f'Z{m}' for m in self.mods)

    # ---- automorphisms: an automorphism is a list, one entry per component:
    #      ('E') a k x k integer matrix mod p ;  ('C') a unit mod q
    def aut_perm(self, aut, transpose=False):
        out = np.zeros_like(self.coords)
        col = 0
        for c, a in zip(self.comps, aut):
            if c[0] == 'E':
                p, k = c[1], c[2]
                M = np.array(a) % p
                if transpose:
                    M = M.T
                out[:, col:col + k] = (self.coords[:, col:col + k] @ M.T) % p
                col += k
            else:
                q = c[1]
                out[:, col] = (self.coords[:, col] * a) % q
                col += 1
        return self.idx(out)

    def identity_aut(self):
        return [np.eye(c[2], dtype=int) if c[0] == 'E' else 1 for c in self.comps]


# ------------------------------------------------------------ component automorphism groups
def gl_elements(k, p):
    """all invertible k x k matrices over GF(p) (small k, p)"""
    els = []
    for entries in itertools.product(range(p), repeat=k * k):
        M = np.array(entries).reshape(k, k)
        if round(np.linalg.det(M)) % p != 0 and _rank_mod_p(M, p) == k:
            els.append(M)
    return els


def _rank_mod_p(M, p):
    A = M.copy() % p
    r = 0
    rows, cols = A.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i, c] % p:
                piv = i
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), -1, p)
        A[r] = (A[r] * inv) % p
        for i in range(rows):
            if i != r and A[i, c]:
                A[i] = (A[i] - A[i, c] * A[r]) % p
        r += 1
    return r


def units(q):
    return [a for a in range(1, q) if gcd(a, q) == 1]


def mat_key(M):
    return tuple(int(x) for x in np.asarray(M).ravel())


def component_group(c):
    """(elements, mult, inverse, identity_index, encode) for an automorphism group of one component"""
    if c[0] == 'E':
        p, k = c[1], c[2]
        els = gl_elements(k, p)
        key = {mat_key(M): i for i, M in enumerate(els)}
        nE = len(els)
        mult = np.zeros((nE, nE), dtype=np.int32)
        for i, A in enumerate(els):
            for j, B in enumerate(els):
                mult[i, j] = key[mat_key((A @ B) % p)]
        ident = key[mat_key(np.eye(k, dtype=int))]
        return els, mult, ident
    q = c[1]
    els = units(q)
    key = {u: i for i, u in enumerate(els)}
    nE = len(els)
    mult = np.zeros((nE, nE), dtype=np.int32)
    for i, a in enumerate(els):
        for j, b in enumerate(els):
            mult[i, j] = key[a * b % q]
    return els, mult, key[1]


class AutGroup:
    """Aut(G) as a direct product of the component automorphism groups (valid when components have
    pairwise coprime orders, or are the full Sylow subgroup for their prime)."""

    def __init__(self, G):
        self.G = G
        self.parts = [component_group(c) for c in G.comps]
        self.sizes = [len(p[0]) for p in self.parts]
        self.order = int(np.prod(self.sizes))

    def mul(self, x, y):
        return tuple(int(self.parts[i][1][x[i], y[i]]) for i in range(len(x)))

    def ident(self):
        return tuple(p[2] for p in self.parts)

    def elem_order(self, x):
        e = self.ident()
        y = x
        k = 1
        while y != e:
            y = self.mul(y, x)
            k += 1
        return k

    def to_aut(self, x):
        return [self.parts[i][0][x[i]] for i in range(len(x))]

    def closure(self, gens):
        e = self.ident()
        S = {e}
        frontier = [e]
        while frontier:
            nf = []
            for s in frontier:
                for g in gens:
                    t = self.mul(s, g)
                    if t not in S:
                        S.add(t)
                        nf.append(t)
            frontier = nf
        return frozenset(S)

    def all_elements(self):
        return itertools.product(*[range(s) for s in self.sizes])


# ------------------------------------------------------------ orbits and classes
def components(n, perms):
    """connected components of the graph x ~ perm(x) over the given permutations"""
    parent = np.arange(n)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for P in perms:
        for x in range(n):
            a, b = find(x), find(int(P[x]))
            if a != b:
                parent[max(a, b)] = min(a, b)
    lab = np.array([find(x) for x in range(n)])
    reps, inv = np.unique(lab, return_inverse=True)
    return inv  # orbit label per element, labels 0..r-1 ordered by smallest member


def slice_structure(G, auts):
    """orbits of <auts> on G (blocks), lag classes (orbits of <auts,-1> on G\\{0}),
    frequency classes (orbits of <auts^T,-1> on G^\\{0})"""
    perms = [G.aut_perm(a) for a in auts]
    permsT = [G.aut_perm(a, transpose=True) for a in auts]
    neg = G.neg_perm()
    orb = components(G.n, perms)
    lag = components(G.n, perms + [neg])
    frq = components(G.n, permsT + [neg])
    return orb, lag, frq
