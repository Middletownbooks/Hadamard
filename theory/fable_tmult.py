#!/usr/bin/env python3
"""Formal multiplication of Goethals-Seidel quadruples by a cyclic factor Z_t.

Question.  Let (A,B,C,D) be +-1 sequences on an abelian group G_v with
PSD_A + PSD_B + PSD_C + PSD_D = 4v at every nontrivial character.  We look for
four sequences Z_1..Z_4 on Z_t x G_v, each fibre {tau} x G_v carrying one of
the 16 "letters"  +-A, +-B, +-C, +-D, +-A*, +-B*, +-C*, +-D*  (X*(g) = X(-g)),
such that  sum_j PSD_{Z_j} = 4tv  at every nontrivial character, FOR EVERY
quadruple of a given TYPE.  A type is a set of linear relations among the
letters (e.g. "A and B symmetric": A* = A, B* = B; or "D = C*").

At a character (rho, psi) the transform of Z_j is  sum_tau rho^tau u_{j,tau},
u = +-y_k or +-conj(y_k), y = transforms of A..D at psi.  The requirement for all
rho in mu_t is: for every shift d != 0 the formal polynomial
    P_d = sum_j sum_tau u_{j,tau+d} * conj(u_{j,tau})
in the commuting symbols y_k, conj(y_k) vanishes on the type, and at d = 0
equals t*|y|^2 (letter counts).

This script enumerates all assignments exhaustively (meet in the middle over
blocks) for small t.  [M] results are exact: the identity is formal.
"""
import sys, itertools, collections
import numpy as np

def symbols_for_type(typ):
    """Return map letter index (k, s) -> canonical symbol id, and names.
    k in 0..3 (A..D), s = 0 plain, 1 conjugate."""
    # canonical symbols: list of names; conj map on canonical symbols
    if typ == 'generic':
        sym = {(k, s): 2 * k + s for k in range(4) for s in range(2)}
        conj = {2 * k + s: 2 * k + 1 - s for k in range(4) for s in range(2)}
        norm_groups = [[0, 1], [2, 3], [4, 5], [6, 7]]
        letters = [(k, s) for k in range(4) for s in range(2)]
    elif typ == 'AB_sym':        # A* = A, B* = B
        sym = {(0, 0): 0, (0, 1): 0, (1, 0): 1, (1, 1): 1, (2, 0): 2, (2, 1): 3, (3, 0): 4, (3, 1): 5}
        conj = {0: 0, 1: 1, 2: 3, 3: 2, 4: 5, 5: 4}
        letters = [(0, 0), (1, 0), (2, 0), (2, 1), (3, 0), (3, 1)]
    elif typ == 'A_sym':
        sym = {(0, 0): 0, (0, 1): 0, (1, 0): 1, (1, 1): 2, (2, 0): 3, (2, 1): 4, (3, 0): 5, (3, 1): 6}
        conj = {0: 0, 1: 2, 2: 1, 3: 4, 4: 3, 5: 6, 6: 5}
        letters = [(0, 0), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (3, 1)]
    elif typ == 'ABC_sym':
        sym = {(0, 0): 0, (0, 1): 0, (1, 0): 1, (1, 1): 1, (2, 0): 2, (2, 1): 2, (3, 0): 3, (3, 1): 4}
        conj = {0: 0, 1: 1, 2: 2, 3: 4, 4: 3}
        letters = [(0, 0), (1, 0), (2, 0), (3, 0), (3, 1)]
    elif typ == 'all_sym':
        sym = {(k, s): k for k in range(4) for s in range(2)}
        conj = {k: k for k in range(4)}
        letters = [(k, 0) for k in range(4)]
    elif typ == 'D_eq_Cstar':    # D = C*: letters A, A*, B, B*, C, C*
        sym = {(0, 0): 0, (0, 1): 1, (1, 0): 2, (1, 1): 3, (2, 0): 4, (2, 1): 5, (3, 0): 5, (3, 1): 4}
        conj = {0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4}
        letters = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]
    elif typ == 'A_sym_D_eq_Cstar':
        sym = {(0, 0): 0, (0, 1): 0, (1, 0): 1, (1, 1): 2, (2, 0): 3, (2, 1): 4, (3, 0): 4, (3, 1): 3}
        conj = {0: 0, 1: 2, 2: 1, 3: 4, 4: 3}
        letters = [(0, 0), (1, 0), (1, 1), (2, 0), (2, 1)]
    elif typ == 'B_eq_Astar_D_eq_Cstar':
        sym = {(0, 0): 0, (0, 1): 1, (1, 0): 1, (1, 1): 0, (2, 0): 2, (2, 1): 3, (3, 0): 3, (3, 1): 2}
        conj = {0: 1, 1: 0, 2: 3, 3: 2}
        letters = [(0, 0), (0, 1), (2, 0), (2, 1)]
    else:
        raise ValueError(typ)
    return sym, conj, letters

def run(t, typ, verbose=True, want=3):
    sym, conj, letters = symbols_for_type(typ)
    nsym = len(conj)
    # monomial index for unordered pairs of canonical symbols
    mono = {}
    for a in range(nsym):
        for b in range(a, nsym):
            mono[(a, b)] = len(mono)
    nm = len(mono)
    # "norm" monomials y_k * conj(y_k) for each original index k
    half = (t - 1) // 2
    L = [(e, l) for e in (1, -1) for l in range(len(letters))]
    # letter weight vector for d = 0: counts per original index k, with multiplicities
    # of the type (e.g. D = C* means C-letters count for both C and D).
    def block_vec(cfg):
        # cfg: tuple of (sign, letter idx) length t
        v = np.zeros((half, nm), dtype=np.int16)
        for d in range(1, half + 1):
            for tau in range(t):
                e1, l1 = cfg[(tau + d) % t]
                e0, l0 = cfg[tau]
                a = sym[letters[l1]]
                b = conj[sym[letters[l0]]]
                if a > b:
                    a, b = b, a
                v[d - 1, mono[(a, b)]] += e1 * e0
        return v.reshape(-1)
    # also need P_{-d} = 0: formal conjugate of P_d; if P_d = 0 formally then conj is 0. fine.
    # count vector: number of letters with canonical symbol class {s, conj s}
    cls = {}
    for s in range(nsym):
        key = min(s, conj[s])
        cls.setdefault(key, len(cls))
    ncls = len(cls)
    def count_vec(cfg):
        c = np.zeros(ncls, dtype=np.int16)
        for e, l in cfg:
            s = sym[letters[l]]
            c[cls[min(s, conj[s])]] += 1
        return c
    # target counts: each original index k appears t times; class weights
    target = np.zeros(ncls, dtype=np.int16)
    for k in range(4):
        s = sym[(k, 0)]
        target[cls[min(s, conj[s])]] += t
    # enumerate block configs up to global sign and cyclic shift (canonical: lexicographically least)
    cfgs = []
    seen = set()
    for cfg in itertools.product(L, repeat=t):
        # canonical under sign and shift
        best = None
        for sh in range(t):
            c1 = tuple(cfg[(i + sh) % t] for i in range(t))
            for sg in (1, -1):
                c2 = tuple((e * sg, l) for e, l in c1)
                if best is None or c2 < best:
                    best = c2
        if best in seen:
            continue
        seen.add(best)
        cfgs.append(best)
    if verbose:
        print(f"t={t} type={typ}: {len(L)} signed letters, {len(cfgs)} block classes, {nm} monomials x {half} shifts", flush=True)
    V = np.array([block_vec(c) for c in cfgs])
    C = np.array([count_vec(c) for c in cfgs])
    n = len(cfgs)
    # pairs i <= j : key = (V_i + V_j, C_i + C_j)
    full = np.concatenate([V, C], axis=1).astype(np.int16)
    tgt = np.concatenate([np.zeros(V.shape[1], dtype=np.int16), target])
    d = collections.defaultdict(list)
    for i in range(n):
        S = full[i] + full[i:]
        for jj in range(S.shape[0]):
            d[S[jj].tobytes()].append((i, i + jj))
    sols = []
    nsol = 0
    for key, prs in d.items():
        s = np.frombuffer(key, dtype=np.int16)
        comp = (tgt - s).astype(np.int16).tobytes()
        if comp in d and key <= comp:
            for p in prs:
                for q in d[comp]:
                    if key != comp or p <= q:
                        nsol += 1
                        if len(sols) < want:
                            sols.append((p, q))
    if verbose:
        print(f"  solutions (unordered pairs of block-pairs, classes mod sign/shift): {nsol}", flush=True)
        names = 'ABCD'
        for (p, q) in sols:
            for i in p + q:
                print('   ', ' '.join(('+' if e > 0 else '-') + names[letters[l][0]] + ('*' if letters[l][1] else ' ') for e, l in cfgs[i]))
            print('    --')
    return nsol, sols, cfgs, letters

if __name__ == '__main__':
    t = int(sys.argv[1])
    for typ in sys.argv[2:]:
        run(t, typ)
