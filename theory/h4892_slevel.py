"""H(4892) stage 2: the S-level maps (compression of the stage-2 PAF vector along a subgroup L of G'), used for
bucketing the 2+2 join.  A[c, l] = number of lags of nonzero class l (h4892_stage2.LAGREP_NZ) in the coset of the
compressed lag class c [P: PAF of the compression along L at hbar = sum of PAF over hbar + L]."""
import numpy as np, h4892_stage2 as T


def smap(Lname):
    # projection G -> G/L, as a key
    def proj(g):
        g1, z = divmod(g, T.Q); a, b = divmod(g1, 2)
        if Lname == 'Z2': return (a, z)
        if Lname == 'Z47': return (b, z)
        if Lname == "G'": return (z,)
    # class of each nonzero lag
    lcls = np.full(T.NG, -1)
    reps = T.LAGREP_NZ
    for g in range(1, T.NG):
        p = T.orbit[g]; lcls[g] = reps.index(T.LAGREP[T.cls[p]])
    # compressed classes: group cosets hbar != 0 by the multiset of classes they contain
    cos = {}
    for g in range(1, T.NG):
        cos.setdefault(proj(g), []).append(g)
    rows = {}
    for key, gs in cos.items():
        r = np.bincount(lcls[gs], minlength=T.NL)
        if all(x == 0 for x in key): continue                     # the zero coset: lag-0 of compression, fixed (energy)
        rows.setdefault(tuple(r), 0); rows[tuple(r)] += 1
    A = np.array(sorted(rows), np.int64)
    return A
