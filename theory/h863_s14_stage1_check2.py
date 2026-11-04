r"""[M] Second, independent decision of stage 1 of the S14 sibling at v = 863 (see theory/h863_s14_stage1_exact.py; the
set-up of the S7/S14 hosts is by OpenAI Codex, the construction Kumar / EKLO).  Shares no code with the first method:
  - per block, the tau-invariant y are enumerated by a recursion over the 11 tau-orbits (cells listed explicitly),
    pruned by the block's three Q-sums (exact reachability of the remaining orbits of each Q-sum's orbit group) and
    an energy cap (cap_i = 2856 - sum_{j != i} min_j, the minima by an exact DP over each group's orbits);
  - the signature is the psi-form of the H-level: for psi in {00, 01, 11} the Z5 sequences Y(u) = sum_q psi(q) y(q, u)
    and their PAF at u0 = 0, 1, 2 (9 values; psi10 is the mirror of psi01).  Targets: psi00 (3452, 0, 0),
    psi01 and psi11 (2764, -688, -688);
  - the decision is a CP-SAT model: one index variable per block, a table constraint tying it to its 9 signature
    values, and the 9 sum equations.  INFEASIBLE from CP-SAT on this finite model is a proof for the listed blocks.
usage: python3 theory/h863_s14_stage1_check2.py [first] [last] [secs] [out.jsonl]"""
import sys, os, json, time
from ortools.sat.python import cp_model

VALS = [1, -1, 13, -13, 15, -15, 27, -27, 29, -29, 41, -41, 43, -43]
# tau-orbits as explicit cell lists (cell = (q, u), q in {'00','01','10','11'})
ORBITS = [[('00', 0)], [('00', 1), ('00', 4)], [('00', 2), ('00', 3)],
          [('11', 0)], [('11', 1), ('11', 4)], [('11', 2), ('11', 3)]] + \
         [[('01', u), ('10', (-u) % 5)] for u in range(5)]
GROUP = {'00': 0, '01': 1, '10': 1, '11': 2}            # Q-sum index: a = z(00), b = z(01) = z(10), c = z(11)
PSI = {'00': {'00': 1, '01': 1, '10': 1, '11': 1}, '01': {'00': 1, '01': -1, '10': 1, '11': -1},
       '11': {'00': 1, '01': -1, '10': -1, '11': 1}}
TARGET = [3452, 0, 0, 2764, -688, -688, 2764, -688, -688]


def contrib(o):
    """per orbit: (multiplicity in each Q-sum, energy multiplier = number of cells)"""
    m = [0, 0, 0]
    for (q, u) in ORBITS[o]:
        if q != '10': m[GROUP[q]] += 1      # z(01) counts the (01, u) cells; the (10, -u) cells repeat them
    return m, len(ORBITS[o])


CONTRIB = [contrib(o) for o in range(11)]


# The three Q-sums involve disjoint orbit groups: a <- orbits 0-2, c <- orbits 3-5, b <- orbits 6-10.  Reachability and
# minimal energy therefore factor over the groups (exact).
GROUPS = [(0, (0, 1, 2)), (2, (3, 4, 5)), (1, (6, 7, 8, 9, 10))]   # (Q-sum index, orbits)


def _group_tables():
    """for each orbit o: the set of sums reachable by the orbits of its group after o (coordinate of that group), and
    for each group: the exact minimal energy of every reachable group sum"""
    RG = {}; MINE = {}
    for g, orbs in GROUPS:
        reach = {orbs[-1] + 1: {0}}
        for k in range(len(orbs) - 1, -1, -1):
            o = orbs[k]; m = CONTRIB[o][0][g]
            reach[o] = {r + m * v for r in reach[o + 1] for v in VALS}
        for o in orbs: RG[o] = (g, reach[o + 1])
        best = {0: 0}
        for o in orbs:
            m, n = CONTRIB[o][0][g], CONTRIB[o][1]
            nb = {}
            for t, e in best.items():
                for v in VALS:
                    tt = t + m * v; ee = e + n * (v * v - 1)
                    if ee < nb.get(tt, 10 ** 9): nb[tt] = ee
            best = nb
        MINE[g] = best
    return RG, MINE


RG, MINE = _group_tables()


def min_energy(z):
    """exact min of sum_cells (y^2 - 1) given the three Q-sums (factorised over the orbit groups)"""
    tot = 0
    for g in range(3):
        if z[g] not in MINE[g]: return None
        tot += MINE[g][z[g]]
    return tot


def enumerate_block(z, cap):
    out = []
    vals = [0] * 11
    def rec(o, rem, e):
        if o == 11:
            if rem == (0, 0, 0): out.append(list(vals))
            return
        m, n = CONTRIB[o]
        for v in VALS:
            ee = e + n * (v * v - 1)
            if ee > cap: continue
            r2 = (rem[0] - m[0] * v, rem[1] - m[1] * v, rem[2] - m[2] * v)
            g, rs = RG[o]
            if r2[g] not in rs: continue
            vals[o] = v; rec(o + 1, r2, ee)
    rec(0, tuple(z), 0)
    return out


def signature(vals):
    y = {}
    for o, cl in enumerate(ORBITS):
        for c in cl: y[c] = vals[o]
    sig = []
    for p in ('00', '01', '11'):
        Y = [sum(PSI[p][q] * y[(q, u)] for q in ('00', '01', '10', '11')) for u in range(5)]
        for u0 in (0, 1, 2):
            sig.append(sum(Y[u] * Y[(u + u0) % 5] for u in range(5)))
    return sig


def decide(z, secs, target=TARGET, etot=2856):
    mins = [min_energy(b) for b in z]
    if any(m is None for m in mins) or sum(mins) > etot:
        return dict(status='EMPTY (Q-sums or energy)', mins=mins, sizes=None)
    lists = [enumerate_block(b, etot - (sum(mins) - mins[i])) for i, b in enumerate(z)]
    sizes = [len(L) for L in lists]
    if min(sizes) == 0: return dict(status='EMPTY (a block list is empty)', mins=mins, sizes=sizes)
    sigs = [[signature(v) for v in L] for L in lists]
    M = cp_model.CpModel()
    S = []
    for i in range(4):
        k = M.NewIntVar(0, sizes[i] - 1, 'k%d' % i)
        lo = [min(s[c] for s in sigs[i]) for c in range(9)]; hi = [max(s[c] for s in sigs[i]) for c in range(9)]
        sv = [M.NewIntVar(lo[c], hi[c], '') for c in range(9)]
        M.AddAllowedAssignments([k] + sv, [[n] + s for n, s in enumerate(sigs[i])])
        S.append((k, sv))
    for c in range(9):
        M.Add(sum(S[i][1][c] for i in range(4)) == target[c])
    sol = cp_model.CpSolver(); sol.parameters.max_time_in_seconds = secs; sol.parameters.num_workers = 4
    st = sol.StatusName(sol.Solve(M))
    out = dict(status=st, mins=mins, sizes=sizes)
    if st in ('OPTIMAL', 'FEASIBLE'):
        out['y_orbits'] = [lists[i][sol.Value(S[i][0])] for i in range(4)]
    return out


if __name__ == '__main__':
    Z = [tuple(tuple(b) for b in z) for z in json.load(open('data/h863_s14_qsolutions.json'))['z']]
    first = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    last = int(sys.argv[2]) if len(sys.argv) > 2 else len(Z)
    secs = float(sys.argv[3]) if len(sys.argv) > 3 else 600
    outp = sys.argv[4] if len(sys.argv) > 4 else None
    for n in range(first, min(last, len(Z))):
        t = time.time(); r = decide(Z[n], secs)
        r.update(index=n, z=Z[n], rows=[a + 2 * b + c for a, b, c in Z[n]], secs=round(time.time() - t, 1))
        s = json.dumps(r); print(s, flush=True)
        if outp:
            with open(outp, 'a') as fo: fo.write(s + '\n')
