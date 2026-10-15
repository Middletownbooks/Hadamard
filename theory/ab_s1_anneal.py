r"""Stage-1 solutions of a GP slice with prescribed Q-compressions, by simulated annealing (search, not exhaustive) [M
for what it finds: every returned y is verified exactly].

For a b >= 3 GP slice (theory/ab_sweep.build) and a quadruple of Q-compressions z_k (one Q-sum per block and coset of
Q in G'), the state is one letter per G'-orbit and block with exactly those Q-sums; a move changes the letters of two
G'-orbits of one block lying over the same Q-cosets (same coset profile up to a factor) so that every Q-sum is kept.
Energy: sum over d in G' of (sum_k PAF_{G'}(y_k)(d) - T1(d))^2; a state of energy 0 is a stage-1 solution with the
prescribed compressions (verified again with ab_s1sd.S1SD.verify).  Restarts with fresh random states.
usage: python3 theory/ab_s1_anneal.py '<cfg json>' '<z json: 4 x c>' [--seed 1] [--restarts 50] [--iters 2000000]
                                      [--T0 2000] [--T1 2] [--out FILE.jsonl] [--want 20]"""
import sys, os, json, time, math
import numpy as np
import numba as nb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_s1sd as SD


@nb.njit(cache=True)
def _paf(y, add, n):
    out = np.zeros(n, np.int64)
    for d in range(n):
        s = 0
        for g in range(n): s += y[g] * y[add[g, d]]
        out[d] = s
    return out


@nb.njit(cache=True)
def _anneal(Y, lab, pts, ptr, grp_orbs, grp_ptr, coef, lett, lptr, add, T1, iters, Tstart, Tend, seed):
    """Y: 4 x n point values; lab: 4 x nO current letter values per orbit; pts/ptr: points of each orbit;
    groups of exchangeable orbits (grp_orbs/grp_ptr) with coefficients coef (multiplicity of the orbit's coset profile);
    letters per orbit lett/lptr"""
    np.random.seed(seed)
    n = Y.shape[1]; nG = len(grp_ptr) - 1
    P = np.zeros((4, n), np.int64)
    for k in range(4): P[k] = _paf(Y[k], add, n)
    S = P.sum(0)
    E = 0
    for d in range(n): E += (S[d] - T1[d]) ** 2
    best = E; Yb = Y.copy()
    ynew = np.zeros(n, np.int64)
    for it in range(iters):
        if E == 0: break
        Tt = Tstart * (Tend / Tstart) ** (it / iters)
        k = np.random.randint(4)
        gi = np.random.randint(nG)
        a0 = grp_ptr[gi]; a1 = grp_ptr[gi + 1]
        if a1 - a0 < 2: continue
        i1 = a0 + np.random.randint(a1 - a0); i2 = a0 + np.random.randint(a1 - a0)
        if i1 == i2: continue
        j1 = grp_orbs[i1]; j2 = grp_orbs[i2]; c1 = coef[i1]; c2 = coef[i2]
        L0 = lptr[j1]; L1 = lptr[j1 + 1]
        a_new = lett[L0 + np.random.randint(L1 - L0)]
        da = a_new - lab[k, j1]
        if da == 0: continue
        num = c1 * da
        if num % c2 != 0: continue
        b_new = lab[k, j2] - num // c2
        ok = False
        for t in range(lptr[j2], lptr[j2 + 1]):
            if lett[t] == b_new: ok = True
        if not ok: continue
        for g in range(n): ynew[g] = Y[k, g]
        for t in range(ptr[j1], ptr[j1 + 1]): ynew[pts[t]] = a_new
        for t in range(ptr[j2], ptr[j2 + 1]): ynew[pts[t]] = b_new
        Pn = _paf(ynew, add, n)
        En = 0
        for d in range(n): En += (S[d] - P[k, d] + Pn[d] - T1[d]) ** 2
        dE = En - E
        if dE <= 0 or np.random.random() < math.exp(-dE / Tt):
            for d in range(n): S[d] += Pn[d] - P[k, d]
            P[k] = Pn
            for g in range(n): Y[k, g] = ynew[g]
            lab[k, j1] = a_new; lab[k, j2] = b_new
            E = En
            if E < best:
                best = E; Yb = Y.copy()
    return best, Yb


class Annealer:
    def __init__(self, cfg, Z):
        self.cfg = cfg; sl = W.build(cfg); self.sl = sl
        Gp = sl.Gp; n = Gp.n; self.n = n
        qc, Qg = W.qcoset_map(sl, cfg['nQ']); self.qc = qc; self.nq = Qg.n
        self.Z = np.array(Z, np.int64); assert self.Z.shape == (4, self.nq)
        self.add = np.array([[Gp.add(g, d) for d in range(n)] for g in range(n)], np.int64)
        nO = len(sl.breps); self.nO = nO
        self.pts_l = [np.nonzero(sl.borb == j)[0] for j in range(nO)]
        self.pts = np.concatenate(self.pts_l).astype(np.int64)
        self.ptr = np.r_[0, np.cumsum([len(p) for p in self.pts_l])].astype(np.int64)
        self.prof = [np.bincount(qc[p], minlength=self.nq) for p in self.pts_l]   # orbit -> points per coset
        # groups: orbits with proportional coset profiles (same support, same normalised profile)
        groups = {}
        for j in range(nO):
            p = self.prof[j]; g = math.gcd(*[int(x) for x in p if x]) if p.any() else 1
            key = tuple((p // g).tolist()); groups.setdefault(key, []).append((j, g))
        self.groups = list(groups.items())
        self.grp_orbs = np.array([j for _, lst in self.groups for j, _ in lst], np.int64)
        self.coef = np.array([c for _, lst in self.groups for _, c in lst], np.int64)
        self.grp_ptr = np.r_[0, np.cumsum([len(lst) for _, lst in self.groups])].astype(np.int64)
        self.lett_l = [np.array(sorted(L), np.int64) for L in sl.letters]
        self.lett = np.concatenate(self.lett_l); self.lptr = np.r_[0, np.cumsum([len(L) for L in self.lett_l])].astype(np.int64)
        self.T1 = sl.T1.astype(np.int64)

    def random_state(self, rng):
        """letters per orbit with the prescribed Q-sums: per group, a uniform DP sample of letters whose
        coefficient-weighted sum meets the group's (normalised) coset sums"""
        lab = np.zeros((4, self.nO), np.int64)
        for k in range(4):
            for key, lst in self.groups:
                key = np.array(key)
                q0 = int(np.nonzero(key)[0][0])
                target = int(self.Z[k, q0]) // int(key[q0])
                if int(self.Z[k, q0]) % int(key[q0]): return None
                # every coset in the support must agree
                for q in np.nonzero(key)[0]:
                    if int(self.Z[k, q]) != target * int(key[q]): return None
                # DP: sum_j c_j * letter_j = target
                suf = [dict() for _ in range(len(lst) + 1)]; suf[len(lst)] = {0: 1}
                for i in range(len(lst) - 1, -1, -1):
                    j, cj = lst[i]; dct = {}
                    for t, m in suf[i + 1].items():
                        for a in self.lett_l[j]:
                            dct[t + cj * int(a)] = dct.get(t + cj * int(a), 0) + m
                    suf[i] = dct
                if suf[0].get(target, 0) == 0: return None
                rem = target
                for i in range(len(lst)):
                    j, cj = lst[i]
                    opts = [int(a) for a in self.lett_l[j] if suf[i + 1].get(rem - cj * int(a), 0)]
                    w = np.array([suf[i + 1][rem - cj * a] for a in opts], float); w /= w.sum()
                    a = opts[rng.choice(len(opts), p=w)]; lab[k, j] = a; rem -= cj * a
        return lab

    def expand(self, lab):
        Y = np.zeros((4, self.n), np.int64)
        for j in range(self.nO): Y[:, self.pts_l[j]] = lab[:, j][:, None]
        return Y

    def run(self, seed, iters, Tstart, Tend):
        rng = np.random.default_rng(seed)
        lab = self.random_state(rng)
        if lab is None: return None, None
        Y = self.expand(lab)
        best, Yb = _anneal(Y, lab.copy(), self.pts, self.ptr, self.grp_orbs, self.grp_ptr, self.coef, self.lett,
                           self.lptr, self.add, self.T1, iters, float(Tstart), float(Tend), int(seed))
        return int(best), Yb


def main():
    a = sys.argv[1:]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    cfg = SD.cfg_from_json(a[0]); Z = json.loads(a[1])
    A = Annealer(cfg, Z)
    P = SD.S1SD(cfg, sl=A.sl)
    seed0 = int(opt('--seed', '1')); R = int(opt('--restarts', '50')); iters = int(opt('--iters', '2000000'))
    T0 = float(opt('--T0', '2000')); T1 = float(opt('--T1', '2')); want = int(opt('--want', '20'))
    outp = opt('--out', None); fo = open(outp, 'a') if outp else None
    found = 0; t0 = time.time()
    for r in range(R):
        best, Yb = A.run(seed0 + r, iters, T0, T1)
        rec = dict(restart=r, seed=seed0 + r, best=best, secs=round(time.time() - t0, 1))
        if best == 0:
            ok = P.verify(Yb); rec['verified'] = ok
            if ok:
                found += 1; rec['y'] = Yb.tolist()
                if fo: fo.write(json.dumps(dict(cfg=a[0], z=Z, y=Yb.tolist())) + "\n"); fo.flush()
        print(json.dumps({k: v for k, v in rec.items() if k != 'y'}), flush=True)
        if found >= want: break


if __name__ == "__main__":
    main()
