r"""Census of K-invariant T-matrices of order t (K of order k) per support split, against the Gaussian [M vs E].

For random supports S (scripts/tmat_kt.py: T1, T2 on S as (P +- R)/2, T3, T4 on S^c), the complete k-tree join
(scripts/kt_generic.kt_join4, up to --cap solutions) counts every T-matrix quadruple with supp T1 u supp T2 = S and
the side lists' sign conventions; each one is re-verified exactly (disjoint supports, sum PAF = t delta).
Prediction for a random S of size |S| (labelled count E = 2^log2E from scripts/tmat_estimate.py):
    E * P(|supp T1 u T2| = |S|) / C(m, |S|), with P from the binomial(m, 1/2) model of the union size,
divided by 4 for the two first-orbit sign conventions of P and R (the join finds (P, R) and (R, P) separately: x 2).
Restart-safe: one JSONL line per S (logs/b917/tmat_census_<t>_<k>.jsonl).
usage: python3 theory/b917_tmat_census.py --t 131 --k 5 --nS 30 --seed 11 [--cap 100000] [--buckets 200]"""
import sys, os, json, time, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
os.chdir(ROOT)
from tmat_kt import orbit_setup, side_list, tm_check
from kt_generic import kt_join4


def count_S(t, orbits, reps, Sset, nb, seed, cap):
    """complete exact 2+2 join: all (P, R) on S against all (P', R') on S^c, PAF(P)+PAF(R)+PAF(P')+PAF(R') = 0 at
    every nonzero class (64-bit linear hash, every match re-verified as T-matrices)"""
    m = len(orbits); Sc = [o for o in range(m) if o not in Sset]
    sgS, PS = side_list(t, orbits, Sset, reps)
    sgC, PC = side_list(t, orbits, Sc, reps)
    rng = np.random.default_rng(seed)
    wv = rng.integers(1, 2 ** 62, size=PS.shape[1], dtype=np.int64)
    hS = PS.astype(np.int64) @ wv; hC = PC.astype(np.int64) @ wv
    HL = (hS[:, None] + hS[None, :]).ravel()                      # left pairs (P, R), index i * nS + j
    o = np.argsort(HL, kind='stable'); HLs = HL[o]; del HL
    nS, nC = len(hS), len(hC)
    found = []
    step = max(1, (1 << 22) // nC)
    for i0 in range(0, nC, step):
        HR = (-(hC[i0:i0 + step][:, None] + hC[None, :])).ravel()
        lo = np.searchsorted(HLs, HR, 'left'); hi = np.searchsorted(HLs, HR, 'right')
        for r in np.nonzero(hi > lo)[0]:
            ic, jc = divmod(int(r), nC); ic += i0
            for p in o[lo[r]:hi[r]]:
                iS, jS = divmod(int(p), nS)
                P1, R1, P2, R2 = sgS[iS], sgS[jS], sgC[ic], sgC[jc]
                T = np.zeros((4, t), np.int64)
                for j, ob in enumerate(Sset):
                    T[0, orbits[ob]] = (int(P1[j]) + int(R1[j])) // 2; T[1, orbits[ob]] = (int(P1[j]) - int(R1[j])) // 2
                for j, ob in enumerate(Sc):
                    T[2, orbits[ob]] = (int(P2[j]) + int(R2[j])) // 2; T[3, orbits[ob]] = (int(P2[j]) - int(R2[j])) // 2
                if tm_check(T): found.append(T)
                if len(found) >= cap: return found
    return found


def main():
    a = sys.argv[1:]
    opt = lambda n, d: a[a.index(n) + 1] if n in a else d
    t = int(opt('--t', '131')); k = int(opt('--k', '5')); nS = int(opt('--nS', '30')); seed = int(opt('--seed', '11'))
    cap = int(opt('--cap', '100000')); nb = int(opt('--buckets', '200'))
    K, orbits, reps = orbit_setup(t, k); m = len(orbits)
    out = os.path.join(ROOT, f'logs/b917/tmat_census_{t}_{k}.jsonl')
    done = set()
    if os.path.exists(out):
        for l in open(out):
            try: d = json.loads(l)
            except Exception: continue
            if 'S' in d: done.add(tuple(d['S']))
    rng = np.random.default_rng(seed)
    for i in range(nS):
        size = m // 2 + int(rng.integers(0, 2)) * (m % 2)
        Sset = sorted(rng.choice(m, size, replace=False).tolist())
        if tuple(Sset) in done: continue
        t0 = time.time()
        sols = count_S(t, orbits, reps, Sset, nb, seed * 1000 + i, cap)
        # distinct T-matrix quadruples (as labelled arrays)
        uniq = {T.astype(np.int8).tobytes() for T in sols}
        rec = dict(t=t, k=k, seed=seed, i=i, S=Sset, size=len(Sset), solutions=len(sols), distinct=len(uniq),
                   secs=round(time.time() - t0, 1))
        with open(out, 'a') as fh: fh.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    import resource
    resource.setrlimit(resource.RLIMIT_AS, (int(3.5 * 2 ** 30), int(3.5 * 2 ** 30)))
    main()
