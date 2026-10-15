r"""Stage 1 of the c = 6 family (theory/b599_fam.py) restricted to low-Q-rank blocks: the r-pattern ansatz [M].

A block y on G' = Q x H' is an r-pattern block if its six restrictions y(q, .) (functions on H' / +-) take at most r
distinct values.  For the bordered compressions (347 class z = c + 2 1_S, monomial class z = c + a delta) every block
compression is constant or two-valued, and the exact low-count solutions of the same-shape rung v = 215 are dominated
by such blocks (report/b599_answer.md section 3).  The ansatz is a restriction, not a slice: it is searched exactly
(every r-pattern block with the given compression, then the exact 4-list join theory/p647_join.join4 on the stage-1
class signatures, every hit re-checked by the full G' PAF), so it finds every stage-1 solution whose four blocks are
r-pattern blocks, and nothing else.
usage: python3 theory/b599_twopat.py --n 9 --q 11 [--r 2] [--classes 0,7,8] [--out FILE] [--maxlist 4000000]"""
import sys, os, json, math, itertools, time
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b599_fam import Fam
import b599_flevels as FL
import p647_join as PJ


def coset_patterns(F, Emax):
    """all J-tuples of letters (one coset's pattern) with their sums and energies"""
    T = np.array(list(itertools.product(F.LET.tolist(), repeat=F.J)), np.int64)
    return T, T @ F.WJ, (T * T) @ F.WJ


def rpattern_blocks(F, z, Emax, r=2, maxlist=4 * 10 ** 6):
    """every block vector with Q-compression z whose six coset patterns take at most r distinct values"""
    T, sm, en = coset_patterns(F, Emax)
    z = np.asarray(z, np.int64)
    vals = sorted(set(z.tolist()))
    J = F.J
    blocks = []
    # assign a pattern index (0..r-1) to every coset; cosets with different z need different patterns
    for assign in itertools.product(range(r), repeat=6):
        if assign[0] != 0: continue
        used = sorted(set(assign))
        if used != list(range(len(used))): continue
        # canonical labelling: first occurrence order
        seen = []; ok = True
        for a_ in assign:
            if a_ not in seen: seen.append(a_)
        if seen != sorted(seen): continue
        sums = {}
        for q, a_ in enumerate(assign):
            if a_ in sums and sums[a_] != z[q]: ok = False; break
            sums[a_] = int(z[q])
        if not ok: continue
        cnt = [assign.count(i) for i in range(len(used))]
        # per pattern: candidates with the required sum and energy
        cand = []
        for i in range(len(used)):
            m = (sm == sums[i]) & (en * cnt[i] <= Emax)
            cand.append(np.nonzero(m)[0])
        tot = 1
        for cnd in cand: tot *= len(cnd)
        if tot == 0: continue
        if tot > maxlist: raise RuntimeError(f"pattern product {tot} > maxlist")
        grids = np.meshgrid(*cand, indexing='ij')
        C = np.stack([g_.ravel() for g_ in grids], 1)                       # (ncombo, npat) pattern indices
        e = sum(en[C[:, i]] * cnt[i] for i in range(len(used)))
        keep = e <= Emax
        if len(used) > 1:
            for i in range(len(used)):
                for j_ in range(i + 1, len(used)): keep &= C[:, i] != C[:, j_]
        C = C[keep]
        if len(C) == 0: continue
        V = np.concatenate([T[C[:, assign[q]]] for q in range(6)], 1)
        blocks.append(V)
    if not blocks: return np.zeros((0, 6 * J), np.int64)
    V = np.unique(np.concatenate(blocks), axis=0)
    return V


def solve_z(F, z, r=2, maxlist=4 * 10 ** 6, log=None):
    T0 = int(F.T1SIG[0])
    Ls = [rpattern_blocks(F, z[k], T0, r, maxlist) for k in range(4)]
    sizes = [len(x) for x in Ls]
    if min(sizes) == 0: return dict(sizes=sizes, solutions=0), []
    mins = [int(F.energy(x).min()) for x in Ls]
    for k in range(4):
        cap = T0 - (sum(mins) - mins[k]); Ls[k] = Ls[k][F.energy(Ls[k]) <= cap]
    sizes = [len(x) for x in Ls]
    if min(sizes) == 0: return dict(sizes=sizes, solutions=0), []
    S = [F.sig(x) for x in Ls]
    sols, st = PJ.join4(S, F.T1SIG.astype(np.int64))
    Ys = []
    for s_ in sols:
        Y = np.array([Ls[k][s_[k]] for k in range(4)])
        assert F.check_stage1(Y) and (F.qcomp(Y) == np.asarray(z)).all()
        Ys.append(Y)
    return dict(sizes=sizes, solutions=len(Ys), work_log2=st.get('work_log2'), join_secs=st.get('secs')), Ys


if __name__ == "__main__":
    a = sys.argv[1:]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    n, q = int(opt('--n', '9')), int(opt('--q', '11')); r = int(opt('--r', '2'))
    outp = opt('--out', ''); maxlist = int(opt('--maxlist', str(4 * 10 ** 6)))
    cls = None if '--classes' not in a else [int(x) for x in opt('--classes', '').split(',')]
    F = Fam(n, q)
    done = set()
    if outp and os.path.exists(outp):
        for l in open(outp):
            try: d = json.loads(l); done.add((d['class'], tuple(d['c'])))
            except Exception: pass
    fo = open(outp, 'a') if outp else None
    tot = 0
    for ci, cl, c, z in FL.class_quads(F):
        if cls is not None and ci not in cls: continue
        if (ci, tuple(c)) in done: continue
        t0 = time.time()
        try:
            res, Ys = solve_z(F, z, r, maxlist)
        except RuntimeError as e:
            res, Ys = dict(error=str(e)), []
        tot += len(Ys)
        rec = dict(n=n, q=q, v=F.v, r=r, **{'class': ci}, c=list(c), rows=[int(x) for x in z.sum(1)], **res,
                   secs=round(time.time() - t0, 1))
        if Ys: rec['examples'] = [Y.tolist() for Y in Ys[:20]]
        s = json.dumps(rec); print(json.dumps({k: v for k, v in rec.items() if k != 'examples'}), flush=True)
        if fo: fo.write(s + "\n"); fo.flush()
    print(json.dumps(dict(stage="done", total_solutions=tot)))
