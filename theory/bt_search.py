r"""Restart-safe production chain search for a bordered-GS structure [M for every core found; rates measured].

Plan: Q-level -> level A (an invariant N_A <= K; lists refined inside each Q-level solution, product over the
Q-cells) -> level B (invariant N_B; distinct level-B vectors inside each level-A fibre by the DP over the level-A
cells) -> expansion of every block's (A, B) fibre (meet in the middle over the A-cells), per-block character filter
-> final exact join on all classes, split first by the level-C compressions (invariant N_C).  Every core is checked
on G (sum of PAFs = target), bordered (theory/bt_border.solve_border) and assembled; H H^T = 4n I is checked; cores
and matrices are written immediately.

Order of work is randomised (seeded) so partial runs sample the structure.  Progress is a JSON-lines log; a restart
skips the (pattern, Q-solution, A, B) items already logged as done.
usage: python3 theory/bt_search.py SPEC.json OUTDIR [--worker i --workers W] [--seed s] [--hours h]
SPEC: {"mods": [...], "Kgens": [[...]], "gens": [matrices] or per-block, "levels": {"A": |N_A|, "B": |N_B|, "C": |N_C|},
       "patterns": [[r1,r2,r3,r4], ...] (signed or sorted; optional), "name": "..."}
"""
import sys, os, json, time, argparse, math, hashlib
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_engine as E, bt_chain as C, bt_border as BB, bt_gates as BG


def build_spec(spec):
    return E.BStruct(tuple(spec['mods']), [tuple(g) for g in spec['Kgens']], spec.get('gamma', 'k'),
                     [np.array(M) if np.array(M).ndim == 2 else [np.array(x) for x in M] for M in spec['gens']]
                     if not spec.get('per_block') else [[np.array(x) for x in bl] for bl in spec['gens']],
                     name=spec.get('name', ''))


def pick_subgroup(st, size, inside=None):
    for N in st.invariant_subgroups():
        if len(N) == size and (inside is None or np.isin(N, inside).all()):
            return N
    raise ValueError(f"no invariant subgroup of order {size}")


class Run:
    def __init__(self, spec, outdir, worker=0, workers=1, seed=1, border=True):
        self.spec = spec
        self.border = border
        self.st = st = build_spec(spec)
        self.out = outdir
        os.makedirs(outdir, exist_ok=True)
        self.worker, self.workers = worker, workers
        self.rng = np.random.default_rng(seed + 1000 * worker)
        lv = spec['levels']
        self.qQ = st.quotient(st.Kel)
        self.qA = st.quotient(pick_subgroup(st, lv['A'], st.Kel))
        self.qB = st.quotient(pick_subgroup(st, lv['B'])) if lv.get('B') else None
        self.qC = st.quotient(pick_subgroup(st, lv['C'])) if lv.get('C') else None
        self.chQ = C.char_matrix(self.qQ); self.chA = C.char_matrix(self.qA)
        self.chB = C.char_matrix(self.qB) if self.qB is not None else None
        self.cvA = [C.coord_values(self.qA, i) for i in range(4)]
        self.PHI = [C.char_orbit_matrix(st, i) for i in range(4)]
        self.logf = os.path.join(outdir, f'log_w{worker}.jsonl')
        self.coref = os.path.join(outdir, f'cores_w{worker}.jsonl')
        self.done = set()
        for f in os.listdir(outdir):
            if f.startswith('log_w') and f.endswith('.jsonl'):
                for L in open(os.path.join(outdir, f)):
                    try:
                        d = json.loads(L)
                        if d.get('event') == 'combo_done':
                            self.done.add(d['item'])
                    except Exception:
                        pass
        self.ncores = 0

    def log(self, **k):
        k['t'] = round(time.time(), 1)
        with open(self.logf, 'a') as f:
            f.write(json.dumps(k) + '\n')

    def q_solutions(self):
        st = self.st
        pats = self.spec.get('patterns') or [list(p) for p in BG.row_patterns(st)]
        out = []
        for rows in pats:
            L = [C.q_level_list(self.qQ, i, rows[i], self.chQ) for i in range(4)]
            if min(len(l[0]) for l in L) == 0:
                continue
            P = [C.paf_reduced(self.qQ, i, L[i][0]) for i in range(4)]
            sol, no = C.join4(P, self.qQ.cls_target)
            for k, s in enumerate(sol):
                out.append((tuple(rows), k, [L[i][0][s[i]] for i in range(4)]))
        return out

    def level_A(self, zs):
        st = self.st
        LA = [C.refine_list(self.qQ, zs[i], self.qA, i, chars=self.chA, cvals=self.cvA[i]) for i in range(4)]
        if min(len(l[0]) for l in LA) == 0:
            return LA, np.zeros((0, 4), np.int64)
        PA = [C.paf_reduced(self.qA, i, LA[i][0]) for i in range(4)]
        solA, noA = C.join4(PA, self.qA.cls_target)
        return LA, solA

    def level_B(self, yA):
        st = self.st
        LB = []
        for i in range(4):
            Y, cnt = C.fibre_dp(st, i, self.qA, yA[i], [], [], self.qB)
            en = np.abs(Y[:, self.qB.borb[i]].astype(float) @ self.chB[0].T) ** 2
            ok = np.all(en <= self.chB[1][None, :] + 1e-6, axis=1)
            LB.append((Y[ok], cnt[ok]))
        if min(len(l[0]) for l in LB) == 0:
            return LB, np.zeros((0, 4), np.int64)
        PB = [C.paf_reduced(self.qB, i, LB[i][0]) for i in range(4)]
        solB, noB = C.join4(PB, self.qB.cls_target)
        return LB, solB

    def finish(self, yA, yB, item):
        st = self.st
        t0 = time.time()
        if self.qB is not None:
            EPS = [C.expand_fibre(st, i, self.qA, yA[i], [self.qB], [yB[i]]) for i in range(4)]
        else:
            EPS = [C.expand_fibre(st, i, self.qA, yA[i], [], []) for i in range(4)]
        n0 = [len(e) for e in EPS]
        EPS = [EPS[i][C.char_filter(st, i, EPS[i], self.PHI[i])] for i in range(4)]
        t1 = time.time()
        fs = C.final_join(st, EPS, self.qC)
        t2 = time.time()
        for f in fs:
            xs = [st.expand(i, EPS[i][f[i]]) for i in range(4)]
            assert st.check_core(xs)
            rec = dict(item=item, eps=[EPS[i][f[i]].astype(int).tolist() for i in range(4)],
                       rows=[int(x.sum()) for x in xs])
            r = BB.solve_border(st, xs) if self.border else None
            if r is not None:
                X, Rh, Ch, k = r
                H = BB.assemble(st, xs, X, Rh, Ch)
                rec['hadamard'] = BB.is_hadamard(H)
                rec['border'] = dict(X=X.tolist(), Rh=Rh.tolist(), Ch=Ch.tolist())
            else:
                rec['hadamard'] = None
            with open(self.coref, 'a') as fo:
                fo.write(json.dumps(rec) + '\n')
            self.ncores += 1
        return dict(fibres=n0, filtered=[len(e) for e in EPS], cores=len(fs),
                    t_expand=round(t1 - t0, 1), t_join=round(t2 - t1, 1))

    def run(self, hours=1e9, max_items=None):
        t_end = time.time() + hours * 3600
        qs = self.q_solutions()
        self.log(event='start', q_solutions=len(qs), spec=self.spec.get('name'))
        order = self.rng.permutation(len(qs))
        k_items = 0
        for qi in order:
            rows, k, zs = qs[qi]
            if hash((rows, k)) % self.workers != self.worker and self.workers > 1:
                pass
            LA, solA = self.level_A(zs)
            self.log(event='levelA', rows=rows, qsol=k, lists=[len(l[0]) for l in LA], solutions=int(len(solA)))
            for ai in self.rng.permutation(len(solA)):
                a = solA[ai]
                yA = [LA[i][0][a[i]] for i in range(4)]
                if self.qB is not None:
                    LB, solB = self.level_B(yA)
                    self.log(event='levelB', rows=rows, qsol=k, A=int(ai), lists=[len(l[0]) for l in LB],
                             combos=int(len(solB)))
                    combos = [(int(bi), [LB[i][0][solB[bi][i]] for i in range(4)]) for bi in self.rng.permutation(len(solB))]
                else:
                    combos = [(0, None)]
                for bi, yB in combos:
                    item = f"{list(rows)}|{k}|{int(ai)}|{bi}"
                    h = int(hashlib.sha1(item.encode()).hexdigest()[:8], 16)
                    if h % self.workers != self.worker or item in self.done:
                        continue
                    res = self.finish(yA, yB, item)
                    self.log(event='combo_done', item=item, **res)
                    k_items += 1
                    if time.time() > t_end or (max_items and k_items >= max_items):
                        self.log(event='stop', items=k_items, cores=self.ncores)
                        return
        self.log(event='exhausted', items=k_items, cores=self.ncores)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('spec'); ap.add_argument('outdir')
    ap.add_argument('--worker', type=int, default=0); ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--hours', type=float, default=1e9)
    ap.add_argument('--max-items', type=int, default=None)
    a = ap.parse_args()
    spec = json.load(open(a.spec))
    Run(spec, a.outdir, a.worker, a.workers, a.seed).run(a.hours, a.max_items)


if __name__ == '__main__':
    main()
