r"""Gate 1, third (supplementary) check of report/prompt_blockswap_863.md: a direct CP-SAT model of stage 1 of the
block-exchanging ansatz (A, B, C, tau C) at v = 863 for one Q tuple.  Shares no code with the two exhaustive methods.

Credit.  Ansatz: the OpenAI model (report/openai_deficit_answer.md, codex/order-3452-compression, 5b26dae).
863 set-up: OpenAI Codex.  Bordered construction: Kumar's form of Theorem 7 of EKLO.

Model: 42 integer variables (11 tau-orbit values of A, 11 of B, 20 cells of C), each restricted to the alphabet
+-{1,13,15,27,29,41,43}; the four Q row sums of each block fixed; every product y(h) y(h') needed by the 20 H-level lags
as an exact AddMultiplicationEquality; all 20 lags imposed with C contributing PAF_C(g) + PAF_C(tau g).  One worker,
seed 863 + i, at most 300 s.  UNKNOWN decides nothing.
usage: python3 theory/blockswap_stage1_cpsat.py real I1,I2,... OUT.jsonl | planted IN_PLANTED.jsonl OUT.jsonl
"""
import sys, os, json, time
from ortools.sat.python import cp_model

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
VALS = sorted(s * v for v in (1, 13, 15, 27, 29, 41, 43) for s in (1, -1))


def c_of(q, u):
    return 5 * q + u


def tau_of(c):
    q, u = divmod(c, 5)
    return 5 * [0, 2, 1, 3][q] + (-u) % 5


def shift(c, g):
    q, u = divmod(c, 5)
    a, b = divmod(g, 5)
    return 5 * (q ^ a) + (u + b) % 5


def build(zs, tgt, seed, secs):
    md = cp_model.CpModel()
    dom = cp_model.Domain.FromValues(VALS)
    blocks = []
    for bi in range(3):
        cells = [None] * 20
        if bi < 2:
            for c in range(20):
                r = min(c, tau_of(c))
                if cells[r] is None:
                    cells[r] = md.NewIntVarFromDomain(dom, f'b{bi}_{r}')
                cells[c] = cells[r]
        else:
            cells = [md.NewIntVarFromDomain(dom, f'c_{c}') for c in range(20)]
        blocks.append(cells)
    for bi, z in enumerate(zs):
        for q in range(4):
            md.Add(sum(blocks[bi][c_of(q, u)] for u in range(5)) == z[q])
    prod = {}

    def P(bi, x, y):
        key = (bi,) + tuple(sorted((x, y), key=lambda t: t))
        vx, vy = blocks[bi][x], blocks[bi][y]
        nk = (bi, min(vx.Index(), vy.Index()), max(vx.Index(), vy.Index()))
        if nk not in prod:
            w = md.NewIntVar(-43 * 43, 43 * 43, f'p{len(prod)}')
            md.AddMultiplicationEquality(w, [vx, vy])
            prod[nk] = w
        return prod[nk]

    def paf(bi, g):
        return sum(P(bi, c, shift(c, g)) for c in range(20))
    for g in range(20):
        md.Add(paf(0, g) + paf(1, g) + paf(2, g) + paf(2, tau_of(g)) == int(tgt[g]))
    sv = cp_model.CpSolver()
    sv.parameters.num_workers = 1
    sv.parameters.random_seed = seed
    sv.parameters.max_time_in_seconds = secs
    t0 = time.time()
    st = sv.Solve(md)
    out = dict(status=sv.StatusName(st), secs=round(time.time() - t0, 1), n_products=len(prod), seed=seed)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        ys = [[sv.Value(v) for v in b] for b in blocks]
        ys.append([ys[2][tau_of(c)] for c in range(20)])
        ok = all(sum(sum(y[c] * y[shift(c, g)] for c in range(20)) for y in ys) == int(tgt[g]) for g in range(20))
        out.update(y=ys, recheck=ok)
    return out


TGT = [2936 if g == 0 else (-516 if g < 5 else 172) for g in range(20)]

if __name__ == '__main__':
    if sys.argv[1] == 'real':
        d = json.load(open(os.path.join(ROOT, 'logs', 'blockswap', 'gate0_m5.json')))
        Z = [t['z'] for t in d['tuples'] if t['survives']]
        for i in [int(x) for x in sys.argv[2].split(',')]:
            r = build(Z[i][:3], TGT, 863 + i, 300)
            r.update(index=i, z=Z[i])
            print(json.dumps(r), flush=True)
            with open(sys.argv[3], 'a') as f:
                f.write(json.dumps(r) + '\n')
    else:
        for line in open(sys.argv[2]):
            p = json.loads(line)
            r = build(p['z'][:3], p['target'], 863 + p['trial'], 300)
            pl = p['plant']
            r.update(trial=p['trial'], index=p.get('index'), plant_found=(r.get('y', [None])[0:3] == [pl['A'], pl['B'], pl['C']]),
                     note='CP-SAT returns one solution; any verified solution of the planted target counts as recovery '
                          'of the planted instance class (the plant is one of the exhaustive hits)')
            print(json.dumps({k: v for k, v in r.items() if k != 'y'}), flush=True)
            with open(sys.argv[3], 'a') as f:
                f.write(json.dumps(r) + '\n')
