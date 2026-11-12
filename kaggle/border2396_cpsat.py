#!/usr/bin/env python3
"""
Border search for the 2396 slice (EKLO Theorem-7-type bordered Goethals-Seidel array with a non-classical
character set).  Self-contained: paste into ONE Kaggle notebook cell, or upload and run `!python <file>`.
Needs only numpy and ortools.  On Kaggle: notebook settings -> Internet ON, then `!pip install ortools` once
(numpy is preinstalled).  Written for middletownbooks/hadamard, report/fable_fresh_answer.md section 7.1.

WHAT IT DECIDES
  Q = Z_2^4 (points q = 0..15, characters chi_m(q) = (-1)^popcount(m & q)), T = {0,1,2,3,4,8,12} (7 characters),
  s = 37, n = 4 * (16 * 37 + 7) = 2396, b = 4 * 7 = 28.
  Each CASE is a quadruple of compressions y_0..y_3 : Q -> Z (cell sums of a core on Q x Z_37) that solves the
  Q-level equations   sum_k |yhat_k(chi)|^2 = 28 on T,  = 2396 off T   (checked at start).
  A BORDER is a triple of +-1 matrices K (b x b), E (b x 4*16), Z (b x 4*16) with
     every row of E and of Z: on each of the four blocks a +-1 function on Q with Walsh support inside T,
     E E^T = 64 I,  Z Z^T = 64 I,  K = E M^T Z^T / 64 a +-1 matrix   (M = Goethals-Seidel array of the y's).
  Then the compressed matrix [[K, sqrt(s) E], [sqrt(s) F, M]], F = -Z^T, has Gram n I (checked at the end with
  integers), and ANY core on Q x Z_37 with these compressions completes to a Hadamard matrix of order 2396.
  A border is NOT a Hadamard matrix: the core itself (312 +-1 orbit signs per block) must still be found.

MODES
  'single'  : each block of a border row is one signed character of T (14 choices per block).
  'general' : each block is any +-1 function with Walsh support in T (30 choices per block; contains 'single').
  'single' is the smaller model (it ran 2 h on 1 worker at 7.3 GB RAM without a verdict); 'general' is larger.
  Memory, measured on case 0 'single': 6.6-6.7 GB peak with 1 worker; 2 workers exceeded 15 GB on a 15 GB machine.
  Kaggle CPU sessions have 4 cores and ~30 GB: MODE='single', WORKERS=2.  Run 'general' only with WORKERS=1 at
  first and watch the RAM gauge.
  INFEASIBLE in 'general' would prove that no border of this coset-constant shape exists for that case.

CASES
  The border problem depends only on the Walsh values of the compressions on T; the 2688 Q-level solutions have
  44 distinct profiles, one representative each below (CASES[i], from solution index SOL_INDEX[i]).
  CASE = 'control' runs a known-good small instance first (Q = Z_2^2, T = 3 characters, s = 3, n = 60; a border
  was found here in about 15 s) to check the installation and the code path.

OUTPUT
  /kaggle/working/border2396/progress.jsonl  one line per case: status, seconds (skips cases already logged:
                                             resume by re-running with the previous output attached)
  /kaggle/working/border2396/border_case<i>.json   K, E, Z, the compressions and all checks, if one is found
"""
import os, sys, json, time, itertools
import numpy as np

CONFIG = dict(
    CASES_TO_RUN = ['control', 0],   # 'control' and/or case indices 0..43, e.g. list(range(44))
    MODE = 'single',                 # 'single' or 'general'
    TIME_LIMIT_S = 10 * 3600,        # per case (Kaggle CPU sessions end after 12 h)
    WORKERS = 2,                     # CP-SAT workers: about 6.7 GB RAM EACH in 'single' mode (measured); 2 fits
                                     # Kaggle's ~30 GB, 3 is borderline, 4 will be killed
    SEED = 1,
    OUT_DIR = '/kaggle/working/border2396' if os.path.isdir('/kaggle/working') else './border2396',
    RESUME_FROM = None,              # e.g. '/kaggle/input/<dataset>/border2396' to skip finished cases
    LOG_SEARCH = False,              # True: CP-SAT prints its search log
)

# ----------------------------------------------------------------------------------------------- data
ORBITS = [[0], [1, 2, 3], [4, 8, 12], [5, 10, 15], [6, 11, 13], [7, 9, 14]]   # sigma-orbits of Q = Z_2^4
T599 = [0, 1, 2, 3, 4, 8, 12]
S599 = 37
CASES = [
    [[-17, 5, 5, -5, -1, 1], [-5, 1, 1, 1, -1, -1], [-11, 5, 5, 7, -5, -7], [-1, 1, 1, -5, 5, -1]],
    [[-17, 5, 5, -5, -1, 1], [-5, 1, 1, 1, -1, -1], [-11, 5, 5, 7, -5, -7], [1, -1, -1, 5, -5, 1]],
    [[-17, 5, 5, -5, -1, 1], [-5, 1, 1, 1, -1, -1], [-1, 1, 1, -5, 5, -1], [11, -5, -5, -7, 5, 7]],
    [[-17, 5, 5, -5, -1, 1], [-5, 1, 1, 1, -1, -1], [1, -1, -1, 5, -5, 1], [11, -5, -5, -7, 5, 7]],
    [[-17, 5, 5, -5, -1, 1], [5, -1, -1, -1, 1, 1], [-11, 5, 5, 7, -5, -7], [-1, 1, 1, -5, 5, -1]],
    [[-17, 5, 5, -5, -1, 1], [5, -1, -1, -1, 1, 1], [-11, 5, 5, 7, -5, -7], [1, -1, -1, 5, -5, 1]],
    [[-17, 5, 5, -5, -1, 1], [5, -1, -1, -1, 1, 1], [-1, 1, 1, -5, 5, -1], [11, -5, -5, -7, 5, 7]],
    [[-17, 5, 5, -5, -1, 1], [5, -1, -1, -1, 1, 1], [1, -1, -1, 5, -5, 1], [11, -5, -5, -7, 5, 7]],
    [[5, -1, -1, -1, 1, 1], [17, -5, -5, 5, 1, -1], [-11, 5, 5, 7, -5, -7], [-1, 1, 1, -5, 5, -1]],
    [[5, -1, -1, -1, 1, 1], [17, -5, -5, 5, 1, -1], [-11, 5, 5, 7, -5, -7], [1, -1, -1, 5, -5, 1]],
    [[5, -1, -1, -1, 1, 1], [17, -5, -5, 5, 1, -1], [-1, 1, 1, -5, 5, -1], [11, -5, -5, -7, 5, 7]],
    [[5, -1, -1, -1, 1, 1], [17, -5, -5, 5, 1, -1], [1, -1, -1, 5, -5, 1], [11, -5, -5, -7, 5, 7]],
    [[-17, 5, 5, -5, -1, 1], [-5, 1, 1, 5, -5, -1], [-7, 1, 1, 7, -1, -7], [-1, 1, 1, -1, -7, 7]],
    [[-17, 5, 5, -5, -1, 1], [-5, 1, 1, 5, -5, -1], [-7, 1, 1, 7, -1, -7], [1, -1, -1, 1, 7, -7]],
    [[-17, 5, 5, -5, -1, 1], [-5, 1, 1, 5, -5, -1], [-1, 1, 1, -1, -7, 7], [7, -1, -1, -7, 1, 7]],
    [[-17, 5, 5, -5, -1, 1], [-5, 1, 1, 5, -5, -1], [1, -1, -1, 1, 7, -7], [7, -1, -1, -7, 1, 7]],
    [[-17, 5, 5, -5, -1, 1], [5, -1, -1, -5, 5, 1], [-7, 1, 1, 7, -1, -7], [-1, 1, 1, -1, -7, 7]],
    [[-17, 5, 5, -5, -1, 1], [5, -1, -1, -5, 5, 1], [-7, 1, 1, 7, -1, -7], [1, -1, -1, 1, 7, -7]],
    [[-17, 5, 5, -5, -1, 1], [5, -1, -1, -5, 5, 1], [-1, 1, 1, -1, -7, 7], [7, -1, -1, -7, 1, 7]],
    [[-17, 5, 5, -5, -1, 1], [5, -1, -1, -5, 5, 1], [1, -1, -1, 1, 7, -7], [7, -1, -1, -7, 1, 7]],
    [[5, -1, -1, -5, 5, 1], [17, -5, -5, 5, 1, -1], [-7, 1, 1, 7, -1, -7], [-1, 1, 1, -1, -7, 7]],
    [[5, -1, -1, -5, 5, 1], [17, -5, -5, 5, 1, -1], [-7, 1, 1, 7, -1, -7], [1, -1, -1, 1, 7, -7]],
    [[5, -1, -1, -5, 5, 1], [17, -5, -5, 5, 1, -1], [-1, 1, 1, -1, -7, 7], [7, -1, -1, -7, 1, 7]],
    [[5, -1, -1, -5, 5, 1], [17, -5, -5, 5, 1, -1], [1, -1, -1, 1, 7, -7], [7, -1, -1, -7, 1, 7]],
    [[-13, 5, 5, -5, -5, 5], [1, -1, -1, 7, 5, -11], [-7, 1, 1, 1, -1, -1], [-1, 1, 1, -1, 5, -5]],
    [[-13, 5, 5, -5, -5, 5], [1, -1, -1, 7, 5, -11], [-7, 1, 1, 1, -1, -1], [1, -1, -1, 1, -5, 5]],
    [[-13, 5, 5, -5, -5, 5], [1, -1, -1, 7, 5, -11], [-1, 1, 1, -1, 5, -5], [7, -1, -1, -1, 1, 1]],
    [[-13, 5, 5, -5, -5, 5], [1, -1, -1, 7, 5, -11], [1, -1, -1, 1, -5, 5], [7, -1, -1, -1, 1, 1]],
    [[-1, 1, 1, -11, -1, 11], [1, -1, -1, 1, -5, 5], [-11, 5, 5, 1, -7, 1], [-5, 1, 1, -1, -1, 1]],
    [[-1, 1, 1, -11, -1, 11], [1, -1, -1, 1, -5, 5], [-11, 5, 5, 1, -7, 1], [5, -1, -1, 1, 1, -1]],
    [[-1, 1, 1, -11, -1, 11], [1, -1, -1, 1, -5, 5], [-5, 1, 1, -1, -1, 1], [11, -5, -5, -1, 7, -1]],
    [[-1, 1, 1, -11, -1, 11], [1, -1, -1, 1, -5, 5], [5, -1, -1, 1, 1, -1], [11, -5, -5, -1, 7, -1]],
    [[-1, 1, 1, -7, -1, 7], [-1, 1, 1, -1, 7, -7], [-17, 5, 5, -1, -1, -1], [-1, 1, 1, 7, -7, -1]],
    [[-1, 1, 1, -7, -1, 7], [-1, 1, 1, -1, 7, -7], [-17, 5, 5, -1, -1, -1], [1, -1, -1, -7, 7, 1]],
    [[-1, 1, 1, -7, -1, 7], [-1, 1, 1, -1, 7, -7], [-1, 1, 1, 7, -7, -1], [17, -5, -5, 1, 1, 1]],
    [[-1, 1, 1, -7, -1, 7], [-1, 1, 1, -1, 7, -7], [1, -1, -1, -7, 7, 1], [17, -5, -5, 1, 1, 1]],
    [[-1, 1, 1, -7, -1, 7], [1, -1, -1, 1, -7, 7], [-17, 5, 5, -1, -1, -1], [-1, 1, 1, 7, -7, -1]],
    [[-1, 1, 1, -7, -1, 7], [1, -1, -1, 1, -7, 7], [-17, 5, 5, -1, -1, -1], [1, -1, -1, -7, 7, 1]],
    [[-1, 1, 1, -7, -1, 7], [1, -1, -1, 1, -7, 7], [-1, 1, 1, 7, -7, -1], [17, -5, -5, 1, 1, 1]],
    [[-1, 1, 1, -7, -1, 7], [1, -1, -1, 1, -7, 7], [1, -1, -1, -7, 7, 1], [17, -5, -5, 1, 1, 1]],
    [[1, -1, -1, 1, -7, 7], [1, -1, -1, 7, 1, -7], [-17, 5, 5, -1, -1, -1], [-1, 1, 1, 7, -7, -1]],
    [[1, -1, -1, 1, -7, 7], [1, -1, -1, 7, 1, -7], [-17, 5, 5, -1, -1, -1], [1, -1, -1, -7, 7, 1]],
    [[1, -1, -1, 1, -7, 7], [1, -1, -1, 7, 1, -7], [-1, 1, 1, 7, -7, -1], [17, -5, -5, 1, 1, 1]],
    [[1, -1, -1, 1, -7, 7], [1, -1, -1, 7, 1, -7], [1, -1, -1, -7, 7, 1], [17, -5, -5, 1, 1, 1]]
]
SOL_INDEX = [0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15, 18, 19, 20, 22, 26, 27, 28, 30, 42, 43, 44, 46, 772, 773, 774, 775, 2116, 2117, 2118, 2119, 2160, 2161, 2162, 2163, 2164, 2165, 2166, 2167, 2172, 2173, 2174, 2175]
CONTROL = dict(a=2, s=3, cells=[[1, 1, 1, -1], [1, 1, 1, 1], [1, 1, 1, 3], [1, -1, 1, 1]])   # y_k(q), q = 0..3

# ----------------------------------------------------------------------------------------------- algebra
def chars(k):
    return np.array([[(-1) ** bin(m & q).count('1') for q in range(k)] for m in range(k)], dtype=np.int64)

def gs_array(ys):
    """Goethals-Seidel array over the elementary abelian 2-group Q (inversion is trivial, so R = I and every
    group-developed block X[g, h] = x(g + h) is symmetric)."""
    k = len(ys[0])
    dev = [np.array([[int(y[g ^ h]) for h in range(k)] for g in range(k)], dtype=np.int64) for y in ys]
    A, B, C, D = dev
    return np.block([[A, B, C, D],
                     [-B, A, D.T, -C.T],
                     [-C, -D.T, A, B.T],
                     [-D, C.T, -B.T, A]])

def setup(case):
    if case == 'control':
        a, s, ys = CONTROL['a'], CONTROL['s'], [np.array(y, dtype=np.int64) for y in CONTROL['cells']]
        k = 2 ** a
    else:
        a, s, k = 4, S599, 16
        ys = []
        for cells in CASES[case]:
            y = np.zeros(k, dtype=np.int64)
            for j, ob in enumerate(ORBITS):
                for q in ob: y[q] = cells[j]
            ys.append(y)
    ch = chars(k)
    walsh = [ch @ y for y in ys]
    energy = [int(sum(w[m] ** 2 for w in walsh)) for m in range(k)]
    # T and n are read off the data and then checked against the Q-level equations
    t_guess = sorted(set(energy))[0] // 4
    T = [m for m in range(k) if energy[m] == 4 * t_guess]
    t = len(T); n = 4 * (k * s + t)
    assert t == t_guess and all(energy[m] == (4 * t if m in T else n) for m in range(k)), 'not a Q-level solution'
    if case != 'control': assert T == T599 and n == 2396
    M = gs_array(ys)
    assert np.array_equal(M @ M.T, np.kron(np.eye(4, dtype=np.int64), sum(np.array([[int(sum(y[q] * y[q ^ g] for q in range(k))) for g in range(k)]]) for y in ys)[0][[[g ^ h for h in range(k)] for g in range(k)]]))
    return dict(a=a, s=s, k=k, ys=ys, ch=ch, T=T, t=t, n=n, b=4 * t, M=M)

def block_functions(P, mode):
    k, ch, T = P['k'], P['ch'], P['T']
    if mode == 'single':
        return np.array([e * ch[m] for e in (1, -1) for m in T], dtype=np.int64)
    out = []
    for bits in range(1 << k):
        f = np.array([1 - 2 * ((bits >> q) & 1) for q in range(k)], dtype=np.int64)
        if all(int(ch[m] @ f) == 0 for m in range(k) if m not in T): out.append(f)
    return np.array(out, dtype=np.int64)

# ----------------------------------------------------------------------------------------------- model
def solve_border(P, mode, tl, workers, seed, log=False):
    from ortools.sat.python import cp_model
    k, b, M = P['k'], P['b'], P['M']
    F = block_functions(P, mode); no = len(F)
    MT = M.T
    tabs = {}
    for c in range(4):
        for c2 in range(4):
            tb = F @ MT[c * k:(c + 1) * k, c2 * k:(c2 + 1) * k] @ F.T          # e-block c  .  M^T  .  z-block c2
            tabs[c, c2] = [int(v) for v in tb.reshape(-1)]
    orth = [int(v) for v in (F @ F.T).reshape(-1)]
    lo = min(min(v) for v in tabs.values()); hi = max(max(v) for v in tabs.values())
    m = cp_model.CpModel()
    E = [[m.NewIntVar(0, no - 1, f'e{r}_{c}') for c in range(4)] for r in range(b)]
    Z = [[m.NewIntVar(0, no - 1, f'z{r}_{c}') for c in range(4)] for r in range(b)]
    for r in range(b):                       # K[r, j] = (1/4k) <E_r, M^T Z_j> = +-1
        for j in range(b):
            terms = []
            for c in range(4):
                for c2 in range(4):
                    tb = tabs[c, c2]
                    if not any(tb): continue
                    ix = m.NewIntVar(0, no * no - 1, ''); m.Add(ix == E[r][c] * no + Z[j][c2])
                    v = m.NewIntVar(lo, hi, ''); m.AddElement(ix, tb, v); terms.append(v)
            sg = m.NewBoolVar('')
            m.Add(sum(terms) == 8 * k * sg - 4 * k)
    for X in (E, Z):                         # rows pairwise orthogonal
        for r in range(b):
            for r2 in range(r + 1, b):
                terms = []
                for c in range(4):
                    ix = m.NewIntVar(0, no * no - 1, ''); m.Add(ix == X[r][c] * no + X[r2][c])
                    v = m.NewIntVar(-k, k, ''); m.AddElement(ix, orth, v); terms.append(v)
                m.Add(sum(terms) == 0)
    # symmetry breaking: rows of E in increasing order of their first block (permuting rows of E permutes rows
    # of K and is harmless), and the first block of E row 0 has positive value at q = 0
    for r in range(b - 1): m.Add(E[r][0] <= E[r + 1][0])
    pos0 = [o for o in range(no) if F[o][0] > 0]
    m.AddAllowedAssignments([E[0][0]], [[o] for o in pos0])
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = workers
    s.parameters.random_seed = seed; s.parameters.log_search_progress = log
    t0 = time.time(); st = s.Solve(m); name = s.StatusName(st)
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return name, time.time() - t0, None
    Em = np.array([np.concatenate([F[s.Value(E[r][c])] for c in range(4)]) for r in range(b)])
    Zm = np.array([np.concatenate([F[s.Value(Z[r][c])] for c in range(4)]) for r in range(b)])
    return name, time.time() - t0, (Em, Zm)

# ----------------------------------------------------------------------------------------------- check
def check_border(P, Em, Zm):
    """Integer checks that share nothing with the model: K is +-1 and Hadamard, rows of E, Z have Walsh support in
    T blockwise, and the compressed matrix [[K, sqrt(s)E],[sqrt(s)F, M]] has Gram n I."""
    k, b, s, n, M, ch, T = P['k'], P['b'], P['s'], P['n'], P['M'], P['ch'], P['T']
    num = Em @ M.T @ Zm.T
    if np.any(num % (4 * k)): return dict(ok=False, why='K not integral')
    K = num // (4 * k)
    ok_pm = bool(np.all(np.abs(K) == 1)); ok_had = bool(np.array_equal(K @ K.T, b * np.eye(b, dtype=np.int64)))
    supp = all(int(ch[mm] @ X[r, c * k:(c + 1) * k]) == 0 for X in (Em, Zm) for r in range(b) for c in range(4)
               for mm in range(k) if mm not in T)
    F = -Zm.T
    top = K @ K.T + s * (Em @ Em.T); cross = K @ F.T + Em @ M.T; bot = s * (F @ F.T) + M @ M.T
    gram = bool(np.array_equal(top, n * np.eye(b, dtype=np.int64)) and not cross.any()
                and np.array_equal(bot, n * np.eye(4 * k, dtype=np.int64)))
    return dict(ok=ok_pm and ok_had and supp and gram, K_pm1=ok_pm, K_hadamard=ok_had, walsh_support_in_T=supp,
                compressed_gram_nI=gram, K=K.tolist())

# ----------------------------------------------------------------------------------------------- driver
def main(cfg=CONFIG):
    try:
        import ortools  # noqa: F401
    except ImportError:
        sys.exit('ortools missing: run  !pip install ortools  (notebook Internet must be ON), then re-run')
    os.makedirs(cfg['OUT_DIR'], exist_ok=True)
    prog = os.path.join(cfg['OUT_DIR'], 'progress.jsonl')
    done = set()
    for src in [cfg.get('RESUME_FROM'), cfg['OUT_DIR']]:
        f = src and os.path.join(src, 'progress.jsonl')
        if f and os.path.exists(f):
            for line in open(f):
                r = json.loads(line)
                if r['mode'] == cfg['MODE'] and r['status'] in ('OPTIMAL', 'FEASIBLE', 'INFEASIBLE'):
                    done.add(str(r['case']))
    for case in cfg['CASES_TO_RUN']:
        if str(case) in done:
            print(f'case {case}: already decided, skipped'); continue
        P = setup(case)
        nf = len(block_functions(P, cfg['MODE']))
        print(f"case {case}: |Q| = {P['k']}, T = {P['T']}, s = {P['s']}, n = {P['n']}, b = {P['b']}, "
              f"{nf} block functions ({cfg['MODE']}); Q-level equations verified", flush=True)
        tl = 300 if case == 'control' else cfg['TIME_LIMIT_S']
        st, secs, res = solve_border(P, cfg['MODE'], tl, cfg['WORKERS'], cfg['SEED'], cfg['LOG_SEARCH'])
        rec = dict(case=case, sol_index=(None if case == 'control' else SOL_INDEX[case]), mode=cfg['MODE'],
                   status=st, seconds=round(secs, 1), time_limit=tl, workers=cfg['WORKERS'],
                   utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
        if res is not None:
            chk = check_border(P, *res)
            rec['checks'] = {kk: v for kk, v in chk.items() if kk != 'K'}
            out = dict(rec, ys=[y.tolist() for y in P['ys']], T=P['T'], s=P['s'], n=P['n'],
                       E=res[0].tolist(), Z=res[1].tolist(), K=chk.get('K'))
            json.dump(out, open(os.path.join(cfg['OUT_DIR'], f'border_case{case}.json'), 'w'))
            print(f"case {case}: {st} in {secs:.0f}s; checks {rec['checks']}", flush=True)
            if case != 'control' and chk['ok']:
                print('*** BORDER FOUND for the 2396 slice: report it (file border_case%s.json) ***' % case, flush=True)
        else:
            print(f'case {case}: {st} after {secs:.0f}s' + (' (no border of this kind exists for this case)'
                  if st == 'INFEASIBLE' else ' (undecided: not a zero)'), flush=True)
        with open(prog, 'a') as f: f.write(json.dumps(rec) + '\n')
        if case == 'control' and not (res is not None and rec['checks']['compressed_gram_nI']):
            sys.exit('control failed: stop and report')

if __name__ == '__main__':
    main()
