"""Low-energy four-list k-tree for skew Goethals-Seidel multiplier slices over any abelian group (Step 3/4).

Pipeline (method of K. Kumar's search/cp_lowe.py + search/listgen, Apache-2.0, re-implemented for general groups;
the join is his search/kt4.c, compiled unchanged):
 1. slice: G = prod Z_n, H = <generators given as integer matrices>; pattern (1; b, c, d) of row sums;
 2. uniform statistics of the PAF class vectors of each block (skew A: uniform pair signs; B, C, D: exact uniform
    sampling at the row sum), conditional means c_k = mu_k + S_k S^+ (0 - sum mu), whitening W_k = S_k^(-1/2) (pinv);
 3. lists: A complete when 2^pairs <= AFULL, else low-energy; B, C, D low-energy members from theory/skew_lg.c
    (Metropolis, energy |W_k (P - c_k)|^2 <= tau_k, tau_k the Gaussian (noncentral chi^2) quantile for L_k / N_k);
 4. the Gaussian estimate inside the lists (means and covariances of the kept keys) [E];
 5. kt4 over NSIG split values of five merge classes (largest spread); every SOL is rebuilt and checked
    (sum PAF = 4v delta by FFT over G, skew A, H-invariance), written as a certificate for scripts/verify_skew.py.
Restart-safe: lists are cached in WORK; kt4 is rerun with the split values not yet logged (--skip).
  python3 theory/skew_kt.py SPEC.json [--L 1e7,1e7,1e7] [--LA 0] [--nsig 3000] [--threads 4] [--seed 1] [--K 1]
SPEC: {"v":.., "mods":[..], "gens":[[[..]]..], "pattern":[1,b,c,d], "tag":".."}
"""
import argparse, json, math, os, subprocess, sys, time, hashlib
import numpy as np
from scipy.stats import ncx2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from skew_lib import Group, Struct, RowSampler

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
SP = os.environ.get('SKEW_WORK', '/tmp/claude-0/-home-user-Hadamard/8cc2c61a-3960-598d-a6cb-a4666d54fa5e/scratchpad')
BIN = os.path.join(SP, 'bin')

def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)

def build_struct(spec):
    mods = spec['mods']; G = Group(mods)
    gens = []
    for M in spec['gens']:
        M = np.array(M, dtype=np.int64).reshape(len(mods), len(mods))
        gens.append(G.index(G.coords @ M.T))
    return Struct(G, gens)

def ntensor(st):
    G = st.G; C, r = st.C, st.r
    N = np.zeros((C, r, r), dtype=np.int32)
    for c, s in enumerate(st.class_rep):
        gs = G.index(G.coords + G.coords[s])
        np.add.at(N[c], (st.orb, st.orb[gs]), 1)
    return N

def keys_of_orbit_signs(st, N, S):
    """S: (n, r) orbit signs -> (n, C) PAF class values, via the quadratic form (exact)."""
    S = S.astype(np.float64)
    return np.rint(np.einsum('ni,cij,nj->nc', S, N, S, optimize=True)).astype(np.int64)

def units_of(st, skew):
    r = st.r
    if skew:
        uo = np.full(r, -1, np.int32); su = np.ones(r, np.int32)
        for u, o in enumerate(st.pairs):
            uo[o] = u; uo[st.onegate[o]] = u; su[o] = 1; su[st.onegate[o]] = -1
        wt = np.array([2 * st.osize[o] for o in st.pairs], np.int32)
    else:
        uo = np.arange(r, dtype=np.int32); su = np.ones(r, np.int32); wt = st.osize.astype(np.int32)
    return uo, su, wt

def whiten(S):
    w, V = np.linalg.eigh(S)
    keep = w > 1e-6 * w.max()
    return (V[:, keep] / np.sqrt(w[keep])) @ V[:, keep].T, int(keep.sum())

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('spec'); ap.add_argument('--L', default='1e7,1e7,1e7'); ap.add_argument('--LA', type=float, default=0)
    ap.add_argument('--AFULL', type=float, default=2 ** 25)
    ap.add_argument('--nsig', type=int, default=3000); ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--K', type=int, default=1)
    ap.add_argument('--samples', type=int, default=20000); ap.add_argument('--T', type=float, default=0.5)
    ap.add_argument('--maxmoves', type=float, default=2e10); ap.add_argument('--gensecs', type=float, default=1800); ap.add_argument('--restart', type=int, default=200000)
    ap.add_argument('--S', default=''); ap.add_argument('--estimate-only', action='store_true')
    ap.add_argument('--maxpairs', type=float, default=3e8); ap.add_argument('--work', default='')
    a = ap.parse_args()
    spec = json.load(open(a.spec)); v = spec['v']; pat = spec['pattern']; tag = spec.get('tag', f'v{v}')
    work = a.work or os.path.join(SP, 'kt', f'{tag}_s{a.seed}'); os.makedirs(work, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    st = build_struct(spec); C, r = st.C, st.r
    assert st.skew_ok and st.gateL
    log(f'v={v} mods={spec["mods"]} |H|={st.Hord} r={r} C={C} pattern={pat} pairs={len(st.pairs)}')
    N = ntensor(st)
    # check the quadratic form against the FFT on random blocks
    Sx = rng.choice(np.array([-1, 1], np.int8), size=(50, r))
    assert np.array_equal(keys_of_orbit_signs(st, N, Sx), st.paf_classes(st.expand(Sx))), 'N tensor mismatch'
    # ---- uniform statistics
    samplers = [None] + [RowSampler(st, b) for b in pat[1:]]
    sizes = [2.0 ** len(st.pairs)] + [s.count for s in samplers[1:]]
    U_S = []
    PA = rng.choice(np.array([-1, 1], np.int8), size=(a.samples, len(st.pairs)))
    U_S.append(st.skew_signs(PA))
    for s in samplers[1:]: U_S.append(s.sample(a.samples, rng))
    K_uni = [keys_of_orbit_signs(st, N, S) for S in U_S]
    mu = [k.mean(0) for k in K_uni]; Sg = [np.cov(k.astype(float), rowvar=False) for k in K_uni]
    Stot = sum(Sg); Sp = np.linalg.pinv(Stot, rcond=1e-9)
    tgt = np.zeros(C)
    cen = [mu[k] + Sg[k] @ Sp @ (tgt - sum(mu)) for k in range(4)]
    Ls = [a.LA or (sizes[0] if sizes[0] <= a.AFULL else 1e7)] + [float(x) for x in a.L.split(',')]
    log('block sizes', [f'{x:.3g}' for x in sizes], 'target lists', [f'{x:.3g}' for x in Ls])
    lists_keys, lists_S = [], []
    for k in range(4):
        cache = os.path.join(work, f'list{k}.npz')
        if os.path.exists(cache):
            z = np.load(cache); lists_keys.append(z['keys']); lists_S.append(z['S']); log(f'block {k}: cached {len(z["keys"])}'); continue
        skew = (k == 0)
        if skew and sizes[0] <= a.AFULL and Ls[0] >= sizes[0]:
            npairs = len(st.pairs)
            idx = np.arange(2 ** npairs, dtype=np.int64)
            P = (1 - 2 * ((idx[:, None] >> np.arange(npairs)) & 1)).astype(np.int8)
            S = st.skew_signs(P)
            keys = np.concatenate([keys_of_orbit_signs(st, N, S[i:i + 50000]) for i in range(0, len(S), 50000)])
            lists_keys.append(keys); lists_S.append(S); np.savez(cache, keys=keys, S=S)
            log(f'block 0: all {len(keys)} skew blocks'); continue
        Wk, rank = whiten(Sg[k])
        lam = float(np.sum((Wk @ (mu[k] - cen[k])) ** 2))
        frac = min(1.0, Ls[k] / sizes[k])
        tau = float(ncx2.ppf(frac, rank, lam)) if frac < 1 else 1e30
        uo, su, wt = units_of(st, skew)
        nstart = 4096
        if skew:
            st_units = rng.choice(np.array([-1, 1], np.int8), size=(nstart, len(st.pairs)))
        else:
            st_units = samplers[k].sample(nstart, rng)
        fin = os.path.join(work, f'gen{k}.bin'); fout = os.path.join(work, f'gen{k}.out')
        with open(fin, 'wb') as f:
            f.write(np.array([r, C, len(wt), int(skew), nstart], np.int32).tobytes())
            f.write(N.astype(np.int32).tobytes()); f.write(Wk.astype(np.float64).tobytes()); f.write(cen[k].astype(np.float64).tobytes())
            f.write(uo.astype(np.int32).tobytes()); f.write(su.astype(np.int32).tobytes()); f.write(wt.astype(np.int32).tobytes())
            f.write(st_units.astype(np.int8).tobytes())
        t0 = time.time()
        cmd = [os.path.join(BIN, 'skew_lg'), fin, fout, str(int(Ls[k])), repr(tau), repr(a.T), str(int(a.maxmoves)),
               str(a.seed * 10 + k), str(a.threads), str(a.restart), str(a.gensecs)]
        log(f'block {k}: generating {int(Ls[k])} of {sizes[k]:.3g} below tau {tau:.3f} (rank {rank}, lambda {lam:.2f})')
        subprocess.run(cmd, check=True)
        with open(fout, 'rb') as f:
            n = int(np.frombuffer(f.read(8), np.int64)[0]); Cc = int(np.frombuffer(f.read(4), np.int32)[0])
            keys = np.frombuffer(f.read(2 * n * Cc), np.int16).reshape(n, Cc).astype(np.int64)
            masks = np.frombuffer(f.read(16 * n), np.uint64).reshape(n, 2)
        U = len(wt)
        bits = np.zeros((n, U), np.int8)
        for u in range(U):
            w, b = (0, u) if u < 64 else (1, u - 64)
            bits[:, u] = ((masks[:, w] >> np.uint64(b)) & np.uint64(1)).astype(np.int8)
        us = (1 - 2 * bits).astype(np.int8)
        S = (st.skew_signs(us) if skew else us)
        chk = keys_of_orbit_signs(st, N, S[:2000])
        assert np.array_equal(chk, keys[:2000]), 'generator keys mismatch'
        if not skew: assert np.all(S.astype(np.int64) @ st.osize == pat[k])
        lists_keys.append(keys); lists_S.append(S); np.savez(cache, keys=keys, S=S)
        log(f'block {k}: kept {n} in {time.time() - t0:.0f}s')
    # ---- Gaussian estimate inside the lists [E]
    m_in = [k.mean(0) for k in lists_keys]; S_in = sum(np.cov(k.astype(float), rowvar=False) for k in lists_keys)
    nvec = st.csize.astype(float)
    Q, _ = np.linalg.qr(np.column_stack([nvec] + [np.eye(C)[:, i] for i in range(C)])); B = Q[:, 1:C]
    S2 = B.T @ S_in @ B; z = B.T @ (tgt - sum(m_in))
    sgn, ld = np.linalg.slogdet(S2)
    logd = -0.5 * float(z @ np.linalg.solve(S2, z)) - 0.5 * (ld + (C - 1) * math.log(2 * math.pi))
    g = int(np.gcd.reduce(st.csize.astype(np.int64)))
    logE_in = (sum(math.log(len(k)) for k in lists_keys) + logd + (C - 1) * math.log(4) + math.log(np.linalg.norm(nvec) / g)) / math.log(10)
    log(f'Gaussian estimate inside the lists {[len(k) for k in lists_keys]}: 10^{logE_in:.2f} [E]')
    json.dump(dict(logE_in=logE_in, lists=[len(k) for k in lists_keys]), open(os.path.join(work, 'estimate.json'), 'w'))
    if a.estimate_only: return
    # ---- kt4
    spread = sum(k.astype(float).var(0) for k in lists_keys)
    Ssel = [int(x) for x in a.S.split(',')] if a.S else sorted(np.argsort(-spread)[:5].tolist())
    assert max(abs(k).max() for k in lists_keys) < 32767
    kin = os.path.join(work, 'kt4_in.bin')
    if not os.path.exists(kin):
        with open(kin, 'wb') as f:
            f.write(np.array([len(k) for k in lists_keys] + [C], np.int32).tobytes())
            f.write(tgt.astype(np.int16).tobytes())
            for k in lists_keys: f.write(k.astype(np.int16).tobytes())
    klog = os.path.join(work, 'kt4.log')
    env = dict(os.environ, OMP_NUM_THREADS=str(a.threads))
    log(f'kt4: merge classes {Ssel}, {a.nsig} split values')
    t0 = time.time()
    with open(klog, 'a') as fh:
        p = subprocess.run([os.path.join(BIN, 'kt4'), kin, ','.join(map(str, Ssel)), str(a.nsig), str(a.K)],
                           stdout=subprocess.PIPE, stderr=fh, text=True, env=env)
    out = p.stdout
    open(os.path.join(work, 'kt4.out'), 'a').write(out)
    log(f'kt4 done in {time.time() - t0:.0f}s: {out.strip().splitlines()[-1] if out.strip() else ""}')
    sols = [list(map(int, l.split()[1:5])) for l in out.splitlines() if l.startswith('SOL')]
    for i, (ia, ib, ic, idd) in enumerate(sols):
        X = [st.expand(lists_S[j][[ix]])[0].astype(np.int64) for j, ix in enumerate((ia, ib, ic, idd))]
        Ptot = sum(st.paf_classes(x[None, :])[0] for x in X)
        sh = st.G.mods
        full = sum(np.rint(np.fft.irfftn(np.abs(np.fft.rfftn(x.reshape(sh).astype(float))) ** 2, s=sh)).astype(np.int64).ravel() for x in X)
        ok = bool(full[0] == 4 * v and np.all(full[1:] == 0))
        log(f'SOL {i}: classes {Ptot.tolist()[:6]}... full PAF ok = {ok}')
        if ok:
            cert = dict(v=v, mods=spec['mods'], gens=spec['gens'],
                        **{n: ''.join('+' if t > 0 else '-' for t in x) for n, x in zip('ABCD', X)})
            cp = os.path.join(work, f'cert_{tag}_{i}.json'); json.dump(cert, open(cp, 'w'))
            hp = os.path.join(ROOT, 'logs', 'skew', 'HITS.jsonl')
            with open(hp, 'a') as f:
                f.write(json.dumps(dict(utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), tag=tag, v=v, cert=cp)) + '\n')
            subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'verify_skew.py'), cp])

if __name__ == '__main__':
    main()
