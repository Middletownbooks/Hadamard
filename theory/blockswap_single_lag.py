r"""[M] Single-lag exact obstructions for the 5 Gate 1 classes of the (A, B, C, tau C) ansatz at v = 863
(report/blockswap_answer.md section 2.4).  For each class representative and each H-level lag g, the exact set of
(excess, lag-g value) pairs of every block (complete energy-capped lists of theory/blockswap_stage1_exact.py; C
contributes PAF_C(g) + PAF_C(tau g)) is joined: is (2856, T(g)) reachable as A + B + C?  An unreachable lag is an
exact, one-coordinate proof that the tuple is empty (a relaxation of Gate 1).  Also reports the residues mod 8, 16
and 32 of the lag-g sums reachable at total excess 2856.
Credit: ansatz by the OpenAI model (codex/order-3452-compression, 5b26dae); set-up OpenAI Codex; Kumar / EKLO Thm 7.
usage: python3 theory/blockswap_single_lag.py  -> logs/blockswap/single_lag.json"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import blockswap_stage1_exact as X
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
REPS = [0, 8, 16, 24, 32]
out = []
for idx in REPS:
    zA, zB, zC, _ = [tuple(z) for z in X.survivors()[idx]]
    mA, mB, mC = X.min_fixed(zA), X.min_fixed(zB), X.min_free(zC)
    LA = X.fixed_list(zA, X.EXC - mB - mC); LB = X.fixed_list(zB, X.EXC - mA - mC)
    SA, SB = X.sig_block(LA), X.sig_block(LB)
    eA = (LA * LA - 1).sum(1); eB = (LB * LB - 1).sum(1)
    capC = X.EXC - mA - mB
    Rs, Es = [], []
    for s in zC:
        r = X.rows_free(s)
        e = np.array([2 * X.excess(x) for x in r], np.int64)
        o = np.argsort(e, kind='stable')
        Rs.append(np.array(r, np.int64)[o]); Es.append(e[o])
    nC = X._enum_C(*Rs, *Es, capC, X.SH, X.TAULAG, np.ones(X.NC, np.int64), np.zeros(0, np.int64),
                   np.zeros((0, 4), np.int32), np.zeros(0, np.int64))
    kC = np.zeros(nC, np.int64); iC = np.zeros((nC, 4), np.int32); eC = np.zeros(nC, np.int64)
    X._enum_C(*Rs, *Es, capC, X.SH, X.TAULAG, np.ones(X.NC, np.int64), kC, iC, eC)
    # C signatures in chunks (sum over lags of PAF_C(g) + PAF_C(tau g))
    SC = np.zeros((nC, X.NC), np.int64)
    for c0 in range(0, nC, 2_000_000):
        Y = np.concatenate([Rs[q][iC[c0:c0 + 2_000_000, q]] for q in range(4)], 1)
        S = X.sig_block(Y)
        SC[c0:c0 + 2_000_000] = S + S[:, X.TAU]
    rec = dict(index=idx, z=[list(zA), list(zB), list(zC)], lags={})
    for g in range(1, X.NC):
        pa = set(zip(eA.tolist(), SA[:, g].tolist())); pb = set(zip(eB.tolist(), SB[:, g].tolist()))
        pc = set(zip(eC.tolist(), SC[:, g].tolist()))
        ab = {}
        for ea, va in pa:
            for eb, vb in pb:
                ab.setdefault(ea + eb, set()).add(va + vb)
        cby = {}
        for ec, vc in pc:
            cby.setdefault(ec, set()).add(vc)
        reach = set()
        for e_ab, vals in ab.items():
            cs = cby.get(X.EXC - e_ab)
            if cs:
                reach |= {u + w for u in vals for w in cs}
        T = int(X.TGT[g])
        rec['lags'][g] = dict(target=T, reachable=T in reach, n_values=len(reach),
                              res8=sorted({v % 8 for v in reach}), res16=sorted({v % 16 for v in reach}),
                              res32=sorted({v % 32 for v in reach}), target_res=[T % 8, T % 16, T % 32])
    rec['unreachable_lags'] = [g for g, v in rec['lags'].items() if not v['reachable']]
    print(idx, rec['z'], 'unreachable single lags:', rec['unreachable_lags'],
          {g: (v['res16'], v['target_res'][1]) for g, v in rec['lags'].items() if g in (1, 5, 6)}, flush=True)
    out.append(rec)
json.dump(out, open(os.path.join(ROOT, 'logs', 'blockswap', 'single_lag.json'), 'w'))
