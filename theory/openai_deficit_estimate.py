"""[E] Synthetic-y diagnostic for (A,B,C,tau C), A,B tau-invariant.
NOT a real-target estimate: no stage-1 certificate is assumed.
python3 theory/openai_deficit_estimate.py --m 5 --N 4000 --seed 863
"""
import argparse
import itertools as it
import json
import math
import numpy as np
from h863_analog import Host, gauss_bits_mean
from openai_deficit_symmetry import SW


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--m', type=int, default=5)
    ap.add_argument('--N', type=int, default=4000)
    ap.add_argument('--seed', type=int, default=863)
    args = ap.parse_args()
    m, N = args.m, args.N
    host = Host(m, 43, 4)
    rng = np.random.default_rng(args.seed)
    tau = [SW[q]*m + (-u)%m for q in range(4) for u in range(m)]
    orb = sorted({tuple(sorted({k, tau[k]})) for k in range(4*m)})
    groups = [orb, orb, [(k,) for k in range(4*m)]]
    independent = [(i, o, len(o)*(2 if i == 2 else 1)) for i, os in enumerate(groups) for o in os]
    excess = host.h_target0() - 16*m
    assert excess % 168 == 0
    # All non-small letters have |y|=13. Select weights summing to excess/168.
    while True:
        order = rng.permutation(len(independent))
        chosen = set(); count = 0
        for idx in order:
            w = independent[idx][2]
            if count+w <= excess//168:
                count += w; chosen.add(int(idx))
            if count == excess//168:
                break
        if count == excess//168:
            break
    E = np.zeros((4,4*m), np.int64); S = np.zeros_like(E)
    c = np.zeros((N,4,4*m,6), np.int8)
    logf = []
    for idx, (i,o,w) in enumerate(independent):
        s = int(rng.choice((-2,2))) if idx in chosen else 0
        e = -s//2 if s else int(rng.choice((-1,1)))
        opts = np.array([[1 if j in a else -1 for j in range(6)]
                         for a in it.combinations(range(6),(s+6)//2)], np.int8)
        draws = opts[rng.integers(len(opts),size=N)]
        for k in o:
            E[i,k], S[i,k], c[:,i,k] = e, s, draws
        logf.append(math.log2(len(opts)))
    E[3] = E[2,tau]; S[3] = S[2,tau]; c[:,3] = c[:,2,tau]
    y = E+7*S
    assert int((y*y).sum()) == host.h_target0()
    x = host.assemble(np.broadcast_to(E,(N,4,4*m)),c)
    V = host.class_paf(x).sum(1)
    smap = host._sum_map(lambda q,u,t:[(q,u,t),(q^3,u,t)],
                         lambda q,u,t:(min(q,q^3),u,t))
    pmap = host._sum_map(lambda q,u,t:[(qq,uu,t) for qq in (q,q^3) for uu in range(m)],
                         lambda q,u,t:(min(q,q^3),t))
    ga, gs, gp = [gauss_bits_mean(v) for v in (V,V@smap.T,V@pmap.T)]
    hd = []
    for i, os in enumerate(groups):
        seen = set(); h = 0
        for o in os:
            oz = tuple(sorted(((k//m)^3)*m+k%m for k in o))
            assert oz in os
            if o == oz or o in seen:
                continue
            seen.update((o,oz))
            h += host.hD_pair(int(S[i,o[0]]),int(S[i,oz[0]]))
        hd.append(h)
    print(json.dumps(dict(m=m,N=N,seed=args.seed,read_at='synthetic mean; not feasible y',
                         log2_fibre=sum(logf),hD_by_generator=hd,hD=sum(hd),
                         C_all=ga,C_S=gs,C_parent=gp,
                         F=gs['bits']-gp['bits'],D_def=ga['bits']-gs['bits']-sum(hd),
                         log2E_synthetic=sum(logf)-ga['bits'])))


if __name__ == '__main__':
    main()
