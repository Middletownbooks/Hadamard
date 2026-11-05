"""Exact integer planted recovery on 9-switch (A,B,C,tau C) faces.
Also checks the polarization coefficients used in the affine-switch barrier.
No true-target core search; no claimed complexity calibration.
"""
import argparse
import itertools as it
import json
from collections import defaultdict
import numpy as np
from h863_analog import Host
from openai_deficit_symmetry import SW, q_gate, joint_gate, W


def independent_joint_check(m, qs):
    """Cartesian-product row enumeration, independent of the reachable-state DP."""
    alphabet = (-43,-41,-29,-27,-15,-13,-1,1,13,15,27,29,41,43)
    cap = 16*m*42-504
    ordinary = defaultdict(set); reflected = defaultdict(set)
    for row in it.product(alphabet,repeat=m):
        en=sum(a*a-1 for a in row)
        if en<=cap:
            ordinary[sum(row)].add(en)
    weights=(1,)+(2,)*((m-1)//2)
    for row in it.product(alphabet,repeat=len(weights)):
        en=sum(w*(a*a-1) for w,a in zip(weights,row))
        if en<=cap:
            reflected[sum(w*a for w,a in zip(weights,row))].add(en)
    def convolution(sets):
        values={0}
        for s in sets:
            values={a+b for a in values for b in s if a+b<=cap}
        return values
    survivors=[]
    for zs in qs:
        terms=[]
        for z in zs[:2]:
            terms.extend((reflected[z[0]],reflected[z[3]],{2*e for e in ordinary[z[1]]}))
        terms.extend({2*e for e in ordinary[a]} for a in zs[2])
        if cap in convolution(terms):
            survivors.append(zs)
    assert survivors==joint_gate(m,qs)
    return len(survivors)


def exact_paf(host, x):
    ans = []
    g = [(q,u,t) for q in range(4) for u in range(host.m) for t in range(host.p)]
    for a,b,c in host.lagreps:
        perm = [((q^a)*host.m+(u+b)%host.m)*host.p+(t+c)%host.p for q,u,t in g]
        ans.append((x.astype(np.int64)*x[...,perm]).sum((-1,-2)))
    return np.stack(ans,-1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--m', type=int, default=5)
    ap.add_argument('--seed', type=int, default=863)
    args = ap.parse_args()
    m = args.m; host = Host(m,43,4); rng = np.random.default_rng(args.seed)
    tau = [SW[q]*m+(-u)%m for q in range(4) for u in range(m)]
    e = rng.choice((-1,1),size=(4,4*m)).astype(np.int8)
    c = rng.choice((-1,1),size=(4,4*m,6)).astype(np.int8)
    for i in (0,1):
        for k in range(4*m):
            if k < tau[k]:
                e[i,tau[k]]=e[i,k]; c[i,tau[k]]=c[i,k]
    # Three disjoint rectangles per independent generator, on q=01, u=0,1.
    moves = []
    for i in range(3):
        for j in (0,2,4):
            move = []
            for u,s in ((0,1),(1,-1)):
                for col,t in ((j,1),(j+1,-1)):
                    k=m+u; c[i,k,col]=s*t; move.append((i,k,col))
                    if i<2:
                        c[i,tau[k],col]=s*t; move.append((i,tau[k],col))
            moves.append(move)
    e[3]=e[2,tau]; c[3]=c[2,tau]
    variants = np.broadcast_to(c,(512,)+c.shape).copy()
    for mask in range(512):
        for b, move in enumerate(moves):
            if mask>>b&1:
                for i,k,j in move:
                    variants[mask,i,k,j]*=-1
        variants[mask,3]=variants[mask,2,tau]
    assert np.all(variants.sum(-1)==c.sum(-1))
    assert np.all(variants.reshape(512,4,4,m,6).sum(-2)==c.reshape(4,4,m,6).sum(-2))
    x = host.assemble(np.broadcast_to(e,(512,)+e.shape),variants)
    V = exact_paf(host,x)
    assert np.array_equal(V,host.class_paf(x).sum(1))
    plant=int(rng.integers(512)); hits=np.flatnonzero(np.all(V==V[plant],axis=1)).tolist()
    assert plant in hits
    # Integer polarization on the cube: this is 4 times the Boolean quadratic coefficient.
    quadratic = [V[0]-V[1<<a]-V[1<<b]+V[(1<<a)|(1<<b)] for a,b in it.combinations(range(9),2)]
    nonzero=sum(bool(np.any(t)) for t in quadratic)
    assert nonzero>0  # nine common-frequency nonzero directions cannot all be orthogonal in R^8
    # Independently verify every enumerated quotient candidate by Walsh energies.
    qs=q_gate(m)
    for zs in qs:
        energies=[sum(sum(W[j][k]*z[k] for k in range(4))**2 for z in zs) for j in range(4)]
        assert energies==[4*host.n,12,12,12]
    print(json.dumps(dict(m=m,seed=args.seed,face_bits=9,candidates=512,
                         plant=plant,hits=hits,all_exact_paf_agree=True,
                         quadratic_nonzero_pairs=nonzero,pairs=len(quadratic),
                         q_candidates_checked=len(qs),joint_survivors=independent_joint_check(m,qs),
                         independent_cartesian_gate_agrees=True)))


if __name__=='__main__':
    main()
