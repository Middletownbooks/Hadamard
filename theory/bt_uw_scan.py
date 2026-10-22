r"""Feasibility of Kumar's Z_2-split (U/W, search/NOTES.md section 19, the route that made 2148) at the census hosts.

For an S-fixed involution z put U = (x + x(. + z))/2 on G/<z> and W = (x - x(. + z))/2 (on a transversal).  An S-orbit
O of G with O + z = O carries U = +-1 and W = 0; a pair of orbits O != O + z is one U-cell, with U in {0, +-1} and one
free sign of W where U = 0.  The class equations split accordingly: the class map c -> c + z (it commutes with
<S, -1>, so it permutes the classes) pairs the equations, and
    U-conditions  v_c + v_{c+z} = t_c + t_{c+z}   (the compression to G/<z>; the characters trivial on z),
    W-conditions  v_c - v_{c+z} = t_c - t_{c+z}   (c != c + z; the characters nontrivial on z).
The lag-0 equation of the compression fixes the number of zero points of U over the four blocks [P]:
    sum_k |U_k^{-1}(0)| = 2|G| - (F(0) + F(z))/4 = |G| + Gamma(z + K)/4      (= n when z is in K).
Printed per host and z:
- cells: U-cells per block (paired / z-stable), e_U and e_W (independent U- and W-equations);
- zc: expected zero cells per block = zero points / (4 x the mean size of the paired cell containing a random point);
  a U-solution leaves W-lists of
  about 2^zc sign vectors per block before the coset-sum constraints, and an exact W-join of about 2^(2 zc) pairs;
- log2 P(W | U) [E]: the probability that an x-quadruple meeting the U-conditions meets the W-conditions too, from the
  lattice Gaussian of the summed class-PAF vectors (blocks uniform, theory/q523_gauss.gauss_bits) applied to all
  equations and to the U-equations only; log2 completions per U-solution = 4 zc + log2 P(W | U).
Calibration: Kumar's 2148 run observed 3.5e-6 completions per U-solution (one in 284,292), with lists of about 200
vectors per block after the constraints; run with --cert on its certificate.
Credit: the split and its use at 2148 are Kumar's (Apache-2.0, github.com/vatskishan03/hadamard); this code is new.
usage: python3 theory/bt_uw_scan.py [--keys k1,k2|@FILE] [--top 2] [--min-logE 40] [--cert CERT.json] [--N 3000]
"""
import sys, os, json, math, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_engine as E, bt_chain as C, bt_census_run as R, bt_census2 as C2
import q523_gauss as QG


def fixed_involutions(st):
    G = st.G
    inv = [g for g in range(1, G.n) if G.add[g, g] == 0]
    return [g for g in inv if all(int(p[g]) == g for bl in st.bperms for p in bl)]


def scan(st, z, N=3000, seed=1):
    G = st.G
    cmap = st.cls[G.add[st.cls_rep, z]]                      # class of rep_c + z
    rows_u, rows_w = [], []
    seen = set()
    for c in range(1, st.ncls):
        if c in seen:
            continue
        c2 = int(cmap[c])
        seen.update((c, c2))
        if c2 == 0:                                          # the pair {0, z}: v_z = F(z), the support count of U
            rows_u.append({c: 1})
        elif c2 == c:
            rows_u.append({c: 1})
        else:
            rows_u.append({c: 1, c2: 1}); rows_w.append({c: 1, c2: -1})
    Fu = np.zeros((len(rows_u), st.ncls), np.int64)
    for i, r in enumerate(rows_u):
        for c, v in r.items():
            Fu[i, c] = v
    # U-cells: orbits O of block 0 with O + z = O (stable) or pairs {O, O + z}
    orb = st.orb[0]; sizes = st.sizes[0]
    partner = orb[G.add[st.reps[0], z]]
    stable = int((partner == np.arange(len(sizes))).sum())
    paired = (len(sizes) - stable) // 2
    psz = [int(sizes[o]) for o in range(len(sizes)) if partner[o] != o]
    mean_sz = float(np.sum(np.square(psz)) / np.sum(psz)) if psz else float('nan')   # size of the cell of a random point
    zero_pts = G.n + int(st.gamma[st.coset[z]]) // 4
    zc = zero_pts / (4 * mean_sz) if psz else 0.0
    rng = np.random.default_rng(seed)
    V = np.zeros((N, st.ncls), np.int64)
    for i in range(4):
        X = rng.choice(np.array([-1, 1], np.int8), size=(N, st.r[i]))
        V += C.class_paf_orbits(st, i, X)
    ga = QG.gauss_bits(V[:, 1:], st.cls_target[1:])
    gu = QG.gauss_bits(V @ Fu.T, Fu @ st.cls_target)
    pwu = gu['bits'] - ga['bits'] if np.isfinite(ga['bits']) and np.isfinite(gu['bits']) else None
    return dict(z=int(z), z_in_K=bool(np.isin(z, st.Kel)), cells=paired + stable, paired=paired, stable=stable,
                e=st.e, e_U=len(rows_u), e_W=len(rows_w), zero_points=zero_pts, mean_cell=round(mean_sz, 2),
                zc=round(zc, 1), log2E_all=round(sum(st.r) - ga['bits'], 2) if np.isfinite(ga['bits']) else None,
                log2P_W_given_U=round(pwu, 1) if pwu is not None else None,
                log2_completions_per_U=round(4 * zc + pwu, 1) if pwu is not None else None,
                rank_all=ga.get('rank'), rank_U=gu.get('rank'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--keys', default=None); ap.add_argument('--top', type=int, default=2)
    ap.add_argument('--min-logE', type=float, default=40); ap.add_argument('--cert', default=None)
    ap.add_argument('--N', type=int, default=3000)
    ap.add_argument('--out', default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                                  'logs/bt/uw_scan.jsonl'))
    a = ap.parse_args()
    hosts = []
    if a.cert:
        cert = json.load(open(a.cert))
        st = E.BStruct(tuple(cert['group']), [tuple(g) for g in cert['subgroup_generators']], 'all',
                       [np.array(M) for M in cert['symmetries']], name='2148')
        hosts.append((dict(key='control_2148', n=537, fam='z3all', r=st.r[0], log2E_Q=None), st))
    else:
        import bt_summary as SU
        D = SU.Data()
        if a.keys:
            ks = open(a.keys[1:]).read() if a.keys.startswith('@') else a.keys
            keys = [k for k in ks.replace('\n', ',').split(',') if k]
            recs = [D.census[k] for k in keys if k in D.census]
        else:
            recs = []
            for (n, fam), L in sorted(D.groups().items()):
                L = [d for d in L if D.status(d) == 'alive' and SU.e_of(d) >= a.min_logE and fam != 'k4split']
                recs += sorted(L, key=lambda d: -SU.e_of(d))[:a.top]
        for d in recs:
            hosts.append((d, C2.build(d) if 'Kgens' in d else R.build(d)))
    with open(a.out, 'a') as fo:
        for d, st in hosts:
            zs = fixed_involutions(st)
            if not zs:
                res = dict(key=d['key'], n=d['n'], fam=d.get('fam', 'k4'), r=d['r'], log2E_Q=d.get('log2E_Q'),
                           note='no S-fixed involution')
                fo.write(json.dumps(res) + '\n'); print(json.dumps(res), flush=True)
                continue
            for z in zs:
                res = dict(key=d['key'], n=d['n'], fam=d.get('fam', 'k4'), r=d['r'], log2E_Q=d.get('log2E_Q'),
                           mods=list(st.G.mods), **scan(st, z, N=a.N))
                fo.write(json.dumps(res) + '\n'); print(json.dumps(res), flush=True)


if __name__ == '__main__':
    main()
