r"""Rigidity screen of census records [P].

Lemma R [P].  Let N <= K be an S-invariant subgroup with |N| >= 5, and suppose S (for a split structure: every S_i)
contains an element tau acting transitively on N \ 0 whose fixed points form a complement M of N (|Fix(tau)| =
|G|/|N|).  Then there is no core.  Proof: with G = M + N, every block is constant on m + (N \ 0) for each m in M
(value c_m) with x(m) = d_m, so for nu_0 in N \ 0, PAF(nu_0) = sum_m [(|N| - 2) c_m^2 + 2 c_m d_m] >= (|N| - 4)|M| > 0,
while the target sum_i PAF_i(nu_0) = -Gamma(0) = -4t < 0.  (Transitivity on N \ 0 forces N elementary abelian:
Z_p with a primitive multiplier, or Z_p^2 with a Singer cycle, as in the 971 hosts.)
usage: python3 theory/bt_thm_rigid.py OUT.jsonl CENSUS.jsonl[,...]   (one line per surviving record: key, dead, N)
"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_engine as E, bt_census_run as R, bt_census2 as C2


def rigid(st, cap=20000):
    G = st.G
    blocks = range(4) if st.split else [0]
    els = {i: E.group_elements(st.bperms[i], G.n, cap=cap) for i in blocks}
    for N in st.invariant_subgroups():
        if len(N) < 5 or not np.isin(N, st.Kel).all():
            continue
        Nset = set(int(x) for x in N)
        nu = int(N[1]) if N[0] == 0 else int(N[0])
        ok_all = True
        for i in blocks:
            found = False
            for p in els[i]:
                if int(np.count_nonzero(p == np.arange(G.n))) * len(N) != G.n:
                    continue
                orb = {nu}; x = nu
                while True:
                    x = int(p[x])
                    if x == nu:
                        break
                    orb.add(x)
                if len(orb) == len(N) - 1 and orb <= Nset:
                    found = True
                    break
            if not found:
                ok_all = False
                break
        if ok_all:
            return int(len(N))
    return None


def main():
    out = sys.argv[1]
    done = set()
    if os.path.exists(out):
        for L in open(out):
            done.add(json.loads(L)['key'])
    with open(out, 'a') as fo:
        for f in sys.argv[2].split(','):
            for L in open(f):
                d = json.loads(L)
                if d['key'] in done or not d.get('verdict', '').startswith('passes'):
                    continue
                st = C2.build(d) if 'fam' in d else R.build(d)
                Nr = rigid(st)
                rec = dict(key=d['key'], n=d['n'], fam=d.get('fam', 'k4'), r=d['r'], log2E_Q=d.get('log2E_Q'),
                           dead=Nr is not None, N=Nr)
                done.add(d['key'])
                fo.write(json.dumps(rec) + '\n'); fo.flush()


if __name__ == '__main__':
    main()
