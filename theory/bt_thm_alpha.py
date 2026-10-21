r"""Theorem alpha screen of census records [P] (theory/bt_census2.alpha_trivial_on_K applied to every element of S).

Theorem alpha [P].  k = 4 family, Q = G/K = Z_2^2, |K| > 1.  If S contains sigma fixing K pointwise and acting with
order 3 on Q, there is no core.  Proof (report/b3_push_answer.md, section 3): WLOG sigma has order 3, so G = K + V
with V = {0, g_1, g_2, g_3} ~ Z_2^2 and sigma = id + alpha; every block is x_K on K and the same w on the three other
cosets.  The equation at d in V \ 0 gives sum_i <x_{i,K}, w_i> = 2 - 4|K|, so x_{i,K} = -w_i for three blocks and
x_K, w agree in exactly one point k_0 for the fourth.  The K-characters phi != 1 then give Re(w^(phi) conj phi(k_0)) = 0,
so w is constant, and sum_i (|x_{i,K}^(phi)|^2 + 3 |w_i^(phi)|^2) = 4 != 4n.
usage: python3 theory/bt_thm_alpha.py OUT.jsonl CENSUS.jsonl[,...]   (one line per record: key, dead)
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_census_run as R, bt_census2 as C2


def main():
    out = sys.argv[1]
    done = {}
    if os.path.exists(out):
        for L in open(out):
            d = json.loads(L); done[d['key']] = d
    with open(out, 'a') as fo:
        for f in sys.argv[2].split(','):
            for L in open(f):
                d = json.loads(L)
                if d['key'] in done or not d.get('verdict', '').startswith('passes'):
                    continue
                st = C2.build(d) if 'fam' in d else R.build(d)
                dead = bool(len(st.Kel) > 1 and C2.alpha_trivial_on_K(st))
                rec = dict(key=d['key'], n=d['n'], fam=d.get('fam', 'k4'), r=d['r'], log2E_Q=d.get('log2E_Q'), dead=dead)
                done[d['key']] = rec
                fo.write(json.dumps(rec) + '\n'); fo.flush()


if __name__ == '__main__':
    main()
