"""Cost model of the flat-pool slot k-tree, from per-rung measurements, for the calibration ladder (rungs, open orders).

Per bucket with k keyed chunk classes (all at their most probable value), with pools of 2^l_A .. 2^l_D blocks:
  pairs_AB = 2^(l_A + l_B - p k),  pairs_CD = 2^(l_C + l_D - p k)      p = bits per side per keyed class
  log2 E_top = l_A + l_B + l_C + l_D - 2 p k - b_r (n_ind - k) - drift(n_ind - k)
where n_ind = (number of lag classes) - 1, b_r = in-bucket rest bits per class, drift(j) = alpha j + beta j^2 the
partial-join excess over the Gaussian.  The sphere (energy-shell) model gives the in-bucket rest variance inflation
1 + (k/2) / (n_ind + 1 - k) over the global Gaussian per-class variance [E, checked at 311: 2.86 against 2.79 measured].
p and the global per-class bits come from the pools' PAF spread: pair-sum std s = sqrt(2) sigma, lattice step 4:
p = log2(sqrt(2 pi) s / 4); b_G = n'_Gauss / n_ind.  Bucket time: t = c_L (2^l_A + 2^l_C) + c_P (pairs_AB + pairs_CD)
seconds (calibrated on the 311 production join), buckets per core = m / E_top (m = tuples a core leaves in a bucket,
2 when slots A and B have equal row sums), with the top-bucket rate held for the first ~10^4 buckets (the shares of
scripts/bucket_rank.py fall by about 2x over them).
usage: python3 scripts/ladder_cost.py RUNGS.json      (a list of dicts: name, n_ind, nprime_gauss, sigma, l (4 log2
       pool sizes), alpha, beta, m, optional k (else best k), optional t_cal [c_L, c_P])"""
import json, sys, math


def model(r, k):
    lA, lB, lC, lD = r['l']; n = r['n_ind']
    p = r.get('p') or math.log2(math.sqrt(2 * math.pi) * math.sqrt(2) * r['sigma'] / 4)
    bG = r['nprime_gauss'] / n
    infl = 1 + (k / 2) / (n + 1 - k)
    br = bG + 0.5 * math.log2(infl)
    j = n - k; drift = r.get('alpha', 0.004) * j + r.get('beta', 0.0035) * j * j
    logE = lA + lB + lC + lD - 2 * p * k - br * j - drift
    pAB, pCD = lA + lB - p * k, lC + lD - p * k
    cL, cP = r.get('t_cal', [6.0e-8, 6.0e-8])
    t = cL * (2 ** lA + 2 ** lC) + cP * (2 ** pAB + 2 ** pCD)
    nb = r.get('m', 2) / 2 ** logE
    return dict(k=k, p=round(p, 3), b_rest=round(br, 3), drift=round(drift, 2), log2E_top=round(logE, 2),
                log2_pairs=[round(pAB, 2), round(pCD, 2)], bucket_s=round(t, 2), log2_buckets_per_core=round(math.log2(nb), 2),
                hours_per_core=round(nb * t / 3600, 2), log2_ops_per_core=round(math.log2(nb * (2 ** lA + 2 ** lC + 2 ** pAB + 2 ** pCD)), 1))


def main():
    for r in json.load(open(sys.argv[1])):
        ks = [r['k']] if 'k' in r else range(6, min(13, r['n_ind']))
        best = min((model(r, k) for k in ks), key=lambda z: z['hours_per_core'])
        print(json.dumps(dict(name=r['name'], **best)))


if __name__ == '__main__':
    main()
