"""Cost model of the flat-pool slot k-tree, from per-rung measurements, for the calibration ladder (rungs, open orders).

Per bucket with k keyed chunk classes (all at their most probable value), with pools of 2^l_A .. 2^l_D blocks:
  pairs_AB = 2^(l_A + l_B - p k),  pairs_CD = 2^(l_C + l_D - p k)      p = bits per side per keyed class
  log2 E_top = l_A + l_B + l_C + l_D - 2 p k - b_r (n_ind - k) - drift(n_ind - k)
where n_ind = (number of lag classes) - 1, b_r = in-bucket rest bits per class and drift(j) the partial-join excess
over the Gaussian.  Both come from the energy shell (each block's PAF vector on a sphere |a|^2 = F):
- variance inflation of the rest inside a bucket, 1 + (k/2) / (n_ind + 1 - k) (sphere model), scaled by 0.55 to
  the in-bucket Gaussian measured at 311 (scripts/bucket_pjoin.py: 2.79 bits per class against the model's 2.86);
- drift: the exact density at 0 of a sum of four uniform sphere points (projected on j coordinates) against the
  Gaussian (scripts/sphere_drift.py) is 0.098 bits per coordinate over all n coordinates [E, matches the measured
  global partial joins at 311 to j = 12]; inside a bucket the measured excess lies between the n = n_ind and the
  n = n_ind - k curves: 0.08 bits per rest class (default rest_drift).
p and the global per-class bits come from the pools' PAF spread: pair-sum std s = sqrt(2) sigma, lattice step 4:
p = log2(sqrt(2 pi) s / 4); b_G = n'_Gauss / n_ind.  Bucket time: t = c_L (2^l_A + 2^l_C) + c_P (pairs_AB + pairs_CD)
seconds (calibrated on the 311 production join), buckets per core = m / E_top (m = tuples a core leaves in a bucket:
reversing one block keeps its PAF, so closure pools hold each core as a cluster of 2^4 x (2 if A, B or C, D have equal
row sums) = 32 tuples; reversal-canonical pools (build_slot_pools*.py --revcanon) leave m = 2 or 1), with the top-bucket rate held for the first ~10^4 buckets (the shares of
scripts/bucket_rank.py fall by about 2x over them).
usage: python3 scripts/ladder_cost.py RUNGS.json      (a list of dicts: name, n_ind, nprime_gauss, sigma, l (4 log2
       pool sizes), m, optional p, rest_drift, k (else best k), t_cal [c_L, c_P])"""
import json, sys, math


def model(r, k):
    lA, lB, lC, lD = r['l']; n = r['n_ind']
    p = r.get('p') or math.log2(math.sqrt(2 * math.pi) * math.sqrt(2) * r['sigma'] / 4) + 0.064   # +0.064: measured at 311 (2.36)
    bG = r['nprime_gauss'] / n
    infl = 1 + (k / 2) / (n + 1 - k)
    br = bG + 0.55 * 0.5 * math.log2(infl)
    j = n - k; drift = r.get('rest_drift', 0.08) * j
    logE = lA + lB + lC + lD - 2 * p * k - br * j - drift
    pAB, pCD = lA + lB - p * k, lC + lD - p * k
    cL, cP = r.get('t_cal', [6.0e-8, 6.0e-8])
    t = cL * (2 ** lA + 2 ** lC) + cP * (2 ** pAB + 2 ** pCD)
    nb = r.get('m', 2) / 2 ** logE
    return dict(k=k, p=round(p, 3), b_rest=round(br, 3), drift=round(drift, 2), log2E_top=round(logE, 2),
                log2_pairs=[round(pAB, 2), round(pCD, 2)], bucket_s=round(t, 2), log2_buckets_per_core=round(math.log2(nb), 2),
                hours_per_core=round(nb * t / 3600, 2), n_prime_true=round(r['nprime_gauss'] + 0.098 * n, 1),
                log2_solution_tuples=round(lA + lB + lC + lD - r['nprime_gauss'] - 0.098 * n, 1), log2_ops_per_core=round(math.log2(nb * (2 ** lA + 2 ** lC + 2 ** pAB + 2 ** pCD)), 1))


def main():
    for r in json.load(open(sys.argv[1])):
        ks = [r['k']] if 'k' in r else range(6, min(13, r['n_ind']))
        best = min((model(r, k) for k in ks), key=lambda z: z['hours_per_core'])
        print(json.dumps(dict(name=r['name'], **best)))


if __name__ == '__main__':
    main()
