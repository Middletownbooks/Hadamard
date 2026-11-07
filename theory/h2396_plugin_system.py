r"""Gate A1 of report/prompt_599_plugin.md: the Paley plug-in system for Legendre quadruples (b = 1) on G = K x F_q.

Credit: the plug-in decomposition, the K-level equations (I), (J), (Q) and the defect identity are K. Kumar's
(vatskishan03/hadamard, search/NOTES.md sections 81, 85, 86, 92; Apache-2.0).  This file re-derives them for the
b = 1 target and checks them; it shares no code with Kumar's programs or with theory/ab_lib.

Set-up [P].  G = K x F_q, H = {1} x S (S = squares of F_q^*).  A block constant on H-orbits is three functions on K:
x0 on K x {0}, x1 on K x S, x2 on K x N.  With e = (x1 + x2)/2, t = (x1 - x2)/2 (exactly one is non-zero at each k)
and s = x0 - e, the block is the group-ring element  X = s (x) delta_0 + t (x) chi + e (x) 1_F  (chi the Legendre
symbol, chi(0) = 0).  In Z[F_q]:  delta*f = f,  1*1 = q 1,  1*chi = 0,  chi^(-1) = eps chi  (eps = chi(-1)),
chi * chi^(-1) = q delta - 1.  Hence (tilde = inversion on K)
    X X^(-1) = (s s~ + q t t~) delta  +  (t s~ + eps s t~) chi  +  (s e~ + e s~ + q e e~ - t t~) 1.
delta, chi, 1 are independent, so for the b = 1 target  T = (4m + 4) delta_G - 4 1_G  (m = q|K|):
    (I)  sum_i (s_i s_i~ + q t_i t_i~)              = (4m + 4) delta_K
    (Q)  sum_i (t_i s_i~ + eps s_i t_i~)             = 0
    (J)  sum_i (s_i e_i~ + e_i s_i~ + q e_i e_i~ - t_i t_i~) = -4 1_K
Stage 1 (the F_q-compression, = (I) + q (J)):  u = s + q e = x0 + (q - 1) e,  sum_i PAF_{u_i} = (4m+4) delta - 4q 1.
Defect identity [P]: (I) at 0:  (q + 1) e_plus + (q - 3) e_minus = 4(|K| - 1)   (e_plus: x0 = e, e_minus: x0 = -e).
Trivial character of K (= the Z_q-level system, the energy identity of report/h4892_answer.md section 5) [P]:
    sum r_i^2 = 4 (r_i = s^_i + q e^_i, the block row sums),  sum (s^_i^2 + q t^_i^2) = 4m + 4,  sum s^_i t^_i = 0.
Any real character psi of K gives the same with psi-sums and right-hand side 4m + 4 (no row condition).

usage: python3 theory/h2396_plugin_system.py check      (random non-solutions: G-level PAF = K-level formula)
       python3 theory/h2396_plugin_system.py identity    (defect allocations for q = 13, |K| = 2..60)"""
import sys, itertools, json, random
import numpy as np


def legendre(q):
    sq = {(x * x) % q for x in range(1, q)}
    return np.array([0] + [1 if y in sq else -1 for y in range(1, q)], dtype=np.int64)


def k_elements(axes):
    return list(itertools.product(*[range(n) for n in axes]))


def g_function(axes, q, x0, x1, x2):
    """the function on G = K x F_q (K = prod Z_axes, elements in itertools.product order), as a dict"""
    chi = legendre(q); els = k_elements(axes); f = {}
    for i, k in enumerate(els):
        for z in range(q):
            f[(k, z)] = int(x0[i]) if z == 0 else (int(x1[i]) if chi[z] == 1 else int(x2[i]))
    return f


def g_paf(axes, q, f):
    """direct G-level periodic autocorrelation, by the definition (double sum)"""
    els = k_elements(axes); out = {}
    for h in els:
        for w in range(q):
            acc = 0
            for k in els:
                kh = tuple((a + b) % n for a, b, n in zip(k, h, axes))
                for z in range(q):
                    acc += f[(k, z)] * f[(kh, (z + w) % q)]
            out[(h, w)] = acc
    return out


def k_conv_inv(axes, a, b):
    """(a b~)(h) = sum_k a(k + h) b(k) on K"""
    els = k_elements(axes); idx = {e: i for i, e in enumerate(els)}; out = np.zeros(len(els), np.int64)
    for hi, h in enumerate(els):
        s = 0
        for ki, k in enumerate(els):
            kh = idx[tuple((x + y) % n for x, y, n in zip(k, h, axes))]
            s += a[kh] * b[ki]
        out[hi] = s
    return out


def k_level(axes, q, x0, x1, x2):
    x0, x1, x2 = (np.asarray(v, np.int64) for v in (x0, x1, x2))
    e = (x1 + x2) // 2; t = (x1 - x2) // 2; s = x0 - e
    eps = int(legendre(q)[q - 1])
    I = k_conv_inv(axes, s, s) + q * k_conv_inv(axes, t, t)
    Qp = k_conv_inv(axes, t, s) + eps * k_conv_inv(axes, s, t)
    J = k_conv_inv(axes, s, e) + k_conv_inv(axes, e, s) + q * k_conv_inv(axes, e, e) - k_conv_inv(axes, t, t)
    return I, Qp, J


def check(trials=3, seed=1):
    rng = random.Random(seed); n_ok = 0; n = 0
    cases = [((2,), 5), ((3,), 5), ((2, 2), 5), ((4,), 13), ((6,), 13), ((3, 2), 13), ((5,), 17), ((2,), 3),
             ((3,), 7), ((4,), 11), ((2, 3), 7)]
    for axes, q in cases:
        chi = legendre(q); els = k_elements(axes); nK = len(els)
        for _ in range(trials):
            x0, x1, x2 = ([rng.choice((1, -1)) for _ in range(nK)] for _ in range(3))
            f = g_function(axes, q, x0, x1, x2); P = g_paf(axes, q, f)
            I, Qp, J = k_level(axes, q, x0, x1, x2)
            ok = all(P[(h, w)] == (I[hi] if w == 0 else 0) + (Qp[hi] * chi[w]) + J[hi]
                     for hi, h in enumerate(els) for w in range(q))
            n += 1; n_ok += ok
            # stage-1 form: the F_q compression of P equals PAF of u = x0 + (q-1) e
            e = (np.array(x1) + np.array(x2)) // 2; u = np.array(x0) + (q - 1) * e
            Pu = k_conv_inv(axes, u, u)
            ok1 = all(sum(P[(h, w)] for w in range(q)) == Pu[hi] for hi, h in enumerate(els))
            n += 1; n_ok += ok1
    return dict(checks=n, agree=n_ok)


def allocations(q, nK):
    """(e_plus, e_minus) with (q + 1) e_plus + (q - 3) e_minus = 4(|K| - 1)"""
    R = 4 * (nK - 1)
    return [(a, (R - (q + 1) * a) // (q - 3)) for a in range(R // (q + 1) + 1) if (R - (q + 1) * a) % (q - 3) == 0]


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'check'
    if cmd == 'check':
        print(json.dumps(check()))
    elif cmd == 'identity':
        q = int(sys.argv[2]) if len(sys.argv) > 2 else 13
        for nK in range(2, 61):
            al = allocations(q, nK)
            print(json.dumps(dict(q=q, K=nK, m=q * nK, v=q * nK + 1, allocations=al,
                                  defects=[a + b for a, b in al])))
