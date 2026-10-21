r"""The k = 8 border problem (q = 8, t = 7, T = all nontrivial characters of Q) for monomial Q-level classes [P, M].

Monomial classes.  At the 751 survivors with log2 E_Q >= 0 (Q = Z_8) the Q-level classes are monomial:
z_k = c_k + a_k delta_{p_k} with {|a_k|} = {2, 2, 2, 4} (sum a_k^2 = 28 = 4t).  Translating, negating or permuting
the blocks of a core gives a core (the PAF sum is unchanged) whose compressions are translated, negated or permuted.
So a monomial class can be completed (after such a change of the core) iff some normal form, e.g.
z = (2, 2, 2, 4) delta_0 up to constants, has a border; the constants do not enter Mc on V_T [P].  The answer depends
only on Q.

1. Kronecker frames.  The frame H_4 (x) H_8[1:] (theory/bt_border.frames_balanced, first frame) has no admissible
   vector at Q = Z_2^3 [P]: on V_T the array is Q_S (x) I with Q_S the 4 x 4 GS matrix of a, so the entries of
   Rh Mc^T s are sum_o C[i, o] phi_o(l) with C = H_4 Q_S^T (C C^T = 112 I) and phi_o(l) in {0, +-4, +-8} the Walsh
   coefficients of block o of s.  Flatness forces phi(l) = (2/7) C^T e for a sign vector e, but
   C^T e = Q_S H_4 e has entries in {+-8, +-16} or {+-4, +-12, +-20}: never in {0, +-14, +-28}.  The block-permuted
   variants were searched [M] (none has an admissible vector, at Z_2^3 and at Z_8).
2. Q = Z_2^3.  The c = 8 Singer-type template of theory/ab_singer.py (the form of Alpoge's 1436 border;
   logs/ab/singer/c8_control) gives a border in the layout used here (Kumar's, theory/bt_border.gs_array) [M],
   checked exactly: X X^T = 28 I, Ch in {+-1}, Ch Ch^T = I_4 (x) Gamma_Q, X Ch^T + Rh Mc^T = 0.
   So every monomial class at Q = Z_2^3 has a border up to core equivalence.
3. Q = Z_8.  The same frames (both template frames, every block sign) have no admissible vector [M].  On V_T the
   Z_8 array is Q_S (x) Pi_+ + Q_S^T (x) Pi_- (Pi_+- = projections onto the balanced functions that are even / odd
   under negation), not Q_S (x) I.  No frame symmetry turns one into the other [P]: block-diagonal column
   permutations on the two sides (and block permutations and signs) map the 16 blocks to sign * sigma for one
   permutation sigma only if the product of the four blocks around every rectangle of the block array is the
   identity, but a rectangle through exactly one diagonal block (a translation) and three off-diagonal blocks
   (reflections g -> p - g) gives a reflection, never the identity on Z_8.  The diagonal symmetries of the Z_8
   array are the permutations commuting with negation (order 96, counted here [M]), so a Singer-type ansatz (a
   symmetry of order 7) is not available.  The border status of the Z_8 classes is open.
usage: python3 theory/bt_border_k8.py [--random 40] [--out logs/bt/border_k8.jsonl]
"""
import sys, os, json, time, itertools, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_border as BB

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATE = os.path.join(ROOT, 'logs/ab/singer/c8_control/template_c8_cyc_a2224_f0_1_layout.npz')
q, t = 8, 7


def tables(kind):
    I = np.arange(q)
    if kind == 'Z2^3':
        return I[:, None] ^ I[None, :], I.copy()
    return (I[:, None] + I[None, :]) % q, (-I) % q


def gamma_mat(add, neg):
    gam = np.full(q, -4, np.int64); gam[0] = 4 * t                 # Gamma = 4 sum_{chi != 1} chi
    return gam[add[neg[:, None], np.arange(q)[None, :]]]


def normal_form():
    z = np.ones((4, q), np.int64)
    z[:, 0] += np.array([2, 2, 2, 4])
    return z


def border_with(Mc, Gm, F, tmax=20.0):
    """border (X, Ch) for the frame F (Rh), checked exactly, or None; also |S_R|"""
    S = BB.admissible_balanced(Mc, F, t, q)
    if len(S) < 4 * t:
        return None, len(S)
    Xc = BB.clique(S, 4 * t, tmax=tmax)
    if Xc is None:
        return None, len(S)
    X = Xc.T.copy()
    num = -(Mc @ F.T @ X)
    if (num % (4 * t)).any():
        return None, len(S)
    Ch = num // (4 * t)
    ok = bool((X @ X.T == 4 * t * np.eye(4 * t, dtype=np.int64)).all() and (np.abs(Ch) == 1).all()
              and (Ch @ Ch.T == np.kron(np.eye(4, dtype=np.int64), Gm)).all()
              and (X @ Ch.T + F @ Mc.T == 0).all()
              and (F @ F.T == 4 * q * np.eye(4 * t, dtype=np.int64)).all()
              and all((F[:, J * q:(J + 1) * q].sum(1) == 0).all() for J in range(4)))
    return ((X, Ch) if ok else None), len(S)


def diagonal_symmetries(Mc):
    """permutations pi of Q with (pi on every block, both sides) fixing Mc on V_T (balanced vectors)"""
    P = np.eye(q, dtype=np.int64) - 1.0 / q                     # projection onto balanced vectors
    PV = np.kron(np.eye(4), P)
    M = PV @ Mc @ PV
    cnt = 0
    for pi in itertools.permutations(range(q)):
        Pm = np.zeros((q, q)); Pm[np.arange(q), pi] = 1
        B = np.kron(np.eye(4), Pm)
        if np.allclose(B @ M @ B.T, M):
            cnt += 1
    return cnt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--random', type=int, default=40)
    ap.add_argument('--out', default=os.path.join(ROOT, 'logs/bt/border_k8.jsonl'))
    a = ap.parse_args()
    T = np.load(TEMPLATE, allow_pickle=True)
    frames = dict(K=T['K'].astype(np.int64), Kc=T['Kc'].astype(np.int64))
    out = []
    z751 = [np.array([[1, 1, 1, 1, 3, 1, 1, 1]] * 3 + [[3, 7, 7, 7, 7, 7, 7, 7]], np.int64),
            np.array([[1, 1, 1, 1, 3, 1, 1, 1], [3, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 3, 1, 1, 1],
                      [3, 7, 7, 7, 7, 7, 7, 7]], np.int64)]           # the two classes of the 751 S10 r105 host
    cases = [('Z2^3', 'normal form', normal_form()), ('Z8', 'normal form', normal_form())]
    cases += [('Z8', '751 class %d' % i, z) for i, z in enumerate(z751)]
    for kind, label, z in cases:
        add, neg = tables(kind)
        Mc = BB.gs_array(list(z), add, neg); Gm = gamma_mat(add, neg)
        t0 = time.time()
        rec = dict(Q=kind, z=label, z_rows=[int(r) for r in z.sum(1)])
        # 1. Kronecker frames: the base frame and random block-permuted variants
        kron = [len(BB.admissible_balanced(Mc, F, t, q)) for F in BB.frames_balanced(q, count=a.random, seed=7)]
        rec.update(kronecker_frames=len(kron), kronecker_base_admissible=kron[0], kronecker_max_admissible=max(kron))
        # 2./3. the Singer-type template frames, every block sign
        best = 0; found = None
        for name, F0 in frames.items():
            for sg in itertools.product((1, -1), repeat=4):
                F = F0 * np.repeat(np.array(sg), q)[None, :]
                res, ns = border_with(Mc, Gm, F)
                best = max(best, ns)
                if res is not None and found is None:
                    found = dict(frame=name, block_signs=list(sg), admissible=ns)
        rec.update(template_max_admissible=best, template_border=found)
        if label == 'normal form':
            rec['diagonal_symmetries'] = diagonal_symmetries(Mc)
        rec['secs'] = round(time.time() - t0, 1)
        print(json.dumps(rec), flush=True)
        out.append(rec)
    with open(a.out, 'w') as fo:
        for r in out:
            fo.write(json.dumps(r) + '\n')


if __name__ == '__main__':
    main()
