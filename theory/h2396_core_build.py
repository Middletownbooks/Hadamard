r"""Build the four +-1 sequences of a Paley plug-in Legendre-quadruple core on G = Z_n x F_q from a stage-1 solution u
(4 lines) and a completion t (one "T ..." line of theory/h2396_stage2.c), map G to cyclic Z_{nq} by the CRT
(g <-> (g mod n, g mod q)), check the core directly (sum PAF = 4m at 0, -4 elsewhere), and write it in the repository's
b = 1 core format (header "v=.. N=.. b=1 m=..", lines A=, B=, C=, D= of +/-) for scripts/verify_h2092.py.
Decoding [P]: u = +-1 -> (s, e) = (u, 0); u = +-q -> (0, u/q); u = +-(q-2) -> (-2u/(q-2), u/(q-2)); then
x0 = s + e, x1 = e + t, x2 = e - t; the block is x0 on Z_n x {0}, x1 on the squares, x2 on the non-squares.
usage: python3 theory/h2396_core_build.py n q UFILE TFILE OUT"""
import sys
import numpy as np

def main(n, q, ufile, tfile, out):
    U = [list(map(int, l.split())) for l in open(ufile) if l.strip()]
    tl = next(l for l in open(tfile) if l.startswith('T '))
    T = [list(map(int, part.split())) for part in tl[2:].split('|')]
    sq = {(z * z) % q for z in range(1, q)}
    m = n * q; seqs = []
    for i in range(4):
        x0, x1, x2 = [], [], []
        for k in range(n):
            v = U[i][k]
            if abs(v) == 1: s, e = v, 0
            elif abs(v) == q: s, e = 0, v // q
            elif abs(v) == q - 2: s, e = -2 * (v // (q - 2)), v // (q - 2)
            else: raise ValueError(v)
            t = T[i][k]
            assert (t == 0) == (e != 0), "t must vanish exactly at the defects"
            x0.append(s + e); x1.append(e + t); x2.append(e - t)
        assert all(abs(a) == 1 for a in x0 + x1 + x2)
        seq = np.zeros(m, np.int64)
        for g in range(m):
            k, z = g % n, g % q
            seq[g] = x0[k] if z == 0 else (x1[k] if z in sq else x2[k])
        seqs.append(seq)
    P = sum(np.array([int((x * np.roll(x, -h)).sum()) for h in range(m)]) for x in seqs)
    assert P[0] == 4 * m and (P[1:] == -4).all(), "not a Legendre-quadruple core"
    rows = [int(x.sum()) for x in seqs]
    with open(out, 'w') as f:
        f.write(f"v={m + 1} N={4 * (m + 1)} b=1 m={m}  Paley plug-in on Z{n} x F{q} (CRT), rows {rows}; theory/h2396_core_build.py\n")
        for name, x in zip("ABCD", seqs): f.write(name + "=" + "".join('+' if a > 0 else '-' for a in x) + "\n")
    print(dict(m=m, v=m + 1, rows=rows, paf_ok=True, out=out))

if __name__ == "__main__":
    n, q = int(sys.argv[1]), int(sys.argv[2])
    main(n, q, sys.argv[3], sys.argv[4], sys.argv[5])
