r"""[M] Independent verifier for G/Z5-level solutions at v = 863 (shares no code with theory/h863_lib.py, h863_g5_solve.py
or the Codex scripts; pure Python, no numpy).

The certificate data/h3452_stage1.json (Z43-compressed solution y) is by OpenAI Codex, directed by the repository owner.

For each solution w (4 blocks x 4 Z2^2-cells x 6 classes) it
  1. rebuilds the cyclotomic classes of <4> in Z43^* from scratch;
  2. derives e and the class-sign row sums from y (y = e + 7 * sum of six class signs, e = +-1);
  3. constructs an explicit 5 x 6 sign table per (block, cell) with the required row sums and column sums w
     (Gale-Ryser greedy: each row puts its +1 entries in the columns with the largest remaining demand), and checks it;
  4. assembles four explicit +-1 sequences on G = Z2 x Z2 x Z5 x Z43 (length 860), S-invariant by construction;
  5. checks the Z43-compression of each sequence against y (all 20 cells), and the summed PAF of the four
     Z5-compressions (length 172 on Z2 x Z2 x Z43) on all 172 lags against 3392 / -60 / 20.
usage: python3 theory/h863_verify_g5.py SOLUTIONS.json [CERT.json]   (SOLUTIONS.json: {"w": [w, ...]})
"""
import json, sys


def classes43():
    sub = []
    a = 1
    for _ in range(7):
        sub.append(a); a = a * 4 % 43
    assert len(set(sub)) == 7 and a == 1
    cls = {}
    reps = []
    for t in range(1, 43):
        if t not in cls:
            for h in sub: cls[t * h % 43] = len(reps)
            reps.append(t)
    return cls, reps


def build_table(rowsums, colsums):
    """explicit 5 x 6 +-1 table with the given row sums and column sums, or None"""
    R, C = len(rowsums), len(colsums)
    rp = [(s + C) // 2 for s in rowsums]            # +1 entries per row
    cp = [(s + R) // 2 for s in colsums]            # +1 entries per column
    if sum(rp) != sum(cp) or any((s + C) % 2 for s in rowsums) or any((s + R) % 2 for s in colsums):
        return None
    rem = cp[:]
    T = [[-1] * C for _ in range(R)]
    for r in sorted(range(R), key=lambda r: -rp[r]):
        order = sorted(range(C), key=lambda c: -rem[c])
        for c in order[:rp[r]]:
            if rem[c] <= 0: return None
            T[r][c] = 1; rem[c] -= 1
    if any(rem): return None
    for r in range(R): assert sum(T[r]) == rowsums[r]
    for c in range(C): assert sum(T[r][c] for r in range(R)) == colsums[c]
    return T


def main(solpath, certpath='data/h3452_stage1.json'):
    cert = json.load(open(certpath))
    Y = cert['y']
    cls, reps = classes43()
    sols = json.load(open(solpath))['w']
    results = []
    for n, w in enumerate(sols):
        X = []                                         # four sequences, index ((q1*2+q2)*5+u)*43 + t
        ok_tables = True
        for i in range(4):
            x = [0] * 860
            for q in range(4):
                rows = []; es = []
                for u in range(5):
                    yv = Y[i][q * 5 + u]
                    e = 1 if yv % 14 == 1 else -1
                    assert (yv - e) % 14 == 0
                    es.append(e); rows.append((yv - e) // 7)
                T = build_table(rows, w[i][q])
                if T is None:
                    ok_tables = False; break
                for u in range(5):
                    base = (q * 5 + u) * 43
                    x[base] = es[u]
                    for t in range(1, 43): x[base + t] = T[u][cls[t]]
            if not ok_tables: break
            X.append(x)
        if not ok_tables:
            results.append(dict(n=n, ok=False, why='no table')); continue
        # 1. Z43-compression equals y
        ok_y = all(sum(X[i][(k * 43) + t] for t in range(43)) == Y[i][k] for i in range(4) for k in range(20))
        # 2. Z5-compression, summed PAF on Z2 x Z2 x Z43, all 172 lags
        Zs = []
        for i in range(4):
            z = [0] * 172
            for q in range(4):
                for t in range(43):
                    z[q * 43 + t] = sum(X[i][(q * 5 + u) * 43 + t] for u in range(5))
            Zs.append(z)
        bad = 0
        for qa in range(4):
            for t0 in range(43):
                v = 0
                for z in Zs:
                    for q in range(4):
                        for t in range(43):
                            v += z[q * 43 + t] * z[(q ^ qa) * 43 + (t + t0) % 43]
                tgt = 3392 if (qa == 0 and t0 == 0) else (-60 if qa == 0 else 20)
                bad += v != tgt
        # S-invariance of the explicit sequences (x(h, 4t) = x(h, t))
        ok_s = all(X[i][k * 43 + t] == X[i][k * 43 + (4 * t) % 43] for i in range(4) for k in range(20) for t in range(43))
        results.append(dict(n=n, ok=ok_y and bad == 0 and ok_s, y_ok=ok_y, s_invariant=ok_s, bad_lags=bad))
    print(json.dumps(dict(cert=certpath, solutions=solpath, checked=len(results),
                          all_ok=all(r['ok'] for r in results), results=results)))


if __name__ == '__main__':
    main(*sys.argv[1:])
