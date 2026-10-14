r"""Finish a 647 x353 core: write it, complete it, verify it, and build the publication bundle pub/H2588/ [M].

1. Read logs/p647/CORE_*.json (four +-1 sequences on Z_646, the order of the slice's blocks).  Re-check sum_k PAF_k
   = (2584, -4, ..., -4) and the multiplier invariance (x_k(353 g) = x_k(g)).
2. Put the block with row sum +-2 first and negate all four if that row is -2 (PAF unchanged), so block A has
   row +2 and the border of pub/H2292/border.json applies (corner = -S T^T / 2).
3. Write data/core_v647_b1_p647_<tag>.txt (A=, B=, C=, D= lines of + and -).
4. Internal completion check: theory/nd_doubly.bordered_gs_cyclic (H H^T = N I in int64).
5. Independent verification: scripts/verify_h2588.py (no code shared with theory/), which also writes
   data/H2588_packed.npz and prints the SHA-256 of the packed rows.
6. pub/H2588/: sequences_Z646.txt, border.json, H2588_normalized.txt.gz (first row and column all +1), verify.py
   (the standalone checker of pub/H2092, pointed at the 2588 files).
usage: python3 theory/p647_finish.py logs/p647/CORE_xxx.json [--tag TAG]"""
import sys, os, json, gzip, hashlib, subprocess, shutil
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def paf(x):
    return np.rint(np.fft.ifft(np.abs(np.fft.fft(x.astype(float))) ** 2).real).astype(np.int64)

def main(core_json, tag=None, log=print):
    d = json.load(open(core_json)); X = np.array(d['core'], np.int64)
    m = X.shape[1]; assert X.shape == (4, 646) and set(np.unique(X).tolist()) <= {-1, 1}
    P = sum(paf(x) for x in X)
    assert P[0] == 4 * m and (P[1:] == -4).all(), "not a core"
    g = np.arange(m)
    assert all((X[k][(353 * g) % m] == X[k]).all() for k in range(4)), "not invariant under x353"
    rows = X.sum(1)
    k2 = int(np.nonzero(np.abs(rows) == 2)[0][0]); order = [k2] + [k for k in range(4) if k != k2]
    X = X[order]
    if X[0].sum() < 0: X = -X
    assert X[0].sum() == 2 and (X[1:].sum(1) == 0).all()
    tag = tag or os.path.splitext(os.path.basename(core_json))[0].replace("CORE_", "")
    core_txt = os.path.join(ROOT, "data", f"core_v647_b1_p647_{tag}.txt")
    with open(core_txt, "w") as f:
        f.write(f"# b = 1 core for H(2588): four +-1 sequences on cyclic Z_646, invariant under x353 (647 x353 slice)\n")
        f.write(f"# source: {os.path.relpath(core_json, ROOT)}; block order {order} of the slice, row sums (2, 0, 0, 0)\n")
        for name, x in zip("ABCD", X): f.write(name + "=" + "".join("+" if v > 0 else "-" for v in x) + "\n")
    log(json.dumps(dict(stage="core_written", path=os.path.relpath(core_txt, ROOT), rows=X.sum(1).tolist())))
    # internal completion check
    from nd_doubly import bordered_gs_cyclic
    H = bordered_gs_cyclic(X)
    log(json.dumps(dict(stage="internal_bordered_gs", hadamard=H is not None, N=None if H is None else int(H.shape[0]))))
    # independent verification
    out = os.path.join(ROOT, "data", "H2588_packed.npz")
    r = subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "verify_h2588.py"), core_txt, out], capture_output=True, text=True)
    log(r.stdout.strip()); assert r.returncode == 0, r.stderr
    res = json.loads(r.stdout.strip().splitlines()[-1]); assert res['hadamard']
    # publication bundle
    pub = os.path.join(ROOT, "pub", "H2588"); os.makedirs(pub, exist_ok=True)
    with open(os.path.join(pub, "sequences_Z646.txt"), "w") as f:
        f.write("# H(2588): b = 1 bordered Goethals-Seidel core, four +-1 sequences on cyclic Z_646 (v = 647),\n")
        f.write("# each invariant under the multiplier 353 (353 = 13 mod 17 of order 4, 353 = 11 mod 19 of order 3).\n")
        f.write("# sum of periodic autocorrelations = -4 at every nonzero shift; row sums 2, 0, 0, 0.\n")
        for name, x in zip("ABCD", X): f.write(name + "=" + "".join("+" if v > 0 else "-" for v in x) + "\n")
    json.dump(dict(S=res['S'], T=res['T'], corner=res['corner']), open(os.path.join(pub, "border.json"), "w"))
    bits = np.load(out)['bits']; N = int(np.load(out)['N'])
    Hm = np.where(np.unpackbits(bits, axis=1)[:, :N] > 0, 1, -1).astype(np.int64)
    Hn = Hm * Hm[:, [0]]; Hn = Hn * Hn[[0], :]
    with gzip.open(os.path.join(pub, "H2588_normalized.txt.gz"), "wt") as f:
        for row in Hn: f.write("".join("+" if v > 0 else "-" for v in row) + "\n")
    src = open(os.path.join(ROOT, "pub", "H2092", "verify.py")).read()
    src = src.replace("H(2092)", "H(2588)").replace("sequences_Z522.txt", "sequences_Z646.txt").replace("H2092_normalized", "H2588_normalized")
    open(os.path.join(pub, "verify.py"), "w").write(src)
    log(json.dumps(dict(stage="pub_bundle", dir=os.path.relpath(pub, ROOT), files=sorted(os.listdir(pub)),
                        sha256_of_packed_rows=res['sha256_of_packed_rows'])))
    return res

if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], a[a.index("--tag") + 1] if "--tag" in a else None, log=lambda s: print(s, flush=True))
