r"""[E] Expected number of cores inside the fibre of a Z43-compressed solution at v = 863 (S7 host, G = Z2^2 x Z5 x Z43,
multiplier x4 of order 7 on Z43, k4 bordered family, t = 3), by the repository's lattice Gaussian
(`bt_gates.stratified_bits`: blocks sampled uniformly in their fibres, remaining class equations read at the true target).
Input: a certificate in the format of `data/h3452_stage1.json` (codex branch, `theory/h3452_compress.py`).
Measured 2026-11-03 on that certificate: log2 fibre 339.11, bits 275.35 / 275.62 (N = 4000 / 16000), rank 59,
Mahalanobis 0.18 / 0.15, so log2 E(fibre) = 63.8 / 63.5.
usage: python3 theory/h863_fibre_estimate.py CERT.json [N ...]"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_engine as E, bt_gates as BG

def main(path, Ns=(4000, 16000)):
    d = json.load(open(path))
    M = np.eye(4, dtype=int); M[-1, -1] = 4
    st = E.BStruct((2, 2, 5, 43), d['spec']['Kgens'], 'k', [M.tolist()])
    qt = st.quotient(st.G.subgroup(d['Ngens']))
    ys = [np.array(y) for y in d['y']]
    lm = sum(BG.FibreSampler(st, i, qt, ys[i]).log2_size for i in range(4))
    for N in Ns:
        r = BG.stratified_bits(st, qt, ys, N=N, seed=7)
        print(json.dumps(dict(N=N, ncls=st.ncls, r=st.r[0], log2_fibre=round(lm, 2), bits=r['bits'],
                              rank=r.get('rank'), mahal=r.get('mahal'), log2E_fibre=round(lm - r['bits'], 2))), flush=True)

if __name__ == '__main__':
    main(sys.argv[1], [int(x) for x in sys.argv[2:]] or (4000, 16000))
