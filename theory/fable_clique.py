#!/usr/bin/env python3
"""k-cliques of a graph given by a boolean adjacency matrix: a small C routine (bitsets, pool-size bound), compiled
on first use with gcc.  Used by fable_border14.py and fable_border15.py.  Exhaustive unless `truncated` is returned."""
import os, subprocess, ctypes
import numpy as np

SRC = r'''
#include <stdint.h>
#include <string.h>
static int W, need, maxout, nout, trunc_; static uint64_t *adj; static long long nodes, maxnodes; static int *out;
static int cur[128];
static void rec(int depth, uint64_t *P){
  nodes++;
  if (depth == need){ if (nout < maxout){ memcpy(out + (long)nout*need, cur, need*sizeof(int)); nout++; } else trunc_ = 1; return; }
  int cnt = 0; for (int w = 0; w < W; w++) cnt += __builtin_popcountll(P[w]);
  if (cnt < need - depth) return;
  if (nodes > maxnodes){ trunc_ = 1; return; }
  uint64_t *P2 = P + W;
  for (int w = 0; w < W; w++){ uint64_t bits = P[w];
    while (bits){
      if (cnt < need - depth) return;
      int b = __builtin_ctzll(bits); bits &= bits - 1; cnt--;
      long v = (long)w*64 + b; int c2 = 0;
      for (int u = 0; u < w; u++) P2[u] = 0;
      P2[w] = bits & adj[v*W + w]; c2 += __builtin_popcountll(P2[w]);
      for (int u = w + 1; u < W; u++){ P2[u] = P[u] & adj[v*W + u]; c2 += __builtin_popcountll(P2[u]); }
      if (c2 < need - depth - 1) continue;
      cur[depth] = (int)v; rec(depth + 1, P2);
      if (trunc_) return;
    }
  }
}
/* ---- coupled search: a need_z-clique Z in graph adjZ such that the set E(Z) = AND of cm[z] over z in Z (a bitset
   over ne vertices) contains a need_e-clique of graph adjE.  Stops at the first solution.  ---- */
static int Wz, We, needz, neede, solz[128], sole[128], foundb; static uint64_t *adjZ, *adjE, *cmz; static long long nodesb, maxnodesb, zcl; static int truncb;
static int cure[128]; static int collect, maxcol, ncol; static int *outall;
static int rec_e(int depth, uint64_t *P){
  nodesb++;
  if (depth == neede){ memcpy(sole, cure, neede*sizeof(int)); return 1; }
  int cnt = 0; for (int w = 0; w < We; w++) cnt += __builtin_popcountll(P[w]);
  if (cnt < neede - depth) return 0;
  if (nodesb > maxnodesb){ truncb = 1; return 0; }
  uint64_t *P2 = P + We;
  for (int w = 0; w < We; w++){ uint64_t bits = P[w];
    while (bits){
      if (cnt < neede - depth) return 0;
      int b = __builtin_ctzll(bits); bits &= bits - 1; cnt--;
      long v = (long)w*64 + b; int c2 = 0;
      for (int u = 0; u < w; u++) P2[u] = 0;
      P2[w] = bits & adjE[v*We + w]; c2 += __builtin_popcountll(P2[w]);
      for (int u = w + 1; u < We; u++){ P2[u] = P[u] & adjE[v*We + u]; c2 += __builtin_popcountll(P2[u]); }
      if (c2 < neede - depth - 1) continue;
      cure[depth] = (int)v; if (rec_e(depth + 1, P2)) return 1;
      if (truncb) return 0;
    }
  }
  return 0;
}
static uint64_t *ework;
static int rec_z(int depth, uint64_t *P, uint64_t *E){
  nodesb++;
  if (depth == needz){ zcl++; int ok = (neede == 0) ? 1 : rec_e(0, E);
    if (!ok) return 0;
    if (!collect) return 1;
    if (ncol < maxcol){ memcpy(outall + (long)ncol*needz, solz, needz*sizeof(int)); ncol++; } else truncb = 1;
    return 0; }
  int cnt = 0; for (int w = 0; w < Wz; w++) cnt += __builtin_popcountll(P[w]);
  if (cnt < needz - depth) return 0;
  if (nodesb > maxnodesb){ truncb = 1; return 0; }
  uint64_t *P2 = P + Wz; uint64_t *E2 = E + We;
  for (int w = 0; w < Wz; w++){ uint64_t bits = P[w];
    while (bits){
      if (cnt < needz - depth) return 0;
      int b = __builtin_ctzll(bits); bits &= bits - 1; cnt--;
      long v = (long)w*64 + b; int c2 = 0, ce = 0;
      for (int u = 0; u < We; u++){ E2[u] = E[u] & cmz[v*We + u]; ce += __builtin_popcountll(E2[u]); }
      if (ce < neede) continue;
      for (int u = 0; u < w; u++) P2[u] = 0;
      P2[w] = bits & adjZ[v*Wz + w]; c2 += __builtin_popcountll(P2[w]);
      for (int u = w + 1; u < Wz; u++){ P2[u] = P[u] & adjZ[v*Wz + u]; c2 += __builtin_popcountll(P2[u]); }
      if (c2 < needz - depth - 1) continue;
      solz[depth] = (int)v; if (rec_z(depth + 1, P2, E2)) return 1;
      if (truncb) return 0;
    }
  }
  return 0;
}
long long coupled(int nz, int ne, int needz_, int neede_, uint64_t *adjZ_, uint64_t *cm_, uint64_t *adjE_,
                  uint64_t *startz, uint64_t *starte, uint64_t *workz, uint64_t *worke, long long maxnodes_,
                  int *outz, int *oute, int *found, int *trunc, long long *zcliques,
                  int collect_, int *outall_, int maxcol_, int *ncol_){
  collect = collect_; outall = outall_; maxcol = maxcol_; ncol = 0;
  Wz = (nz + 63)/64; We = (ne + 63)/64; if (We == 0) We = 1; needz = needz_; neede = neede_; adjZ = adjZ_; adjE = adjE_; cmz = cm_;
  nodesb = 0; maxnodesb = maxnodes_; truncb = 0; zcl = 0;
  memcpy(workz, startz, Wz*sizeof(uint64_t)); memcpy(worke, starte, We*sizeof(uint64_t));
  int f = rec_z(0, workz, worke);
  if (f){ memcpy(outz, solz, needz*sizeof(int)); memcpy(oute, sole, neede*sizeof(int)); }
  *found = f; *trunc = truncb; *zcliques = zcl; *ncol_ = ncol; return nodesb;
}
long long kcliques(int n, int need_, uint64_t *adj_, uint64_t *start, uint64_t *work, int *out_, int maxout_,
                   long long maxnodes_, int *nout_, int *trunc){
  W = (n + 63)/64; need = need_; adj = adj_; out = out_; maxout = maxout_; maxnodes = maxnodes_; nodes = 0; nout = 0; trunc_ = 0;
  memcpy(work, start, W*sizeof(uint64_t));
  if (need == 0){ nout = 1; } else rec(0, work);
  *nout_ = nout; *trunc = trunc_; return nodes;
}
'''
_lib = None
def _load():
    global _lib
    if _lib is None:
        d = os.path.join(os.path.dirname(os.path.abspath(__file__)), '__pycache__'); os.makedirs(d, exist_ok=True)
        c = os.path.join(d, 'fable_clique.c'); so = os.path.join(d, 'fable_clique.so')
        if not os.path.exists(so) or open(c).read() != SRC if os.path.exists(c) else True:
            open(c, 'w').write(SRC)
            subprocess.check_call(['gcc', '-O2', '-shared', '-fPIC', '-o', so, c])
        _lib = ctypes.CDLL(so)
        _lib.kcliques.restype = ctypes.c_longlong
    return _lib

def kcliques(A, need, cand=None, maxout=100000, maxnodes=10**10):
    """A: n x n bool (symmetric, zero diagonal); cand: iterable of allowed vertices (default all).
    Returns (list of cliques as sorted index lists, nodes, truncated)."""
    lib = _load(); n = A.shape[0]
    if need == 0: return [[]], 0, False
    if n == 0: return [], 0, False
    W = (n + 63) // 64
    Ab = np.zeros((n, W * 64), dtype=np.uint8); Ab[:, :n] = A
    adj = np.packbits(Ab.reshape(n, W, 64), axis=2, bitorder='little').view(np.uint64).reshape(n, W).copy()
    m = np.zeros(W * 64, dtype=np.uint8)
    if cand is None: m[:n] = 1
    else: m[np.fromiter(cand, dtype=np.int64)] = 1
    start = np.packbits(m.reshape(W, 64), axis=1, bitorder='little').view(np.uint64).reshape(W).copy()
    work = np.zeros((need + 2) * W, dtype=np.uint64)
    out = np.zeros(max(1, maxout) * need, dtype=np.int32)
    nout = ctypes.c_int(0); tr = ctypes.c_int(0)
    P = ctypes.POINTER
    nodes = lib.kcliques(ctypes.c_int(n), ctypes.c_int(need), adj.ctypes.data_as(P(ctypes.c_uint64)),
                         start.ctypes.data_as(P(ctypes.c_uint64)), work.ctypes.data_as(P(ctypes.c_uint64)),
                         out.ctypes.data_as(P(ctypes.c_int)), ctypes.c_int(maxout), ctypes.c_longlong(maxnodes),
                         ctypes.byref(nout), ctypes.byref(tr))
    res = out[:nout.value * need].reshape(nout.value, need).tolist()
    return res, int(nodes), bool(tr.value)

def _bits(A, ncols):
    n = A.shape[0]; W = max(1, (ncols + 63) // 64)
    Ab = np.zeros((n, W * 64), dtype=np.uint8); Ab[:, :ncols] = A
    return np.packbits(Ab.reshape(n, W, 64), axis=2, bitorder='little').view(np.uint64).reshape(n, W).copy(), W

def coupled(AZ, CM, AE, needz, neede, candz=None, cande=None, maxnodes=10**11, collect=0):
    """A needz-clique Z of graph AZ (nz x nz bool) such that the common neighbours of Z in the bipartite matrix CM
    (nz x ne bool; inside cande) contain a neede-clique of graph AE (ne x ne bool).  First solution only.
    Returns (zlist or None, elist, nodes, truncated, number of complete Z examined)."""
    lib = _load(); nz, ne = CM.shape
    adjZ, Wz = _bits(AZ, nz); cm, We = _bits(CM, ne); adjE, We2 = _bits(AE, ne) if ne else (np.zeros((1, 1), dtype=np.uint64), 1)
    mz = np.zeros((1, nz), dtype=bool); mz[0, (list(candz) if candz is not None else slice(None))] = True
    me = np.zeros((1, max(ne, 1)), dtype=bool)
    if ne: me[0, (list(cande) if cande is not None else slice(None))] = True
    sz, _ = _bits(mz, nz); se, _ = _bits(me[:, :max(ne, 1)], max(ne, 1))
    workz = np.zeros((needz + 2) * Wz, dtype=np.uint64); worke = np.zeros((needz + neede + 4) * We, dtype=np.uint64)
    outz = np.zeros(max(needz, 1), dtype=np.int32); oute = np.zeros(max(neede, 1), dtype=np.int32)
    f = ctypes.c_int(0); tr = ctypes.c_int(0); zc = ctypes.c_longlong(0); ncol = ctypes.c_int(0)
    outall = np.zeros(max(1, collect) * max(needz, 1), dtype=np.int32)
    P = ctypes.POINTER; u64 = lambda a: a.ctypes.data_as(P(ctypes.c_uint64)); i32 = lambda a: a.ctypes.data_as(P(ctypes.c_int))
    lib.coupled.restype = ctypes.c_longlong
    nodes = lib.coupled(ctypes.c_int(nz), ctypes.c_int(ne), ctypes.c_int(needz), ctypes.c_int(neede), u64(adjZ), u64(cm), u64(adjE),
                        u64(sz), u64(se), u64(workz), u64(worke), ctypes.c_longlong(maxnodes), i32(outz), i32(oute),
                        ctypes.byref(f), ctypes.byref(tr), ctypes.byref(zc),
                        ctypes.c_int(1 if collect else 0), i32(outall), ctypes.c_int(collect), ctypes.byref(ncol))
    if collect:      # every Z-clique whose common neighbourhood contains an E-clique (up to `collect` of them)
        return outall[:ncol.value * needz].reshape(ncol.value, needz).tolist(), None, int(nodes), bool(tr.value), int(zc.value)
    if f.value: return outz[:needz].tolist(), oute[:neede].tolist(), int(nodes), bool(tr.value), int(zc.value)
    return None, None, int(nodes), bool(tr.value), int(zc.value)

if __name__ == '__main__':     # self-test against a plain Python enumeration
    import itertools
    rng = np.random.default_rng(0)
    for n, p, k in ((40, 0.5, 5), (70, 0.3, 4), (130, 0.2, 3), (25, 0.8, 7)):
        A = rng.random((n, n)) < p; A = np.triu(A, 1); A = A | A.T
        ref = sorted(list(c) for c in itertools.combinations(range(n), k) if all(A[i, j] for i, j in itertools.combinations(c, 2)))
        got, nodes, tr = kcliques(A, k, maxout=10**6)
        assert not tr and sorted(got) == ref, (n, k, len(got), len(ref))
        cand = list(range(0, n, 2)); ref2 = [c for c in ref if all(v % 2 == 0 for v in c)]
        got2, _, _ = kcliques(A, k, cand=cand, maxout=10**6); assert sorted(got2) == ref2
        print(n, k, len(ref), 'cliques: C routine agrees', flush=True)
    for trial in range(30):
        nz, ne, kz, ke = 30, 26, 3, 3
        AZ = rng.random((nz, nz)) < 0.4; AZ = np.triu(AZ, 1); AZ = AZ | AZ.T
        AE = rng.random((ne, ne)) < 0.4; AE = np.triu(AE, 1); AE = AE | AE.T
        CM = rng.random((nz, ne)) < (0.45 if trial % 2 else 0.6)
        ref = False
        for z in itertools.combinations(range(nz), kz):
            if not all(AZ[i, j] for i, j in itertools.combinations(z, 2)): continue
            com = [e for e in range(ne) if all(CM[i, e] for i in z)]
            if any(all(AE[i, j] for i, j in itertools.combinations(c, 2)) for c in itertools.combinations(com, ke)): ref = True; break
        zs, es, nodes, tr, zc = coupled(AZ, CM, AE, kz, ke)
        assert (zs is not None) == ref and not tr, (trial, ref, zs)
        if zs is not None:
            assert all(AZ[i, j] for i, j in itertools.combinations(zs, 2)) and all(AE[i, j] for i, j in itertools.combinations(es, 2)) and all(CM[i, e] for i in zs for e in es)
    print('coupled routine agrees with brute force on 30 random instances', flush=True)
    for trial in range(10):
        nz, ne, kz, ke = 24, 20, 3, 2
        AZ = rng.random((nz, nz)) < 0.4; AZ = np.triu(AZ, 1); AZ = AZ | AZ.T
        AE = rng.random((ne, ne)) < 0.4; AE = np.triu(AE, 1); AE = AE | AE.T
        CM = rng.random((nz, ne)) < 0.5
        ref = []
        for z in itertools.combinations(range(nz), kz):
            if not all(AZ[i, j] for i, j in itertools.combinations(z, 2)): continue
            com = [e for e in range(ne) if all(CM[i, e] for i in z)]
            if any(all(AE[i, j] for i, j in itertools.combinations(c, 2)) for c in itertools.combinations(com, ke)): ref.append(list(z))
        got, _, nodes, tr, zc = coupled(AZ, CM, AE, kz, ke, collect=100000)
        assert sorted(got) == sorted(ref) and not tr
    print('collect mode agrees with brute force on 10 random instances', flush=True)
