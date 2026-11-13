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
