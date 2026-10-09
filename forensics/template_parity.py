"""Minimal number of deviations from a four-arc support pattern forced by parity alone (b=1, m = 2h).
For fixed supports S_k = {sigma_k = +1} on Z_h, PAF_{c_k}(t) = |S_k ∩ (S_k + t)| (mod 2), and the b=1 conditions need
   sum_k |S_k ∩ (S_k + t)| even for t = 1..n (target -2),   |S_A| odd, |S_B|, |S_C|, |S_D| even (row sums (+-1,0,0,0)).
(The negacyclic half then has the right parities automatically: |S^c ∩ (S^c+t)| = h - 2|S| + |S ∩ (S+t)|.)
Candidates per block: an arc of length L (L near h/2) with 0, 1 or 2 flipped positions.  We find the least total number
of flips over the four blocks that satisfies every parity condition."""
import sys, itertools
import numpy as np
def ov_parity(s, n):
    return np.array([int((s * np.roll(s, -t)).sum()) & 1 for t in range(1, n + 1)], dtype=np.uint8)
def key(vec):
    return np.packbits(vec).tobytes()
for h in (83, 89):
    n = (h - 1) // 2; Ls = range(n - 6, n + 8)
    cands = {0: [], 1: [], 2: []}     # (vector, |S| parity, description)
    for L in Ls:
        base = np.zeros(h, int); base[:L] = 1
        cands[0].append((ov_parity(base, n), L & 1, (L,)))
        for y in range(h):
            s = base.copy(); s[y] ^= 1
            cands[1].append((ov_parity(s, n), int(s.sum()) & 1, (L, y)))
    # 2 flips in one block: only needed if 1 total flip is impossible; build lazily
    def search(flips_per_block):
        """flips_per_block: tuple of 4 flip counts (block A first); returns one solution or None"""
        lists = [cands[f] for f in flips_per_block]
        # parity of |S|: A odd, others even
        want = [1, 0, 0, 0]
        L0 = [c for c in lists[0] if c[1] == want[0]]; L1 = [c for c in lists[1] if c[1] == want[1]]
        L2 = [c for c in lists[2] if c[1] == want[2]]; L3 = [c for c in lists[3] if c[1] == want[3]]
        right = {}
        for c2 in L2:
            for c3 in L3:
                right.setdefault(key(c2[0] ^ c3[0]), (c2, c3))
        for c0 in L0:
            for c1 in L1:
                hit = right.get(key(c0[0] ^ c1[0]))
                if hit: return (c0[2], c1[2], hit[0][2], hit[1][2])
        return None
    print(f'h={h} (v={2*h+1}), arc lengths {Ls.start}..{Ls.stop - 1}:')
    print('   0 flips (four exact arcs):', search((0, 0, 0, 0)))
    for pat in sorted(set(itertools.permutations((1, 0, 0, 0)))):
        r = search(pat); print('   1 flip  in block', 'ABCD'[pat.index(1)], ':', r)
    for pat in sorted(set(itertools.permutations((1, 1, 0, 0)))):
        r = search(pat); print('   2 flips in blocks', ''.join('ABCD'[i] for i in range(4) if pat[i]), ':', 'feasible, e.g. ' + str(r) if r else None)
    # which single flips can occur at all (with three exact arcs in the other blocks)?
    pure = {}
    for c in cands[0]: pure.setdefault(c[1], []).append(c)
    ok = set()
    for c in cands[1]:
        for blk in range(4):
            want = [1, 0, 0, 0]
            if c[1] != want[blk]: continue
            others = [want[i] for i in range(4) if i != blk]
            lists = [pure.get(w, []) for w in others]
            target = c[0]
            found = any(not (target ^ a[0] ^ b[0] ^ d[0]).any() for a in lists[0] for b in lists[1] for d in lists[2])
            if found: ok.add(c[2])
    desc = sorted({('centre of the arc' if y == (L - 1) / 2 else 'centre of the complement' if y - L == (h - L - 1) / 2 else f'other (L={L}, y={y})') for L, y in ok})
    print(f'   h={h}: positions of a lone parity-legal flip relative to its arc [0, L):', desc)
