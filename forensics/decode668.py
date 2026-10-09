"""Decode the announced order-668 matrix recorded in google-deepmind/formal-conjectures PR #4900, first commit 91ebb29
(FormalConjectures/Wikipedia/Hadamard668.lean). Rows are `packedRow a b c d` with BitVec concatenation a ++ b ++ c ++ d
(a = most significant); the least-significant bit is the first entry of the row; set bit = +1."""
import re, sys
import numpy as np
src = open('Hadamard668_announced.lean').read()
rows = re.findall(r'packedRow\s+(0x[0-9a-f]+)#167\s+(0x[0-9a-f]+)#167\s+(0x[0-9a-f]+)#167\s+(0x[0-9a-f]+)#167', src)
assert len(rows) == 668, len(rows)
H = np.zeros((668, 668), dtype=np.int64)
for r, (a, b, c, d) in enumerate(rows):
    val = (int(a, 16) << 501) | (int(b, 16) << 334) | (int(c, 16) << 167) | int(d, 16)
    bits = np.array([(val >> j) & 1 for j in range(668)])
    H[r] = 2 * bits - 1
ok = np.array_equal(H @ H.T, 668 * np.eye(668, dtype=np.int64))
print('decoded 668x668, H H^T = 668 I:', ok)
np.save('H668_announced.npy', H)
# normalised? first row / first column
print('first row all +1:', bool((H[0] == 1).all()), ' first column all +1:', bool((H[:, 0] == 1).all()))
print('corner 4x4:\n', H[:4, :4])
# border rows constant on the four 166-blocks of columns 4..667?
for r in range(4):
    blk = [set(H[r, 4 + 166 * q: 4 + 166 * (q + 1)].tolist()) for q in range(4)]
    print('border row', r, 'block values', blk)
for cidx in range(4):
    blk = [set(H[4 + 166 * q: 4 + 166 * (q + 1), cidx].tolist()) for q in range(4)]
    print('border col', cidx, 'block values', blk)
