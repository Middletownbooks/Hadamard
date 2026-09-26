# D-optimal transfer audit: registered predictions (written before any run)

Single thresholds; an unevaluable bar counts as FAIL.

P1  Exhaustive count of two-circulant pairs (A,B) with A A^T + B B^T = (2v-2) I + 2 J (all 2^(2v) ordered pairs, no
    normalisation) equals the prompt's figures exactly: v=7 -> 392, v=9 -> 1944, v=11 -> 0.   PASS iff all three equal.
P2  Every pair found in P1 gives R = [[A, B], [-B^T, A^T]] with R R^T = diag(L, L), L = (2v-2) I + 2 J, and
    |det R| = (2n-2)(n-2)^((n-2)/2) (Ehlich-Wojtas, n = 2v), checked in exact integer arithmetic.   PASS iff all.
P3  Negative tests: flipping one entry of a verified pair makes the P1 verifier AND the P2 matrix check fail; a
    deliberately wrong EW constant (n-1 instead of n-2 in the base) makes the det check fail.   PASS iff all fail.
P4  For odd v <= 13: PAF(s) = v (mod 4) for every +-1 sequence and lag, and sum_{s>=1} PAF(s) = a^2 - v.
    PASS iff no counterexample (exhaustive).
P5  At v=7 and v=9, among pairs with a^2+b^2 = 4v-2, the smallest non-zero sum over s=1..(v-1)/2 of
    (PAF_A(s)+PAF_B(s)-2)^2 equals 32.   PASS iff equal at both.
P6  4v-2 is a sum of two squares for exactly one of {487, 523, 571, 659, 739, 643}, namely 643.   PASS iff so.
P7  Gaussian bits-per-equation model b(k,v) = 0.5*log2(2*pi*k*v) - 2 reproduces the prompt's Monte-Carlo b within
    0.15 bits at all eight (k, v) points of its table, and gives b = 4.87 +- 0.15 at (k=4, v=523).   PASS iff all.

## Added after P1-P7 ran (registered before running P8)
P8  For every verified pair at v=9, PSD_A(chi) + PSD_B(chi) = 2v-2 exactly at every non-trivial character chi of Z_v
    (so the per-block PSD bound is PSD <= 2v-2, the analogue of <= 4v on the Hadamard route).  PASS iff all, to 1e-9.
