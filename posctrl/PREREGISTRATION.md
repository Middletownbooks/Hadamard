# Blind positive control for the c=6 compression-and-lift pipeline: registered predictions (before any run)

Rung v = 95 (h = 15 = 3*5, G = Z2 x Z3 x Z15), 347 corner/border template, blocks invariant under a multiplier u.
Slices: p in {3, 5} x u in {7, 17, 11, 19, 29} (orders 4 and 2).  Per slice: exact stage A (every shift, every
compressed block, exact 4-list join -> ALL compressions), then exact lift of every compression (full lift lists,
exact join), every hit assembled into the 380 x 380 matrix and checked H H^T = 380 I in integer arithmetic.
Single thresholds; an unevaluable bar counts as FAIL.

PC1  At least one verified Hadamard matrix of order 380 is found over the 10 slices.        PASS iff >= 1.
PC2  The forced-residue filters (per-class residue_ok, then pairwise joint at M=16 and 64) never reject a
     compression that has >= 1 exact lift.                                                  PASS iff 0 false rejections.
PC3  For every liftable compression and every yc class occurring among its exact solutions, the fast three-sided
     lift (wc_lift) returns exactly the same set of solutions in that class as the exact join. PASS iff all identical.
PC4  Negative test: the object check rejects a solution with one entry flipped.             PASS iff rejected.
Reported, no threshold: per slice #compressions, #liftable, lift rate, #solutions.

## Addendum (registered before running): larger rungs with high-order multipliers
v=95 is only exhaustively searchable for order-4 multipliers: two p=3 slices are provably empty (0 compressions);
p=5,u=17 has 837,168 compressions (being lifted in full); the other slices are too large (stage A or lift).
Added slices, all exhaustive:
- v=131 (h=21): p in {3, 7} x u in {5, 11, 19} (order 6)
- v=215 (h=35): p in {5, 7} x u in {17, 23, 37, 73} (order 12), and p=7 x u in {11, 19, 59, 79, 101, 149} (order 6)
PC1b  At least one verified Hadamard matrix (order 524 or 860) is found over these slices.   PASS iff >= 1.
PC2, PC3 and PC4 apply unchanged to every slice.
