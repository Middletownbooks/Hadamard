### Screen totals

| values | count |
|---|---:|
| some slice alive | 286 |
| every slice dead (exact gates) | 179 |
| no slice with 24-130 orbits | 158 |
| skipped (Kumar lead / Schumer) | 7 |

| (slice, pattern) records | count |
|---|---:|
| alive | 4405 |
| dead Z_d | 612 |
| dead R | 425 |
| dead L | 257 |
| no slice in range | 158 |
| skipped: Kumar lead | 3 |
| skipped: Kumar running | 2 |
| skipped: Schumer 6644? | 1 |
| skipped: Schumer 7852? | 1 |

### Exact gates of the leading slices [M]

- **549<22>**: 4 row patterns, 4 dead; [1, 3, 31, 35]: ax0:Z3 dead (0); [1, 7, 25, 39]: ax0:Z3 dead (0) ...
- **549<16>**: 4 row patterns, 4 dead; [1, 3, 31, 35]: ax0:Z3 dead (0); [1, 7, 25, 39]: ax0:Z3 dead (0) ...
- **361torus(4,7)**: 4 row patterns, 0 dead; [1, 5, 7, 37]: ax0:Z19 pass (43), ax1:Z19 pass (2338173); [1, 7, 13, 35]: ax0:Z19 pass (45), ax1:Z19 pass (2207583) ...
- **325singer3x3**: 5 row patterns, 0 dead; [1, 5, 7, 35]: ax0:Z5 n/a, ax1:Z5 n/a, ax2:Z13 pass (44952); [1, 7, 17, 31]: ax0:Z5 n/a, ax1:Z5 n/a, ax2:Z13 pass (46473) ...
- **549 unip [('unip', [[1, 1], [0, 1]], 9)]**: 9 row patterns, 3 dead; [1, 1, 13, 45]: ax0:Z3 n/a, ax1:Z3 dead (0), ax2:Z61 undecided; [1, 1, 13, 45]: ax1:Z3 dead (0), Z3^2 dead (0) ...
- **549 unip [('unip', [[1, 1], [0, 1]], 12)]**: 4 row patterns, 4 dead; [1, 3, 31, 35]: ax1:Z3 dead (0), Z3^2 dead (0); [1, 7, 25, 39]: ax1:Z3 dead (0), Z3^2 dead (0) ...
- **549 unip [('unip', [[1, 1], [0, 1]], 1), ('id', None, 12)]**: 4 row patterns, 4 dead; [1, 3, 31, 35]: ax1:Z3 dead (0), Z3^2 dead (0); [1, 7, 25, 39]: ax1:Z3 dead (0), Z3^2 dead (0) ...
- **549 unip [('unip', [[1, 1], [0, 1]], 1), ('id', None, 13)]**: 8 row patterns, 0 dead; [1, 1, 13, 45]: ax1:Z3 pass (48), ax0:Z3 n/a, ax2:Z61 undecided, Z3^2 pass (16797); [1, 3, 31, 35]: ax1:Z3 pass (38), ax0:Z3 n/a, ax2:Z61 undecided, Z3^2 pass (20763) ...

### Ranked screen table (top 30 by k-tree cost) [E]

Cost = log10 of pair work per solution in a four-list k-tree with lists (min(N_A, 2^24), 5e7, 3e7, 3e7), from the analytic in-list estimate (skew_lib.inlist_estimate). On this scale the solved controls are 613: 10.2, 1123: 11.4, 761: 12.9.

| # | v | slice | \|H\| | orbits | classes | rows | log10 E (slice) | log10 E in lists | cost | gates | next step |
|---:|---:|---|---:|---:|---:|---|---:|---:|---:|---|---|
| 1 | 549 | G=Z549, H=[[22]] | 15 | 41 | 20 | [1, 3, 31, 35] | 8.31 | 2.76 | 12.22 | dead at Z3 [M] (all patterns) | none (dead) |
| 2 | 361 | G=Z19xZ19, H=[('torus', [[4, 0], [0, 7]])] | 9 | 45 | 22 | [1, 5, 7, 37] | 11.86 | 2.8 | 12.25 | alive: both Z19 compressions pass [M] | searched (§5) |
| 3 | 325 | G=Z5xZ5xZ13, H=[('singer', [[1, 1], [2, 3]], 1), ('id', None | 9 | 45 | 22 | [1, 7, 17, 31] | 11.1 | 2.72 | 12.33 | alive: Z13 passes; Z5^2 undecided (cap) | searched (§5) |
| 4 | 433 | G=Z433, H=[[27]] | 9 | 49 | 24 | [1, 1, 19, 37] | 14.91 | 0.83 | 14.41 | prime: no Z_d gate | k-tree beyond budget (cost >= 14) |
| 5 | 1483 | G=Z1483, H=[[28]] | 39 | 39 | 19 | [1, 1, 1, 77] | 6.97 | 0.38 | 14.58 | prime: no Z_d gate | k-tree beyond budget (cost >= 14) |
| 6 | 589 | G=Z589, H=[[102]] | 15 | 45 | 22 | [1, 5, 31, 37] | 8.53 | 0.17 | 14.88 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 7 | 549 | G=Z3xZ3xZ61, H=[('torus', [[1, 0], [0, 1]], 12)] | 15 | 45 | 22 | [1, 3, 31, 35] | 5.98 | -1.5 | 16.55 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 8 | 487 | G=Z487, H=[[41]] | 9 | 55 | 27 | [1, 17, 17, 37] | 16.17 | -1.44 | 16.68 | prime: no Z_d gate | k-tree beyond budget (cost >= 14) |
| 9 | 817 | G=Z817, H=[[68]] | 21 | 45 | 22 | [1, 7, 37, 43] | 6.54 | -1.7 | 16.75 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 10 | 725 | G=Z5xZ5xZ29, H=[('singer', [[1, 1], [2, 3]], 7)] | 21 | 45 | 22 | [1, 21, 33, 37] | 4.46 | -1.79 | 16.84 | gates pass | k-tree beyond budget (cost >= 14) |
| 11 | 589 | G=Z589, H=[[64]] | 15 | 49 | 24 | [1, 5, 11, 47] | 8.17 | -1.95 | 17.19 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 12 | 489 | G=Z489, H=[[40]] | 9 | 57 | 28 | [1, 15, 19, 37] | 15.66 | -2.09 | 17.33 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 13 | 481 | G=Z481, H=[[9]] | 9 | 57 | 28 | [1, 1, 31, 31] | 14.83 | -2.35 | 17.59 | gates pass | k-tree beyond budget (cost >= 14) |
| 14 | 723 | G=Z723, H=[[94]] | 15 | 51 | 25 | [1, 29, 31, 33] | 11.82 | -2.45 | 17.69 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 15 | 481 | G=Z481, H=[[16]] | 9 | 57 | 28 | [1, 5, 23, 37] | 14.77 | -2.47 | 17.71 | gates pass | k-tree beyond budget (cost >= 14) |
| 16 | 949 | G=Z949, H=[[16], [74]] | 27 | 45 | 22 | [1, 19, 25, 53] | 3.89 | -3.04 | 18.09 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 17 | 523 | G=Z523, H=[[19]] | 9 | 59 | 29 | [1, 19, 19, 37] | 16.95 | -3.24 | 18.48 | prime: no Z_d gate | k-tree beyond budget (cost >= 14) |
| 18 | 597 | G=Z597, H=[[61]] | 11 | 57 | 28 | [1, 1, 19, 45] | 14.53 | -3.3 | 18.54 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 19 | 817 | G=Z817, H=[[11]] | 21 | 49 | 24 | [1, 7, 37, 43] | 6.15 | -3.65 | 18.89 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 20 | 655 | G=Z655, H=[[176]] | 13 | 55 | 27 | [1, 3, 3, 51] | 11.68 | -4.09 | 19.33 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 21 | 1565 | G=Z1565, H=[[16]] | 39 | 45 | 22 | [1, 3, 3, 79] | 4.58 | -4.41 | 19.45 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 22 | 1201 | G=Z1201, H=[[96]] | 25 | 49 | 24 | [1, 1, 49, 49] | 9.86 | -4.35 | 19.59 | prime: no Z_d gate | k-tree beyond budget (cost >= 14) |
| 23 | 1203 | G=Z1203, H=[[25]] | 25 | 51 | 25 | [1, 1, 47, 51] | 9.15 | -5.09 | 20.33 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 24 | 793 | G=Z793, H=[[22]] | 15 | 57 | 28 | [1, 29, 31, 37] | 11.43 | -5.58 | 20.82 | gates pass | k-tree beyond budget (cost >= 14) |
| 25 | 921 | G=Z921, H=[[64]] | 17 | 57 | 28 | [1, 33, 35, 37] | 11.93 | -5.82 | 21.06 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 26 | 361 | G=Z19xZ19, H=[('torus', [[1, 0], [0, 4]])] | 9 | 57 | 28 | [1, 5, 7, 37] | 4.56 | -6.3 | 21.54 | gates pass | k-tree beyond budget (cost >= 14) |
| 27 | 1301 | G=Z1301, H=[[44]] | 25 | 53 | 26 | [1, 1, 51, 51] | 10.25 | -6.33 | 21.57 | prime: no Z_d gate | k-tree beyond budget (cost >= 14) |
| 28 | 1103 | G=Z1103, H=[[17]] | 19 | 59 | 29 | [1, 37, 39, 39] | 12.72 | -7.37 | 22.61 | prime: no Z_d gate | k-tree beyond budget (cost >= 14) |
| 29 | 899 | G=Z899, H=[[25]] | 21 | 55 | 27 | [1, 15, 39, 43] | 4.58 | -7.59 | 22.83 | gates undecided/none | k-tree beyond budget (cost >= 14) |
| 30 | 1255 | G=Z1255, H=[[16]] | 25 | 55 | 27 | [1, 1, 47, 53] | 7.6 | -8.06 | 23.3 | gates undecided/none | k-tree beyond budget (cost >= 14) |
