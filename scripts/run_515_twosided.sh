#!/bin/bash
# two-sided yc stage over all joint-feasible compressions of (h=85, p=5, u=$U); 2 at a time, 2 numba threads each
# usage: U=161 scripts/run_515_twosided.sh
cd /home/user/Hadamard; mkdir -p logs/twosided
U=${U:?set U}; XC=data/config_screen_h85_p5_u${U}_pass.npy
N=$(python3 -c "import numpy as np; print(len(np.load('$XC')))")
seq 0 $((N-1)) | xargs -P2 -I{} bash -c "NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 python3 -u scripts/bge3_twosided_run.py 85 5 $U 5 $XC {} 2e7 100000 > logs/twosided/u${U}_sol{}.log 2>&1"
touch logs/twosided/done_u${U}
