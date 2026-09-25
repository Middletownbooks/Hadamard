#!/bin/bash
# sampled stage B on the 12 joint-residue-feasible compressions of (h=85, p=5, u=59); 2 at a time
cd /home/user/Hadamard; mkdir -p logs/stageB_u59
seq 0 11 | xargs -P2 -I{} bash -c "XCFILE=data/config_screen_h85_p5_u59_pass.npy python3 -u scripts/bge3_stageB_sampled.py 85 5 59 {} 4000000 400 {} > logs/stageB_u59/sol{}.log 2>&1"
touch logs/stageB_u59/done
