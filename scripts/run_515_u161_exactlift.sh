#!/bin/bash
# exact two-sided lift on all new yc classes of u=161 sols 3 (new), 4, 5 (split in halves); 2 jobs at a time
cd /home/user/Hadamard; D=data/twosided_h85_p5_u161_q5
job(){ # name sol classes_file
  JOIN_BATCH=2e8 NUMBA_NUM_THREADS=2 CLASSES_FILE=$3 RESUME_LOG=logs/twosided/lift_$1.log.prev \
    python3 -u scripts/bge3_twosided_lift.py 85 5 161 5 data/config_screen_h85_p5_u161_pass.npy $2 >> logs/twosided/lift_$1.log 2>&1; }
export -f job
printf "s3new 3 ${D}_sol3_classes_new.npy\ns4 4 ${D}_sol4_classes.npy\ns5a 5 ${D}_sol5_classes_a.npy\ns5b 5 ${D}_sol5_classes_b.npy\n" | \
  xargs -P2 -L1 bash -c 'job "$0" "$1" "$2"'
touch logs/twosided/done_u161_exactlift
