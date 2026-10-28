#!/bin/sh
# Restart every live b3-push worker after a container restart (each runner resumes from its output file; the census
# runners cache their structure enumerations next to their output, so a restart skips the enumeration).
# The kernel schedules by autogroup (one group per session): plain nice only orders processes inside one session, and
# every session gets an equal CPU share.  So each worker gets its own session (setsid), and its group's share is set
# through /proc/self/autogroup: census at 0, gates and screens at 5-15, the slow split census at 19.
cd /home/user/Hadamard
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1      # the matrices are small; idle BLAS threads only add load
CEN=$(ls logs/bt/census/k4_P4_w*.jsonl logs/bt/census2/*.jsonl | tr '\n' ',' | sed 's/,$//')
run() {   # run PRIO LOG [ULIMIT_KB|-] CMD...
  prio=$1; log=$2; lim=$3; shift 3
  pre=""; [ "$lim" != "-" ] && pre="ulimit -v $lim; "
  setsid nohup sh -c "echo $prio > /proc/self/autogroup 2>/dev/null; ${pre}exec nice -n $prio $*" >> "$log" 2>&1 < /dev/null &
}
# census (643 with at most 4,096 endomorphisms finished at 13:05; 839 k8 finished at 16:15; 643 up to 16,384 endomorphisms at 16:35)
run 19 logs/bt/census/split_k4.log - python3 theory/bt_census_split.py logs/bt/census/split_k4.jsonl logs/bt/census/k4_P4_w0.jsonl,logs/bt/census/k4_P4_w1.jsonl,logs/bt/census/k4_P4_w2.jsonl,logs/bt/census/k4_P4_w3.jsonl,logs/bt/census/k4_P4_w5.jsonl --rmin 16 --rmax 50
# screens [P] and exact gates [M]
run 10 logs/bt/census/thm_rigid.log - python3 theory/bt_thm_rigid.py logs/bt/census/thm_rigid.jsonl $CEN
run 10 logs/bt/closure.log - sh logs/bt/closure.sh            # tiers 1-3, 10 s then 60 s per level, 3 runners
# (the window gates, logs/bt/gates_run_window.jsonl, finished at 15:21)
# timed exact gates on the top 20 alive 643 hosts (levels no runner has tried yet)
# census extension: 739 on Z2^3 x Z4 (2^17 endomorphisms) and 643 on Z2 x Z4 x Z16 (2^15); the enumeration is cached once done
# k4one at the three orders the first k4one queue left out (12 <= r <= 90, as at the other orders)
# (closure2, tier 2 from the most expensive end, was stopped on Oct 27: 46 structures, 0 kills; closure.sh covers tier 2)
