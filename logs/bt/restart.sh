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
run 10 logs/bt/closure.log - sh logs/bt/closure.sh
run 5 logs/bt/census2/close_k8c.log 6000000 python3 theory/bt_close.py logs/bt/census2/close.jsonl logs/bt/census2/k8c_839.jsonl
# (the window gates, logs/bt/gates_run_window.jsonl, finished at 15:21)
# timed exact gates on the top 20 alive 643 hosts (levels no runner has tried yet)
run 5 logs/bt/gates_top643.log 7000000 python3 theory/bt_gates_levels.py logs/bt/gates_run_top643.jsonl $CEN --levels auto --keys @logs/bt/top643_keys.txt --tmax 900 --seen-all
# census extension: 739 on Z2^3 x Z4 (2^17 endomorphisms) and 643 on Z2 x Z4 x Z16 (2^15); the enumeration is cached once done
run 0 logs/bt/census/k4_P4_w10.log - python3 theory/bt_census_run.py logs/bt/census/k4_P4_w10.jsonl --ns 643 --rmin 8 --rmax 160 --cap 32768 --min-cap 16384
# k4one at the three orders the first k4one queue left out (12 <= r <= 90, as at the other orders)
run 5 logs/bt/census2/close_k4one739.log 6000000 python3 theory/bt_close.py logs/bt/census2/close.jsonl logs/bt/census2/k4one_739.jsonl

