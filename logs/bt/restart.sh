#!/bin/sh
# Restart every live b3-push worker after a container restart (each runner resumes from its output file; the census
# runners cache their structure enumerations next to their output, so a restart skips the enumeration).
# The kernel schedules by autogroup (one group per session): plain nice only orders processes inside one session, and
# every session gets an equal CPU share.  So each worker gets its own session (setsid), and its group's share is set
# through /proc/self/autogroup: census at 0, gates and screens at 5-15, the slow split census at 19.
cd /home/user/Hadamard
CEN=$(ls logs/bt/census/k4_P4_w*.jsonl logs/bt/census2/*.jsonl | tr '\n' ',' | sed 's/,$//')
run() {   # run PRIO LOG [ULIMIT_KB|-] CMD...
  prio=$1; log=$2; lim=$3; shift 3
  pre=""; [ "$lim" != "-" ] && pre="ulimit -v $lim; "
  setsid nohup sh -c "echo $prio > /proc/self/autogroup 2>/dev/null; ${pre}exec nice -n $prio $*" >> "$log" 2>&1 < /dev/null &
}
# census (643 with at most 4,096 endomorphisms finished at 13:05)
run 0 logs/bt/census/k4_P4_w4.log - python3 theory/bt_census_run.py logs/bt/census/k4_P4_w4.jsonl --ns 643 --rmin 8 --rmax 160 --cap 16384 --min-cap 4096 --part 0/3
run 0 logs/bt/census/k4_P4_w4.log - python3 theory/bt_census_run.py logs/bt/census/k4_P4_w4.jsonl --ns 643 --rmin 8 --rmax 160 --cap 16384 --min-cap 4096 --part 1/3
run 0 logs/bt/census/k4_P4_w4.log - python3 theory/bt_census_run.py logs/bt/census/k4_P4_w4.jsonl --ns 643 --rmin 8 --rmax 160 --cap 16384 --min-cap 4096 --part 2/3
run 0 logs/bt/census2/k8c_839.log - python3 theory/bt_census2.py logs/bt/census2/k8c_839.jsonl --fam k8 --ns 839 --rmin 12 --rmax 110 --cap 16384 --part 0/2
run 0 logs/bt/census2/k8c_839.log - python3 theory/bt_census2.py logs/bt/census2/k8c_839.jsonl --fam k8 --ns 839 --rmin 12 --rmax 110 --cap 16384 --part 1/2
run 19 logs/bt/census/split_k4.log - python3 theory/bt_census_split.py logs/bt/census/split_k4.jsonl logs/bt/census/k4_P4_w0.jsonl,logs/bt/census/k4_P4_w1.jsonl,logs/bt/census/k4_P4_w2.jsonl,logs/bt/census/k4_P4_w3.jsonl,logs/bt/census/k4_P4_w5.jsonl --rmin 16 --rmax 50
# screens [P] and exact gates [M]
run 10 logs/bt/census/thm_rigid.log - python3 theory/bt_thm_rigid.py logs/bt/census/thm_rigid.jsonl $CEN
run 10 logs/bt/closure.log - sh logs/bt/closure.sh
run 15 /dev/null - sh logs/bt/close_loop.sh
# timed exact gates on the alive structures with r <= 80 and log2 E_Q >= -5 (the window of the known successes)
run 10 logs/bt/gates_window.log 8000000 python3 theory/bt_gates_levels.py logs/bt/gates_run_window.jsonl $CEN --levels auto --keys @logs/bt/window_keys.txt --tmax 1800 --maxlen 3e7
# timed exact gates on the top 20 alive 643 hosts (levels no runner has tried yet)
run 5 logs/bt/gates_top643.log 7000000 python3 theory/bt_gates_levels.py logs/bt/gates_run_top643.jsonl $CEN --levels auto --keys @logs/bt/top643_keys.txt --tmax 900 --seen-all
