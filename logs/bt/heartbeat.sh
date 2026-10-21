#!/bin/bash
# Heartbeat for the b3-push session (run under the Monitor tool): one status line every INT seconds (default 1500),
# plus an immediate line when a theory/bt_* worker ends (a job finished or died).  Every line wakes the session and
# costs usage, so the interval is long and worker starts are not reported.
# usage: bash logs/bt/heartbeat.sh [INT]
INT=${1:-1500}
cd /home/user/Hadamard || exit 1
workers() {
  # the bracket keeps this script's own command line from matching
  ps -eo args | grep -E "^python3 theor[y]/bt_" | sed -E 's#^python3 theory/(bt_[a-z0-9_]+)\.py ([^ ]*).*#\1:\2#; s#logs/bt/##' | sort
}
prev=$(workers)
last=0
while true; do
  sleep 20
  cur=$(workers)
  now=$(date +%s)
  if [ "$cur" != "$prev" ]; then
    gone=$(comm -23 <(echo "$prev") <(echo "$cur") | tr '\n' ' ')
    [ -n "$gone" ] && echo "ENDED $(date -u +%H:%M) ${gone} | running: $(echo "$cur" | grep -c .)"
    prev=$cur
  fi
  if [ $((now - last)) -ge "$INT" ]; then
    up=$(awk '{print int($1/60)}' /proc/uptime)
    load=$(cut -d' ' -f1 /proc/loadavg)
    lines=$(cat logs/bt/census/k4_P4_w*.jsonl logs/bt/census2/*.jsonl 2>/dev/null | wc -l)
    dirty=$(git status --short logs/bt 2>/dev/null | wc -l)
    echo "HB $(date -u +%H:%M) up=${up}m workers=$(echo "$cur" | grep -c .) load=${load} census_lines=${lines} uncommitted=${dirty}"
    last=$now
  fi
done
