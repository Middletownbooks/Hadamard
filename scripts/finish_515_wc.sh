#!/bin/bash
# unattended: when the fast-lift queue finishes, record results and commit/push them
cd /home/user/Hadamard
while [ ! -f logs/wc/done ]; do sleep 300; done
python3 scripts/record_wc.py > logs/wc/final_summary.txt 2>&1
git add data/wc_lift_515_u161.jsonl data/solution_v515_wc_*.txt 2>/dev/null
git commit -qm "515 u=161: fast exact lift queue finished (all yc classes of sols 3,4,5,12,13,14); $(cat logs/wc/final_summary.txt)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01Gh3Zhrg9pjR7nGAEcRMoMR"
for i in 1 2 3 4; do git push -q -u origin claude/focused-planck-j1645f && break; sleep $((2**i)); done
