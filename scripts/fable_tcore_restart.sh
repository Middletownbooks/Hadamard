#!/bin/bash
cd /home/user/Hadamard
while ps -eo args | grep -E "^python3 theory/fable_tcore.py" | grep -v grep > /dev/null; do sleep 20; done
exec scripts/fable_tcore_queue.sh
