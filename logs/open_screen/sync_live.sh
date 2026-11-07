#!/bin/bash
W=${W:-/tmp/claude-0/-home-user-Hadamard/23d0cbad-e05e-56bb-b7b0-b45b94a4a068/scratchpad/live}
cd "$(dirname "$0")/../.."
cp $W/deficit_s*.jsonl $W/b0_prime_direct.jsonl logs/open_screen/
[ -f $W/z2gate.jsonl ] && cp $W/z2gate.jsonl logs/open_screen/
[ -f $W/stage1_exact.jsonl ] && cp $W/stage1_exact.jsonl logs/open_screen/stage1_exact_step5.jsonl
true
