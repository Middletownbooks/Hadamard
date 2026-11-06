#!/bin/bash
W=${W:-/tmp/claude-0/-home-user-Hadamard/23d0cbad-e05e-56bb-b7b0-b45b94a4a068/scratchpad/live}
cd "$(dirname "$0")/../.."
cp $W/deficit_s*.jsonl $W/b0_prime_direct.jsonl logs/open_screen/
