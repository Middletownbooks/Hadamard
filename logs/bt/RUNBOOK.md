# b3 push: runbook for the monitoring session

Written 2026-10-28 for a session that keeps the b3 push running without changing it. The task is
`report/prompt_b3_push.md`; the write-up is `report/b3_push_answer.md`. Your job: keep the workers alive, commit their
logs, and refresh the generated report tables. Code changes, new computations and new claims in the report belong to a
high-effort session; when in doubt, stop and report.

## What is running

| Worker | Process | Writes |
|---|---|---|
| Exact closure pass | `sh logs/bt/closure.sh`, which starts up to three `python3 theory/bt_gates_levels.py` runners | `logs/bt/gates_run_closure.jsonl`, `logs/bt/closure.log` |
| Split census | `python3 theory/bt_census_split.py` (lowest priority, slow) | `logs/bt/census/split_k4.jsonl`, `.log` |

Both resume from their output files after any restart. `logs/bt/restart.sh` starts them; nothing else needs starting.

## Normal and abnormal signs

Normal:
- Heartbeat `HB ... workers=4`: three closure runners and the split census. A restarted heartbeat prints its first
  `HB` line after 15 minutes, not at once.
- `workers=2` while a single-runner retry runs (memory retry or large-list retry, 11 GiB). This can last a few hours.
- `ENDED ... bt_gates_levels` events: a pass ended, or a runner hit its memory cap and was restarted (exit code 75; a
  new PID appears within seconds).
- New lines in `logs/bt/closure.log` such as `closure round 1 tier 2 tmax 60 done`.

Abnormal: report these to the user and change nothing:
- A new `Traceback` in `logs/bt/closure.log`. The newest known one is at line 7865 (Oct 27, already fixed). Check with
  `grep -n Traceback logs/bt/closure.log | tail -1`.
- An out-of-memory kill: `dmesg | grep 'Killed process'` prints a line.
- `workers=1` or `workers=0` on two heartbeats in a row while `closure.log` does not end with `closure pass finished`.
- `bt_summary.py --fill` fails, or `git diff report/` changes text outside the generated blocks (`<!-- ... -->` markers).
- A merge conflict or a rejected push.

## Every hour (the scheduled check-in)

1. `uptime`. If the machine is up less than 10 minutes, the container restarted: do "After a container restart" below.
2. If the heartbeat Monitor is not running, start it again: command `bash /home/user/Hadamard/logs/bt/heartbeat.sh 900`,
   timeout 1800000. It expires every 30 minutes; restart it each time it expires.
3. Refresh the tables: `cd /home/user/Hadamard && timeout 900 python3 ./theory/bt_summary.py --fill`.
4. Commit and push (block below).
5. Tell the user in two or three lines: workers running, closure stage
   (`grep 'closure round' logs/bt/closure.log | tail -1`), and the `Over N structures: ...` line of the report
   (`grep '^Over' report/b3_push_answer.md`), which holds the deeper-level kill count and the alive count.

Step 4 of the check-in prompt ("cost and gates on new survivors") means only this table refresh. Do not start gates.

## After a container restart

1. Check that every JSONL line parses (a restart can cut a line in half):

   ```
   cd /home/user/Hadamard && python3 - <<'EOF'
   import json, glob
   bad = []
   for f in sorted(glob.glob('logs/bt/**/*.jsonl', recursive=True)):
       for i, ln in enumerate(open(f, 'rb').read().split(b'\n')):
           if ln.strip():
               try: json.loads(ln)
               except Exception: bad.append((f, i + 1))
   print('bad lines:', bad)
   EOF
   ```

   If it prints any bad line, stop and report it. Do not edit or delete log lines.
2. `cd /home/user/Hadamard && sh logs/bt/restart.sh`
3. Start the heartbeat Monitor again (as above).
4. Commit and push.

## Commit and push

```
cd /home/user/Hadamard && git add -A logs/bt report/b3_push_answer.md && git commit -q -m "b3 push: worker logs

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01HJEKq2rZaG1TEC71abxas2" && git fetch -q origin claude/focused-planck-j1645f && git merge -q --no-edit origin/claude/focused-planck-j1645f && git push -u origin HEAD:claude/focused-planck-j1645f
```

The stop hook asks for this after most turns; run the same block. Push only to `claude/focused-planck-j1645f`: never to
`main`, no pull requests, no force pushes. Keep the two trailer lines exactly as above. A session other than this one
replaces the `Claude-Session` line with its own link.

## Never

- Kill a process this session did not start. Use PIDs only, never `pkill -f`. The workers above count as started by
  the session that ran `restart.sh`.
- Edit `logs/bt/closure.sh`, `logs/bt/restart.sh`, `logs/bt/closure_keys.py` or anything in `theory/`. A running
  `sh` script must not be edited, and the settings in it were chosen deliberately.
- Start other computations, or run heavy Python next to the workers. All processes share one 13.4 GiB memory group.
  `bt_summary.py --fill` peaks at 1.3 GB and is safe; nothing larger is.
- Edit or delete lines in any `*.jsonl` log.
- Write new claims into the report. Only `bt_summary.py --fill` changes it. If a family's alive count reaches 0, the
  generated block in §5.1 says "closed"; tell the user and leave the handwritten lists to a high-effort session.

## When the closure pass finishes

`logs/bt/closure.log` ends with `closure pass finished (round N)`. Refresh the tables, commit and push, and tell the user
that the pass is done and the final report pass (task 12) is due. Leave the split census running.

## Expected timeline (from Oct 28, 04:40 UTC, without long outages)

| Stage | Estimate |
|---|---|
| Tier 1 large-list retry | under 1 h |
| Tier 2 (k4one at 643, 955, 971): 60 s retry, memory retry, large-list retry | about 6–9 h |
| Tier 3 (6,149 structures): 10 s pass | about 34 h (±2×) |
| Tier 3: 60 s retry and memory retry | likely under a day |
| Round 2 (re-sweep; ends when it adds no record) | a few hours |

About 2–3 days of machine time in all.

## Where the results stand (Oct 28)

- (A) No Hadamard matrix; more running time will not change this (report §0 and §8).
- (B) The census is complete: 170,061 structures, 9,654 alive at 04:28.
- (C) Families closed exactly include k4one at 739, which the closure pass closed on Oct 28 (report §0 and §5.1).
