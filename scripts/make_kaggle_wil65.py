"""Generate a self-contained Kaggle notebook script for the Williamson n=65 search (stage D sharded).
Usage: make_kaggle_wil65.py > kaggle/wil65.py"""
import os
R = os.path.dirname(os.path.abspath(__file__))
src = {f: open(os.path.join(R, f)).read() for f in ("wil_common.py", "wil_join_nb.py", "wil_search.py")}
print('''# Williamson matrices of order 65: exact search, one Kaggle session = a list of shards of one row-sum case.
# Pre-registration: williamson/PREREGISTRATION.md in the repo. Every stage is exact; nothing is sampled.
#
# SETTINGS (edit these, then Run All):
CASE   = 0            # row-sum case 0..4: (9,9,7,7) (11,9,7,3) (11,11,3,3) (13,9,3,1) (15,5,3,1)
K      = 32           # shards per case (keep the same K for every session of a case)
SHARDS = [0, 1, 2, 3] # shard ids this session runs (0..K-1); NPROC run at once
NPROC  = 4
SELFTEST = True       # n=45 control first (must find 2 orbit reps / 12 quadruples)
#
# Output: /kaggle/working/wil65/wil_n65_c<CASE>_s<i>of<K>.json per shard (download them all).
# Any found quadruple is printed as a line starting with SOLUTION.
import os, sys, json, subprocess, time
D = "/kaggle/working/wil65"; os.makedirs(D, exist_ok=True)
SRC = ''' + repr(src) + '''
for f, t in SRC.items(): open(os.path.join(D, f), "w").write(t)
env = dict(os.environ, PYTHONPATH=D, NUMBA_NUM_THREADS="1")
if SELFTEST:
    out = subprocess.run([sys.executable, os.path.join(D, "wil_search.py"), "45", "5", "9", "--out", os.path.join(D, "selftest45.json")],
                         env=env, capture_output=True, text=True, cwd=D)
    r = json.load(open(os.path.join(D, "selftest45.json")))
    print("selftest n=45:", r["N_found"], r["N_expanded"], r["secs"], "s")
    assert r["N_found"] == 2 and r["N_expanded"] == 12, "SELFTEST FAILED - do not use this session"
todo = list(SHARDS); running = []
t0 = time.time()
while todo or running:
    while todo and len(running) < NPROC:
        i = todo.pop(0); o = os.path.join(D, f"wil_n65_c{CASE}_s{i}of{K}.json")
        if os.path.exists(o): print("skip existing", o); continue
        lg = open(o.replace(".json", ".log"), "w")
        running.append((i, o, subprocess.Popen([sys.executable, "-u", os.path.join(D, "wil_search.py"), "65", "5", "13",
                        "--case", str(CASE), "--shard", f"{i}/{K}", "--out", o], env=env, stdout=lg, stderr=subprocess.STDOUT, cwd=D)))
    time.sleep(30)
    for x in list(running):
        i, o, pr = x
        if pr.poll() is not None:
            running.remove(x)
            if pr.returncode or not os.path.exists(o): print(f"shard {i}: FAILED rc={pr.returncode} (see log)"); continue
            r = json.load(open(o)); c = r["cases"][0]
            print(f"shard {i}: done ledger={c['ledger']} of {c['combos_this_shard']} work={c['work']} sols={c['sols']} secs={r['secs']}", flush=True)
            for s in r["solutions"]: print("SOLUTION", json.dumps(s))
print("session done", round(time.time() - t0), "s")
''')
