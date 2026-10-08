"""Driver for generalized bordered-GS options: cost every pre-checked option, screen the cheapest, rank.

Input: options JSONL from bge3_options.py.  For each option with precheck=True (and c >= 4):
  1. bge3_cost_psi.py (Gaussian cost, symmetry groups preserving the target), smallest r first, STOPAFTER positives;
  2. for the TOP cheapest configurations per v: bge3_normscreen_psi.py at every proper divisor t of h,
     then bge3_screen_psi.py at t = 1 (exact DP, can overflow -> 'undecided');
  3. JSONL rows: option, best configuration, n, log2E, list-method bound (n - log2E)/2, screen verdicts.
usage: bge3_general_run.py options.jsonl out.jsonl [--par 4] [--top 3] [--rmax 160] [--stopafter 3] [--only v,v]"""
import sys, json, os, subprocess, itertools
from concurrent.futures import ThreadPoolExecutor

a = sys.argv[1:]; opts = {"--par": "4", "--top": "3", "--rmax": "160", "--stopafter": "3", "--only": "", "--maxcfg": "40", "--resume": ""}
for k in opts:
    if k in a: opts[k] = a[a.index(k) + 1]
inp, outp = a[0], a[1]
only = {int(x) for x in opts["--only"].split(",") if x}

def psi_arg(o):
    if o["b"] == o["c"] - 1: return "all"
    return ";".join(",".join(str(x) for x in k) for k in o["psi"])

def run(cmd, timeout, env=None):
    r = subprocess.run(["timeout", str(timeout)] + cmd, capture_output=True, text=True, env=env)
    return [json.loads(l) for l in r.stdout.splitlines() if l.startswith("{")], r.stderr[-400:]

def cost(o):
    env = dict(os.environ, STOPAFTER=opts["--stopafter"], MAXCFG=opts.get("--maxcfg", "40"), OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    recs, err = run(["python3", "scripts/bge3_cost_psi.py", str(o["v"]), ",".join(map(str, o["Q"])), psi_arg(o),
                     "20", opts["--rmax"]], 7200, env)
    rows = [r for r in recs if "log2E" in r and r["log2E"] > 0]
    return o, rows, (None if recs else err)

options = [json.loads(l) for l in open(inp)]
options = [o for o in options if o["precheck"] and o["c"] >= 4 and (not only or o["v"] in only)]
done_keys = set()
for f in [x for x in opts["--resume"].split(",") if x]:
    for l in open(f):
        try: r = json.loads(l)
        except Exception: continue
        if r.get("stage") == "cost": done_keys.add((r["v"], tuple(r["Q"]), r["psi"]))
options = [o for o in options if (o["v"], tuple(o["Q"]), psi_arg(o)) not in done_keys]
print(json.dumps({"stage": "start", "options": len(options)}), flush=True)
results = []
with ThreadPoolExecutor(int(opts["--par"])) as ex:
    for o, rows, err in ex.map(cost, options):
        best = min(rows, key=lambda r: r["n"]) if rows else None
        rec = dict(stage="cost", v=o["v"], c=o["c"], Q=o["Q"], b=o["b"], h=o["h"], psi=psi_arg(o), err=err,
                   n_positive=len(rows), best=best,
                   bound=round((best["n"] - best["log2E"]) / 2, 1) if best else None)
        results.append((o, rows)); print(json.dumps(rec), flush=True)

# screen the cheapest configurations per v
byv = {}
for o, rows in results:
    for r in rows: byv.setdefault(o["v"], []).append((r["n"], o, r))
def screen(item):
    n, o, r = item
    h = o["h"]; ts = [str(t) for t in range(1, h) if h % t == 0]
    base = [str(o["v"]), ",".join(map(str, o["Q"])), psi_arg(o), str(r["d"]), str(r["r"]), str(r["L"])]
    norm, e1 = run(["python3", "scripts/bge3_normscreen_psi.py"] + base + ts, 3600)
    norm_ok = any(all(g["tests"].values()) for g in norm) if norm else None
    dp, e2 = (run(["python3", "scripts/bge3_screen_psi.py"] + base + ["1"], 3600, dict(os.environ, MAXST="1e7"))
              if norm_ok else ([], None))
    dpv = [g["tests"].get("1", g["tests"].get(1)) for g in dp if not g.get("partial")]
    verdict = ("dead (norm)" if norm_ok is False else
               "dead (DP t=1)" if dpv and all(x is False for x in dpv) else
               "alive" if any(x is True for x in dpv) else "undecided")
    return dict(stage="screen", v=o["v"], c=o["c"], Q=o["Q"], b=o["b"], h=o["h"], psi=psi_arg(o), d=r["d"], r=r["r"],
                L=r["L"], n=r["n"], log2E=r["log2E"], bound=round((r["n"] - r["log2E"]) / 2, 1),
                norm=[g["tests"] for g in norm], dp=dpv, verdict=verdict)
items = [it for v in byv for it in sorted(byv[v], key=lambda x: x[0])[: int(opts["--top"])]]
with open(outp, "w") as f, ThreadPoolExecutor(int(opts["--par"])) as ex:
    for res in ex.map(screen, items):
        print(json.dumps(res), flush=True); f.write(json.dumps(res) + "\n")
