r"""Gate A2: the calibration table of the q = 13 plug-in analogues (K = Z_n, n = |K|), assembled from the logs [M, E].
Columns: n, v = 13 n + 1, defect allocations (theory/h2396_plugin_system.py identity), the repository estimator
(logs/h2396/gauss_analogues.jsonl: stage-1 log2E, E per y for a synthetic y, total), the exact stage-1 count
(logs/h2396/s1enum, stage1_tool), the exact number of plug-in cores (stage 2 on every stage-1 class), and the
sampled instances (CP-SAT level 1, exact level 2; logs/h2396/s1cp, sampler_n*) with the estimator at the real y.
usage: python3 theory/h2396_a2_table.py [--md]"""
import sys, os, json, glob, math, re
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
L = lambda *p: os.path.join(root, 'logs', 'h2396', *p)

def allocs(n, q=13):
    R = 4 * (n - 1)
    return [(a, (R - (q + 1) * a) // (q - 3)) for a in range(R // (q + 1) + 1) if (R - (q + 1) * a) % (q - 3) == 0]

def gauss():
    out = {}
    for l in open(L('gauss_analogues.jsonl')):
        d = json.loads(l); key = json.loads(d['key']); n = (key[1] - 1) // 13
        s1 = d['stage1'].get('stage1_log2E'); ey = d.get('E_per_y')
        out[n] = dict(s1=s1, Ey=ey, tot=None if s1 is None or ey is None else round(s1 + ey, 1))
    out[46] = dict(s1=114.69, Ey=-7.18, tot=107.5, note='energy-only stage 1 (logs/open_screen/deficit_s*.jsonl)')
    return out

def exact():
    ex = {6: dict(s1=0, cores=0, how='stage 1 empty: repo tool A+B, h2396_s1enum, h2396_s1cp'),
          8: dict(s1=0, cores=0, how='stage 1 empty: repo tool A+B, h2396_s1enum, h2396_s1cp')}
    for f in glob.glob(L('s1enum', 'summary.log')) + glob.glob(L('s1enum', 'run_n*.log')):
        for line in open(f):
            m = re.search(r'n=(\d+) .*stage1_ordered_count=(\d+)', line)
            if m: ex.setdefault(int(m.group(1)), {})['s1'] = int(m.group(2))
    for f in glob.glob(L('s1enum', 'stage2_n*.log')):
        n = int(re.search(r'stage2_n(\d+)', f).group(1))
        for line in open(f):
            if line.startswith('WEIGHTED_TOTAL'): ex.setdefault(n, {})['cores'] = int(line.split()[1])
    # --zqpre runs: a zero count of stage-1 solutions in the Z_q-alive profiles means no core [P + M]
    for f in glob.glob(L('s1enum', 'run_n*_zq*.log')):
        for line in open(f):
            m = re.search(r'n=(\d+) .*stage1_ordered_count_in_Zq_alive_profiles=(\d+)', line)
            if m and int(m.group(2)) == 0: ex.setdefault(int(m.group(1)), {})['cores'] = 0
    return ex

def sampled():
    """distinct sampled stage-1 instances per n with their exact completion counts"""
    out = {}; seen = set()
    first = L('s1cp', 'n18_y0.txt')                     # the first CP-SAT instance: COUNT 0 (h2396_stage2)
    if os.path.exists(first):
        seen.add(open(first).read().strip()); out.setdefault(18, []).append(0)
    for f in glob.glob(L('sampler_n*', 'inst.jsonl')):
        n = int(re.search(r'sampler_n(\d+)', f).group(1))
        for line in open(f):
            d = json.loads(line)
            if 'stage2' in d and d['stage2']:
                txt = open(os.path.join(os.path.dirname(f), f"inst_{d['key']}.txt")).read().strip()
                if txt in seen: continue
                seen.add(txt); out.setdefault(n, []).append(int(d['stage2'].split()[1]))
    return out

if __name__ == "__main__":
    G = gauss(); X = exact(); S = sampled()
    rows = []
    for n in sorted(set(G) | set(X) | set(S)):
        g = G.get(n, {}); x = X.get(n, {}); s = S.get(n, [])
        rows.append(dict(n=n, v=13 * n + 1, allocations=allocs(n), gauss_stage1_log2E=g.get('s1'), gauss_E_per_y=g.get('Ey'),
                         gauss_total_log2E=g.get('tot'), exact_stage1_ordered=x.get('s1'), exact_cores=x.get('cores'),
                         sampled_instances=len(s), sampled_completions=sum(c for c in s if c is not None) if s else None,
                         sampled_list=s))
    if '--md' in sys.argv:
        print('| n = \\|K\\| | v | allocations (e+, e-) | Gaussian stage 1 (log2) | Gaussian E per y | Gaussian total (log2) | exact stage 1 (ordered) | exact cores | sampled y: completions each |')
        print('|---|---|---|---|---|---|---|---|---|')
        for r in rows:
            f = lambda v: '–' if v is None else v
            st1 = r['exact_stage1_ordered']; st1 = '–' if st1 is None else (f"{st1} (2^{math.log2(st1):.1f})" if st1 else '0')
            smp = '–' if not r['sampled_instances'] else f"{r['sampled_instances']} y: " + ", ".join(str(c) for c in r['sampled_list'])
            print(f"| {r['n']} | {r['v']} | {r['allocations'] or 'none'} | {f(r['gauss_stage1_log2E'])} | {f(r['gauss_E_per_y'])} | {f(r['gauss_total_log2E'])} | {st1} | {f(r['exact_cores'])} | {smp} |")
    else:
        for r in rows: print(json.dumps(r))
