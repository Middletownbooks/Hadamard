r"""Aggregated exact compression screens (theory/s1z_zqscreen.py)
over the configurations of logs/s1z/ts/configs.json, most interesting first (stage 1 above 2^-4 per symmetry class by
the z-screen, then by stage-2 cost); restart-safe per multiset.  Runs in-process (killing it stops everything).
usage: python3 theory/s1z_agg_batch.py [--part i/W] [--memgb 4.5]
(all levels of the configuration, the full Q x Z_q level last; a level whose group product exceeds the cap leaves the
multiset undecided)"""
import sys, os, json, resource
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_s1sd as SD
from s1z_zqscreen import run_config
from s1z_lib import log2_symmetry

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
a = sys.argv[1:]
resource.setrlimit(resource.RLIMIT_AS, (int(float(a[a.index('--memgb') + 1]) if '--memgb' in a else 4.5) * 2 ** 30,) * 2)
part, nW = (int(x) for x in (a[a.index('--part') + 1] if '--part' in a else '0/1').split('/'))
cfgs = json.load(open(os.path.join(ROOT, 'logs/s1z/ts/configs.json')))
zs = {}
p = os.path.join(ROOT, 'logs/s1z/zscreen_skipped.jsonl')
if os.path.exists(p):
    for l in open(p):
        r = json.loads(l); zs[json.dumps(r['cfgkey'])] = r
skip = {'739_c4_q23_s2_H8_U[]', '719_c6_q17_s2_H7_U[[neg,[1]]]'}
def prio(o):
    c = o['cfg']; key = json.dumps([c['v'], c['c'], c['q'], c['s'], c['H1'], c['U']])
    S1 = zs.get(key, {}).get('log2_stage1_total')
    sym = log2_symmetry(SD.cfg_from_json(json.dumps(c)))['log2_sym']
    rich = S1 is not None and S1 - sym > -4
    return (0 if rich else 1, o['cost'])
todo = sorted([o for o in cfgs if o['tag'] not in skip], key=prio)
for i, o in enumerate(todo):
    if i % nW != part: continue
    out = os.path.join(ROOT, f"logs/s1z/agg/agg_{o['tag']}.jsonl")
    print('config', o['tag'], flush=True)
    run_config(SD.cfg_from_json(json.dumps(o['cfg'])), out, o['levels'])
