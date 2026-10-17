r"""Per-character Z_q-level exact screens (theory/s1z_zqscreen.CharAggScreen) over every compression multiset of the
configurations of logs/s1z/ts/configs.json (or given cfg JSON) [M]: a multiset is dead when, for some character chi of
Q that is real on the U-invariant compressions, no choice of letters and lifts with the multiset's compressions and
sum_k E_k = T1(0) meets sum_k PAF_{Z_q}(xbar_k^chi) = T^chi.
usage: python3 theory/s1z_char_batch.py OUT.jsonl [--memgb 3.5] [--tags t1,t2]"""
import sys, os, json, time, resource
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import ab_s1sd as SD
from s1z_lib import ZOrb
from s1z_zqscreen import CharAggScreen, real_characters

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
a = sys.argv[1:]; outp = a[0]
memgb = float(a[a.index('--memgb') + 1]) if '--memgb' in a else 3.5
resource.setrlimit(resource.RLIMIT_AS, (int(memgb * 2 ** 30),) * 2)
cfgs = json.load(open(os.path.join(ROOT, 'logs/s1z/ts/configs.json')))
if '--tags' in a: cfgs = [o for o in cfgs if o['tag'] in a[a.index('--tags') + 1].split(',')]
done = set()
if os.path.exists(outp):
    for l in open(outp): done.add(json.loads(l)['tag'])
for o in cfgs:
    if o['tag'] in done: continue
    t0 = time.time(); rec = dict(tag=o['tag'])
    try:
        cfg = SD.cfg_from_json(json.dumps(o['cfg'])); Z = ZOrb(cfg)
        C = Z.compressions(); cls, idx = Z.classes(C); qd = Z.class_quads(cls); cache = {}
        for q in qd: q['n0'] = Z.n0(cls, q['m'], cache)
        qd = [q for q in qd if q['n0'] > 0]
        chars = real_characters(Z)
        scr = [(k, CharAggScreen(Z, v)) for k, v in chars]; caches = [dict() for _ in scr]
        per = []
        for mi, q in enumerate(qd):
            zq = [cls[i]['rep'] for i in q['m']]; res = []
            for (k, A), cc in zip(scr, caches):
                try: ok, h = A.alive(zq, cc); res.append(int(ok))
                except (RuntimeError, MemoryError) as e: res.append(-1)
            per.append(dict(i=mi, m=list(q['m']), ordered=q['ordered'], alive=res))
        rec.update(chars=[list(k) for k, _ in chars], multisets=len(qd),
                   dead=sum(1 for p in per if 0 in p['alive']), undecided=sum(1 for p in per if 0 not in p['alive'] and -1 in p['alive']),
                   per=per)
    except (Exception, MemoryError) as e:
        rec['error'] = repr(e)[:200]
    rec['secs'] = round(time.time() - t0, 1)
    with open(outp, 'a') as f: f.write(json.dumps(rec) + "\n")
    print(json.dumps({k: v for k, v in rec.items() if k != 'per'}), flush=True)
