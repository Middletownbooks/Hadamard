r"""Coset-Parseval screen of multiplier slices [P for the criterion, M for its evaluation].

Criterion (theorem of theory/b599x4_parseval.py).  For a slice group S acting on G, a nontrivial subgroup N <= H
(H = G for b = 0, 1; the index-c subgroup for b >= 3), and C_N = the union of the S-orbits O with O + N != O:
every character nontrivial on N needs sum_k |x_k^|^2 = 4v, but sees only x on C_N, so
    sum_k sum_{chi|N != 1} |x_k^(chi)|^2 <= 4 |G| |C_N| (1 - 1/|N|^2)  (|N| odd),   <= 4 |G| |C_N|  (|N| even)
against the required |G| (1 - 1/|N|) 4v.  The slice has no core when
    |C_N| (1 + 1/|N|) < v  (|N| odd)      or      |C_N| < (1 - 1/|N|) v  (|N| even).
The screen tests the subgroups N <= H of prime order: an orbit that is stable under N is stable under every subgroup
of N, so the prime-order N have the smallest C_N; they are the case met by multipliers u = 1 (mod p^(e-1)) on an axis
Z_(p^e), e >= 2, whose unit orbits are cosets of p^(e-1) Z_(p^e).

Families screened: theory/ab_cbatch.configs (c = 4, 6, 8 at the open orders), theory/ab_sweep.configs_bge3 (maxU = 2)
and configs_b0, the b = 1 multiplier sweeps logs/msweep_record/m*.jsonl (cyclic Z_m, multiplier subgroups), the 112
configurations of logs/s1z/table.json, and the slices of known cores as negative controls (they must pass).
usage: python3 theory/b599x4_cosetscreen.py OUT.jsonl [--summary FILE]"""
import sys, os, json, math, itertools
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W
from ab_lib import AG, closure, factorint

OPEN = [571, 599, 643, 659, 719, 739, 751, 823, 839, 859, 863, 883, 907, 917, 919, 947, 955, 971, 991]


def orbit_labels(G, Slist):
    perms = [G.act(u) for u in Slist]
    orb = -np.ones(G.n, np.int64); k = 0
    for g in range(G.n):
        if orb[g] >= 0: continue
        for P in perms: orb[P[g]] = k
        k += 1
    return orb, k


def prime_subgroups(G, Hmask, maxsub=4000):
    """subgroups N <= H of prime order (as element lists), and the cyclic subgroups of order p^2 when cheap"""
    out = []; seen = set()
    Hel = np.nonzero(Hmask)[0]
    for g in Hel:
        if g == 0: continue
        o = G.order(int(g))
        if o not in factorint(o) or factorint(o)[o] != 1:
            continue                                    # only prime order generators here
        els = {0}; x = int(g)
        while x != 0: els.add(x); x = int(G.add(x, g))
        key = frozenset(els)
        if key in seen: continue
        seen.add(key); out.append(sorted(els))
        if len(out) > maxsub: break
    return out


def screen(G, Slist, Hmask, v):
    orb, r = orbit_labels(G, Slist)
    best = None
    for N in prime_subgroups(G, Hmask):
        stable_pt = np.ones(G.n, bool)
        for n in N[1:]:
            sh = G.idx(G.co + G.co[n])
            stable_pt &= (orb[sh] == orb)
        # an orbit is N-stable iff all its points map into it under every n
        bad_orbits = np.zeros(r, bool); np.logical_or.at(bad_orbits, orb, ~stable_pt)
        C = int(bad_orbits[orb].sum()); nN = len(N)
        lhs = C * (1 + 1 / nN) if nN % 2 else C / (1 - 1 / nN)
        margin = v / lhs if lhs else math.inf               # > 1: dead
        rec = dict(N_order=nN, N_gen=[int(x) for x in G.co[N[1]]], C=C, margin=round(margin, 4), dead=bool(lhs < v))
        if best is None or rec['margin'] > best['margin']: best = rec
    return dict(orbits=int(r), G=int(G.n), best=best, dead=bool(best and best['dead']))


def screen_cfg(cfg):
    sl = W.build(cfg)
    Hmask = sl.Hmask if cfg['b'] >= 3 else np.ones(sl.G.n, bool)
    return screen(sl.G, sl.S, Hmask, cfg['v'])


def key_bge3(cfg):
    return json.dumps([cfg['v'], cfg['c'], list(cfg['Q']), list(cfg['H1']), cfg['q'], cfg['s'],
                       [list(u) if u[0] != 'neg' else ['neg', list(u[1])] for u in cfg['U']]])


def families():
    from ab_cbatch import configs as cb_configs
    seen = set()
    for v in OPEN:
        for c in (4, 6, 8):
            for cfg in cb_configs(v, c):
                k = key_bge3(cfg)
                if k in seen: continue
                seen.add(k); yield 'cbatch', cfg
        for cfg in W.configs_bge3(v, cs=(4, 6, 8), maxU=2):
            cfg = dict(cfg, Q=list(cfg['Q']), H1=list(cfg['H1']), mods=list(cfg['mods']))
            k = key_bge3(cfg)
            if k in seen: continue
            seen.add(k); yield 'sweep_bge3', cfg
        for cfg in W.configs_b0(v, maxU=2):
            cfg = dict(cfg, Q=[], H1=list(cfg['H1']), mods=list(cfg['mods']))
            k = key_bge3(cfg)
            if k in seen: continue
            seen.add(k); yield 'sweep_b0', cfg


def b1_records(path='logs/msweep_record'):
    for fn in sorted(os.listdir(path)):
        if not fn.endswith('.jsonl'): continue
        lines = [json.loads(l) for l in open(os.path.join(path, fn))]
        m = lines[0]['m']
        for r in lines[1:]:
            if 'gens' in r: yield m, r


def table_cfgs(path='logs/s1z/table.json'):
    for t in json.load(open(path)):
        v, c, q, s, H1, U = t['key']
        Q = list(t['Q']); h = (v + 1) // c - 1
        yield t, dict(v=v, c=c, b=c - 1, h=h, m=c * h, Q=Q, H1=list(H1), q=q, s=s, f=(q - 1) // s,
                      mods=Q + list(H1) + [q], nQ=len(Q),
                      U=[tuple(u) if u[0] != 'neg' else ('neg', tuple(u[1])) for u in U])


def controls():
    """slices of known cores (must not be killed)"""
    out = []
    c347 = json.load(open('logs/ab/cfg_347_s6.json'))
    c347['U'] = [tuple(u) if u[0] != 'neg' else ('neg', tuple(u[1])) for u in c347['U']]
    out.append(('347 (published b = 5 core, H(1388))', c347))
    return out


def b1_screen(m, gens, v):
    G = AG((m,))
    S = closure(G, [(g,) for g in gens])
    return screen(G, S, np.ones(G.n, bool), v)


if __name__ == "__main__":
    a = sys.argv[1:]; outp = a[0]
    summ = a[a.index('--summary') + 1] if '--summary' in a else None
    fo = open(outp, 'w'); tally = defaultdict(lambda: [0, 0]); killed = []
    # negative controls first
    for name, cfg in controls():
        r = screen_cfg(cfg); r.update(family='control', name=name)
        fo.write(json.dumps(r) + '\n'); assert not r['dead'], name
    # known b = 1 and b = 0 cores' slices (cyclic multiplier slices): 523 x181, 647 x353, 883 order-21, 919
    for name, m, gens, v in (('523 x181 (H(2092))', 522, [181], 523), ('647 x353 (H(2588))', 646, [353], 647)):
        r = b1_screen(m, gens, v); r.update(family='control', name=name)
        fo.write(json.dumps(r) + '\n'); assert not r['dead'], name
    for fam, cfg in families():
        try:
            r = screen_cfg(cfg)
        except Exception as e:
            r = dict(error=repr(e)[:200])
        r.update(family=fam, key=json.loads(key_bge3(cfg)))
        fo.write(json.dumps(r) + '\n')
        if 'dead' in r:
            tally[fam][0] += 1; tally[fam][1] += int(r['dead'])
            if r['dead']: killed.append(r['key'])
    for m, rec in b1_records():
        v = m + 1
        r = b1_screen(m, rec['gens'], v)
        r.update(family='msweep_b1', key=[v, rec['gens'], rec.get('order')], log2E_record=rec.get('log2E'))
        fo.write(json.dumps(r) + '\n')
        tally['msweep_b1'][0] += 1; tally['msweep_b1'][1] += int(r['dead'])
        if r['dead']: killed.append(r['key'])
    tab = []
    for t, cfg in table_cfgs():
        r = screen_cfg(cfg); r.update(family='table112', key=t['key'], verdict_before=t.get('verdict'))
        fo.write(json.dumps(r) + '\n'); tab.append(r)
        tally['table112'][0] += 1; tally['table112'][1] += int(r['dead'])
    fo.close()
    S = dict(tally={k: dict(screened=v[0], killed=v[1]) for k, v in tally.items()},
             table112_killed=[dict(key=r['key'], verdict_before=r['verdict_before'], best=r['best']) for r in tab if r['dead']])
    print(json.dumps(S, indent=1))
    if summ: json.dump(S, open(summ, 'w'), indent=1)
