r"""Open-order screen, step 4: merge the screens and estimates into logs/open_screen/table.jsonl (one line per order,
family and best slice) and a markdown table [M gates, E estimates].

Inputs (all restart-safe logs of this screen):
  logs/open_screen/screens_s*.jsonl   exact screens A, B per slice (theory/open_screen_enum.py)
  logs/open_screen/deficit_s*.jsonl   stage-1 estimate, stage-2 H2/C2, split deficits (theory/open_screen_deficit.py)
  logs/open_screen/b0_prime_direct.jsonl   direct (unstaged) b = 0 multiplier slices at prime v (scripts/b0_open_screen.py)
  logs/open_screen/tmat_cw.jsonl      T(t) slices at the prime cofactors t of composite v (scripts/tmat_primev.py)
  report/zq_census_answer.md section 3 (Z_q all-character survivors; transcribed below)
Per staged slice:  E_tot = stage-1 log2E + E per y  (log2 expected cores in the slice);  cost = best split's
  D_def + log2(2^S-step + 2^(h_D/2)).  Regime (report/h863_lift_answer.md section 5): D_def <= 26, h_D <= 40,
  E per y >= 10, C2 <= 150.  Best slice per (order, family): tier 0 = in regime and E_tot >= 0, 1 = in regime, 2 = E_tot >= 0, 3 = rest; then lowest cost.
usage: python3 theory/open_screen_table.py [--md FILE]"""
import sys, os, json, glob, math
from collections import defaultdict, Counter
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
L = os.path.join(ROOT, 'logs/open_screen')
OPEN = json.load(open(os.path.join(ROOT, 'data/open_orders_lt10000.json')))['open']

# report/zq_census_answer.md section 3: best alive structure per order (S, r, e, log2 E_Q) and the census status
ZQ = {669: ('z3all', 'S6 r118', 111.41, 9, 'N37/N6/N3/N2 over list or time caps'),
      1115: ('z5all', 'S18 r90', 1.58, 3, 'N74, N37, N6, N3, N2 over caps; CP-SAT 60 s undecided'),
      1169: ('z7all', None, None, 0, 'closed [M]'), 1257: ('z3all', 'S45 r54', -31.43, 2, 'alive, estimate < 0'),
      1437: ('z3all', None, None, 0, 'closed [M]'), 1473: ('z3all', None, None, 0, 'closed [P, M]'),
      1713: ('z3all', 'S36 r90', -108.76, 9, 'alive, estimate < 0'),
      1719: ('z3all', 'S30 r105', -31.96, 361, 'alive, estimates < 0 (z9all, z9x3all: no finite estimate)'),
      1795: ('z5all', None, None, 0, 'closed [M]'), 1893: ('z3all', None, None, 0, 'closed [P, M]'),
      1929: ('z3all', None, None, 0, 'closed [P, M]'),
      1969: ('z11all', 'S22 r110', 23.73, 1, 'Q-level join capped; CP-SAT 600 s undecided'),
      1981: ('z7all', 'S69 r54', -36.95, 3, 'alive, estimate < 0'), 2095: ('z5all', 'S45 r90', -68.86, 3, 'alive, estimate < 0'),
      2171: ('z13all', None, None, 0, 'closed [M]'), 2215: ('z5all', 'S48 r90', -115.09, 6, 'alive, estimate < 0'),
      2315: ('z5all', 'S60 r108', None, 3, 'alive, no finite estimate'),
      2327: ('z13all', 'S33 r90', 34.95, 1, 'N89, N2 undecided at 600 s CP-SAT'),
      2429: ('z7all', None, None, 0, 'closed [M]')}


def load(pat):
    out = []
    for f in sorted(glob.glob(os.path.join(L, pat))):
        for l in open(f):
            try: out.append(json.loads(l))
            except Exception: pass
    return out


def regime(d):
    b = d.get('best') or {}
    if not b or d.get('C2') is None: return None
    t = dict(D_def=b['D_def'] <= 26, h_D=b['h_D'] <= 40, E_per_y=d['E_per_y'] >= 10, C2=d['C2'] <= 150)
    return t


def slice_desc(s):
    U = s['U']; us = ''
    if U:
        g = U[0]
        us = ' U=-1 on axes ' + str(g[1]) if g[0] == 'neg' else f" U=<{s.get('Ugen', [None])[0]}> (order {g[1]}) on axis {g[0]}"
    Gp = 'x'.join(f'Z{n}' for n in list(s['Q']) + list(s['H1'])) or '1'
    return f"G'={Gp}, q={s['q']}, s={s['s']}, f={s['f']}{us}, r={s['r']}"


def main():
    scr = {r['key']: r for r in load('screens_s*.jsonl')}
    dfc = {}
    for r in load('deficit_s*.jsonl'):
        if r['key'] not in dfc or 'error' in dfc[r['key']] or 'timeout' in dfc[r['key']]: dfc[r['key']] = r
    z2 = {r['key']: r for r in load('z2gate.jsonl')}
    s1x = {}
    for r in load('stage1_exact*.jsonl'):
        k = json.dumps(r['key']).replace(' ', '')
        A, B, C = r.get('A', {}), r.get('B', {}), r.get('C', {})
        nsol = [x.get('solutions') for x in (A, B) if 'solutions' in x]
        empty = (len(nsol) >= 1 and all(n == 0 for n in nsol) and (len(nsol) == 2 or C.get('status') == 'INFEASIBLE'))
        pl_ok = all(all(v in (True, 'FEASIBLE') for kk, v in p_.items() if kk != 'T0') for p_ in r.get('planted', []))
        if any(n and n > 0 for n in nsol) or C.get('status') == 'FEASIBLE': s1x[k] = 'stage 1 non-empty [M]'
        elif empty and pl_ok: s1x[k] = 'stage 1 empty [M]'
        else: s1x.setdefault(k, 'stage 1 undecided')
    norm = lambda k: json.dumps(json.loads(k)).replace(' ', '')
    from open_screen_enum import all_configs, cfgkey
    enum = defaultdict(list)
    for c in all_configs(): enum[(c['v'], c['fam'])].append(cfgkey(c))
    rows = []
    for N in OPEN:
        v = N // 4
        for fam in ('b0', 'b1', 'k4'):
            keys = enum.get((v, fam), [])
            if not keys and not ((fam == 'b0' and v % 2) or fam == 'b1' or (fam == 'k4' and v % 4 == 3)): continue
            st = Counter(); alive = []
            for k in keys:
                s = scr.get(k)
                if s is None: st['not run'] += 1; continue
                if s.get('A') is False: st['A dead'] += 1
                elif s.get('A') is None: st['A undecided' if 'A_reason' in s or 'timeout' in s or 'error' in s else 'A undecided'] += 1; alive.append(k)
                elif s.get('B') is False: st['B dead'] += 1
                elif z2.get(k, {}).get('Z2') is False: st['Z2 dead'] += 1
                elif z2.get(k, {}).get('wscreen') is False: st['w-screen dead'] += 1
                elif s1x.get(norm(k)) == 'stage 1 empty [M]': st['stage 1 empty [M]'] += 1
                else: st['alive (B pass)' if s.get('B') else 'alive (B capped)'] += 1; alive.append(k)
            cand = []
            for k in alive:
                d = dfc.get(k)
                if not d or 'best' not in d: continue
                s1 = (d.get('stage1') or {}).get('stage1_log2E')
                s1 = -math.inf if s1 is None else s1
                Et = s1 + d['E_per_y'] if d.get('E_per_y') is not None else -math.inf
                rg = regime(d); ok = all(rg.values()) if rg else False
                tier = 0 if (ok and Et >= 0) else 1 if ok else 2 if Et >= 0 else 3
                cand.append((tier, 0, d['best']['cost'], k, d, Et, rg))
            rec = dict(N=N, v=v, fam=fam, enumerated=len(keys), gates=dict(st),
                       estimated=len(cand), note=None if keys else 'no GP slice with 24 <= r <= 160')
            if cand:
                cand.sort(key=lambda t: t[:3])
                _, _, _, k, d, Et, rg = cand[0]; s = scr[k]; b = d['best']
                st1 = d.get('stage1') or {}
                rec.update(slice=slice_desc(s), key=k, A=s.get('A'), B=s.get('B'), Z2=z2.get(k, {}).get('Z2'),
                           wscreen=z2.get(k, {}).get('wscreen'), stage1_exact=s1x.get(norm(k)),
                           stage1_method=st1.get('stage1_method', 'ab_gauss'), row_quads=s.get('A_row_quads', [])[:3],
                           stage1_log2Ny=st1.get('log2Ny'), stage1_log2E=st1.get('stage1_log2E'),
                           H2=d['H2'], C2=d['C2'], E_per_y=d['E_per_y'], log2E_tot=None if not math.isfinite(Et) else round(Et, 2),
                           split=b['L'], C_S=b['C_S'], H_S=b['H_S'], D_def=b['D_def'], h_D=b['h_D'], Sstep=b['Sstep'],
                           cost=b['cost'], regime=rg, in_regime=bool(rg and all(rg.values())),
                           tier=cand[0][0], n_tier0=sum(1 for c in cand if c[0] == 0),
                           n_in_regime=sum(1 for c in cand if c[0] in (0, 1)),
                           n_Etot_pos=sum(1 for c in cand if c[0] in (0, 2)),
                           next_gate='exact stage 1 (2+2 join over the lag-0/row-feasible y lists)',
                           next_gate_log2cost=None if st1.get('log2Ny') is None else round(st1['log2Ny'] / 2, 1))
            rows.append(rec)
    # direct b = 0 at prime v
    best = {}
    for r in load('b0_prime_direct.jsonl'):
        if not r.get('feasible') or r.get('n_bits') is None: continue
        t = r['t']
        if t not in best or r['n_bits'] < best[t]['n_bits']: best[t] = r
    pv = {int(x) for x in open(os.path.join(L, 'prime_v.txt')).read().split(',')} if os.path.exists(os.path.join(L, 'prime_v.txt')) else set()
    seen = Counter(r.get('t') for r in load('b0_prime_direct.jsonl'))
    for v in sorted(pv):
        r = best.get(v)
        rows.append(dict(N=4 * v, v=v, fam='b0direct', slices_screened=seen.get(v, 0),
                         slice=None if not r else f"d={r['d']} ({r.get('twist')}), orbits={r['orbits']}",
                         n_bits=None if not r else r['n_bits'], log2E=None if not r else r['log2E'],
                         note='unstaged multiplier slice; 883 was built at n = 111.9 (k-tree), 919 at 91.1 (3-cycle)'))
    # Cooper-Wallis
    tm = defaultdict(list)
    for r in load('tmat_cw.jsonl'): tm[r['t']].append(r)
    from sympy import factorint, isprime
    for N in OPEN:
        v = N // 4
        if isprime(v): continue
        fs = factorint(v)
        routes = []
        for t in fs:
            if t < 100: continue
            w = v // t
            feas = [r for r in tm.get(t, []) if r.get('rowsum_feasible') and r.get('n_bits') is not None]
            bestt = min(feas, key=lambda r: r['n_bits']) if feas else None
            routes.append(dict(t=t, w=w, T_slices=len(tm.get(t, [])), best_T=None if not bestt else
                               dict(k=bestt['k'], orbits=bestt['orbits'], n_bits=bestt['n_bits'], log2E=bestt['log2E'])))
        rows.append(dict(N=N, v=v, fam='cw', routes=routes,
                         note='T(t) with Williamson-type of order w = v/t; the split with a large Williamson order is not screened'))
    for v, (fam, st, E, nal, status) in ZQ.items():
        rows.append(dict(N=4 * v, v=v, fam='zq', zq_family=fam, best_alive=st, log2E_Q=E, alive=nal, status=status,
                         source='report/zq_census_answer.md section 3'))
    with open(os.path.join(L, 'table.jsonl'), 'w') as fo:
        for r in rows: fo.write(json.dumps(r, default=float) + "\n")
    return rows


def md(rows):
    out = ['| 4v | v | family | slices (enumerated / gates) | best slice | stage-1 log2E | E/y | C2 | split | D_def | h_D | cost | log2 E_tot | regime |',
           '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    for r in rows:
        if r['fam'] not in ('b0', 'b1', 'k4'): continue
        g = ', '.join(f'{k} {n}' for k, n in sorted(r['gates'].items()))
        if 'slice' not in r:
            out.append(f"| {r['N']} | {r['v']} | {r['fam']} | {r['enumerated']} / {g or '-'} | {r.get('note') or 'none estimated'} | | | | | | | | | |"); continue
        rg = 'yes' if r['in_regime'] else 'no: ' + ','.join(k for k, ok in r['regime'].items() if not ok)
        out.append(f"| {r['N']} | {r['v']} | {r['fam']} | {r['enumerated']} / {g} | {r['slice']} | {r['stage1_log2E']} | {r['E_per_y']} | {r['C2']} | {r['split']} | {r['D_def']} | {r['h_D']} | {r['cost']} | {r['log2E_tot']} | {rg} |")
    return "\n".join(out)


if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    rows = main()
    if '--md' in sys.argv:
        open(sys.argv[sys.argv.index('--md') + 1], 'w').write(md(rows) + "\n")
    print(len(rows), Counter(r['fam'] for r in rows))
