"""Slot pools for scripts/kt_join_slots.py from multiplier-orbit representatives (scripts/flat_probe.py output).

Slot A: the canonical blocks (one per multiplier orbit) of row sum a, at most LA of them.  Slots B, C, D: uniform
samples of size L from the rotation closure of the orbit representatives of their row sums (the closure of a
rho-flat set is rho-flat: rotations keep F and the row sum, forensics in scripts/flat_probe.py).  Representatives
from several probe/generation directories are merged (np.unique), so longer generation only adds orbits.
With --closure (no symmetry breaking): every slot is its row sum's rotation closure, the whole closure when it has at
most L blocks, else a uniform sample of L.  Each core then sits in the pools in all e rotations (landing in different
buckets), which multiplies hits per bucket by about e for about 1.4x the bucket time: better for a first hit.
With --revcanon: keep one block of each reversal pair {x, x(-i)} (the least of w and its rotation by e/2).  Reversing
one block keeps its PAF, so in closure pools every core sits in its bucket as a cluster of 2^4 (x A<->B) identical-key
tuples; reversal-canonical slots keep one tuple of each cluster (per A<->B) at half the pool size [P].
usage: python3 scripts/build_slot_pools.py V D a,b,c,d OUTDIR L LA DIR1 [DIR2 ...] [--closure] [--revcanon]"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from flat_probe import canonical, closure_sample, rotate


def full_closure(reps, e):
    """all rotations of the orbit representatives (deduplicated: periodic orbits are rare)"""
    return np.unique(np.concatenate([rotate(reps, s, e) for s in range(e)]))


def main():
    a = sys.argv[1:]; clos = '--closure' in a; rc = '--revcanon' in a; a = [x for x in a if x not in ('--closure', '--revcanon')]
    v, d = int(a[0]), int(a[1]); prof = [int(x) for x in a[2].split(',')]; out = a[3]; L = int(a[4]); LA = int(a[5]); dirs = a[6:]
    e = (v - 1) // d; os.makedirs(out, exist_ok=True); rng = np.random.default_rng(2024)
    reps = {}
    for s in set(prof):
        parts = [np.load(os.path.join(dd, f'orbits_{s}.npy')) for dd in dirs if os.path.exists(os.path.join(dd, f'orbits_{s}.npy'))]
        reps[s] = np.unique(np.concatenate(parts)) if parts else np.zeros(0, np.uint64)
    info = {}
    for slot, s in zip('ABCD', prof):
        if clos:
            x = full_closure(reps[s], e) if len(reps[s]) * e <= L else closure_sample(reps[s], e, L, rng)
        elif slot == 'A':
            x = reps[s] if len(reps[s]) <= LA else np.sort(reps[s][rng.choice(len(reps[s]), LA, replace=False)])
        else:
            x = closure_sample(reps[s], e, L, rng)
        if rc: x = x[x <= rotate(x, e // 2, e)]
        np.save(os.path.join(out, f'slot_{slot}.npy'), x)
        info[slot] = dict(rowsum=s, size=len(x), log2=round(float(np.log2(max(len(x), 1))), 2), orbits_available=len(reps[s]))
    json.dump(dict(v=v, d=d, profile=prof, slots=info, sources=dirs, closure=clos, revcanon=rc), open(os.path.join(out, 'slots.json'), 'w'), indent=1)
    print(json.dumps(info))


if __name__ == '__main__':
    main()
