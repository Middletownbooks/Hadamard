"""W-word version of scripts/build_slot_pools.py: slot pools for scripts/kt_join_slots_w.py from multiplier-orbit
representatives as (N, W) masks (scripts/flat_gen_w.py output).

Slot A: the canonical blocks (one per multiplier orbit) of row sum a, at most LA of them (LA = 0: slot A is a
closure sample like the others).  Slots B, C, D: uniform samples of size L from the rotation closure of the orbit
representatives of their row sums.  Symmetry breaking (LA > 0) only pays when the samples cover most of each closure;
with low coverage, equal closure samples in all four slots give more tuples per unit of join work (the closure of a
rho-flat set is rho-flat: rotations keep F and the row sum, forensics in scripts/flat_probe.py).  Representatives
from several probe/generation directories are merged (np.unique), so longer generation only adds orbits.
With --closure: every slot is the whole rotation closure of its row sum when that has at most L blocks (else a
sample of L), as scripts/build_slot_pools.py --closure; slots with equal row sums share one pool (symlinks).
usage: python3 scripts/build_slot_pools_w.py V D a,b,c,d OUTDIR L LA DIR1 [DIR2 ...] [--closure]"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import flatw as FW


def closure_sample_w(reps, e, n, rng):
    i = rng.integers(0, len(reps), n); s = rng.integers(0, e, n); out = np.empty((n, reps.shape[1]), np.uint64)
    for s_ in range(e):
        m = s == s_
        if m.any(): out[m] = FW.rotate_w(np.ascontiguousarray(reps[i[m]]), s_, e)
    return FW.unique_w(out)


def full_closure_w(reps, e):
    return FW.unique_w(np.concatenate([FW.rotate_w(np.ascontiguousarray(reps), s, e) for s in range(e)]))


def main():
    a = sys.argv[1:]; clos = '--closure' in a; a = [x for x in a if x != '--closure']
    v, d = int(a[0]), int(a[1]); prof = [int(x) for x in a[2].split(',')]; out = a[3]; L = int(a[4]); LA = int(a[5]); dirs = a[6:]
    e = (v - 1) // d; os.makedirs(out, exist_ok=True); rng = np.random.default_rng(2024)
    reps = {}
    for s in set(prof):
        parts = [np.load(os.path.join(dd, f'orbits_{s}.npy')) for dd in dirs if os.path.exists(os.path.join(dd, f'orbits_{s}.npy'))]
        reps[s] = FW.unique_w(np.concatenate([p.reshape(len(p), -1) for p in parts])) if parts else np.zeros((0, 1), np.uint64)
    info = {}; first = {}
    for slot, s in zip('ABCD', prof):
        if clos and s in first:                  # closure mode: equal row sums share one pool (a symlink; the disk join
            dst = os.path.join(out, f'slot_{slot}.npy')   # then shares the prep of equal pools in equal key order)
            if os.path.lexists(dst): os.remove(dst)
            os.symlink(f'slot_{first[s]}.npy', dst); info[slot] = dict(info[first[s]], same_as=first[s]); continue
        first.setdefault(s, slot)
        if clos and len(reps[s]) * e <= L:
            x = full_closure_w(reps[s], e)
        elif slot == 'A' and LA > 0 and not clos:
            x = reps[s] if len(reps[s]) <= LA else reps[s][np.sort(rng.choice(len(reps[s]), LA, replace=False))]
        else:
            x = closure_sample_w(reps[s], e, L, rng)
        np.save(os.path.join(out, f'slot_{slot}.npy'), x)
        info[slot] = dict(rowsum=s, size=len(x), log2=round(float(np.log2(max(len(x), 1))), 2), orbits_available=len(reps[s]))
    json.dump(dict(v=v, d=d, profile=prof, slots=info, sources=dirs, closure=clos), open(os.path.join(out, 'slots.json'), 'w'), indent=1)
    print(json.dumps(info))


if __name__ == '__main__':
    main()
