#!/usr/bin/env python3
"""Exact per-lag binary product projection with replica concurrence.

Each occurrence of a sign in a correlation equation gets its own real copy.
The local projector satisfies every lag exactly by choosing the cheapest
agreements/disagreements; the linear projector averages all copies of a sign.
No Fourier projection, learned model, or per-block solution list is used.
"""
import argparse
import itertools
import json
from pathlib import Path
import time

import numpy as np

from openai_fresh_sphere import Model, exact_paf


class ReplicaModel:
    def __init__(self, model):
        self.model = model
        mods = model.mods
        if model.kind == 'b0':
            ids0 = model.labels
            step = len(model.orbits)
            ulist = model.orbits[1]  # orbit containing 1
        else:
            ids0 = (np.arange(20)[:, None]*7 + model.labels).reshape(mods)
            step = 140
            ulist = model.orbits[1]
        vertex = np.stack([ids0 + b*step for b in range(4)])
        self.vertex = vertex
        self.nvars = 4*step
        seen = {tuple(0 for _ in mods)}
        lags = []
        edge_ids = []
        plus = []
        constants = []
        axes = tuple(range(1, 1+len(mods)))
        for lag in itertools.product(*(range(m) for m in mods)):
            if lag in seen:
                continue
            lags.append(lag)
            for u in ulist:
                for sign in (-1, 1):
                    a = [(sign*t) % m for t,m in zip(lag, mods)]
                    a[-1] = (a[-1]*u) % mods[-1]
                    seen.add(tuple(a))
            left = vertex.ravel()
            right = np.roll(vertex, lag, axes).ravel()
            # Same vertex at both ends is the constant product +1.
            keep = left != right
            ids = np.stack((left[keep], right[keep]), axis=-1)
            target = int(model.target[lag])-int((~keep).sum())
            assert (len(ids)+target) % 2 == 0
            k = (len(ids)+target)//2
            assert 0 <= k <= len(ids)
            edge_ids.append(ids)
            plus.append(k)
            constants.append(int((~keep).sum()))
        # Group equal-sized factors to use batched partitioning.
        groups = {}
        for lag, ids, k, const in zip(lags, edge_ids, plus, constants):
            groups.setdefault((len(ids), k), []).append((lag, ids, const))
        self.groups = []
        allids = []
        cursor = 0
        for (nedges, k), rows in sorted(groups.items()):
            ids = np.stack([r[1] for r in rows])
            n = ids.size
            self.groups.append((cursor, cursor+n, ids.shape, k, [r[0] for r in rows]))
            allids.extend(ids.ravel())
            cursor += n
        self.rows = []
        if model.kind == 'k4_863':
            for ids, target in zip(vertex.reshape(80, 43), model.y.ravel()):
                k = (43+int(target))//2
                self.rows.append((cursor, cursor+43, k))
                allids.extend(ids)
                cursor += 43
        elif model.profile:
            for b, target in enumerate(model.profile):
                ids = vertex[b].ravel()
                self.rows.append((cursor, cursor+len(ids), (len(ids)+target)//2))
                allids.extend(ids)
                cursor += len(ids)
        self.ids = np.array(allids, dtype=np.int32)
        self.counts = np.bincount(self.ids, minlength=self.nvars)
        assert np.all(self.counts > 0)
        self.ncopies = cursor
        self.lags = lags

    def pa(self, x):
        a = np.empty_like(x)
        for lo, hi, shape, k, _ in self.groups:
            xx = x[:, lo:hi].reshape((x.shape[0],)+shape)
            aa = a[:, lo:hi].reshape(xx.shape)
            sp = xx[..., 0]+xx[..., 1]
            sm = xx[..., 0]-xx[..., 1]
            gain = np.abs(sp)-np.abs(sm)
            choose = np.zeros(gain.shape, dtype=bool)
            if k == shape[1]:
                choose[:] = True
            elif k:
                # Exact cardinality projection: largest k gains win.
                selected = np.argpartition(gain, shape[1]-k, axis=-1)[..., -k:]
                np.put_along_axis(choose, selected, True, axis=-1)
            first = np.where(np.where(choose, sp, sm) >= 0, 1., -1.)
            aa[..., 0] = first
            aa[..., 1] = np.where(choose, first, -first)
        for lo, hi, k in self.rows:
            aa = a[:, lo:hi]
            aa[:] = -1
            if k == hi-lo:
                aa[:] = 1
            elif k:
                selected = np.argpartition(x[:, lo:hi], hi-lo-k, axis=-1)[:, -k:]
                np.put_along_axis(aa, selected, 1., axis=-1)
        return a

    def means(self, x):
        return np.stack([np.bincount(self.ids, weights=row, minlength=self.nvars)
                         for row in x])/self.counts

    def pb(self, x):
        return self.means(x)[:, self.ids]

    def expand(self, x):
        return x[..., self.vertex]

    def assert_local(self, a):
        assert np.all(np.abs(a) == 1)
        for lo, hi, shape, k, _ in self.groups:
            aa = a[:, lo:hi].reshape((a.shape[0],)+shape)
            assert np.all((aa[..., 0]*aa[..., 1] == 1).sum(-1) == k)
        for lo, hi, k in self.rows:
            assert np.all((a[:, lo:hi] == 1).sum(-1) == k)


def run(args):
    model = Model(args.kind, args.p, args.d, args.profile, args.cert)
    r = ReplicaModel(model)
    rng = np.random.default_rng(args.seed)
    # Independent replica starts: no known core is read or injected.
    x = rng.normal(size=(args.trials, r.ncopies))
    r.assert_local(r.pa(x))
    active = np.ones(args.trials, dtype=bool)
    best = np.full(args.trials, np.inf)
    hits = []
    t0 = time.perf_counter()
    stop = 'iterations'
    for it in range(1, args.iters+1):
        a = r.pa(x)
        b = r.pb(2*a-x)
        x += args.beta*(b-a)
        if it % args.check == 0:
            votes = r.means(a)
            signs = np.where(votes >= 0, 1, -1)
            seqs = r.expand(signs)
            f = np.fft.fftn(seqs, axes=model.axes)
            err = np.rint(np.fft.ifftn(np.sum(np.abs(f)**2, axis=1), axes=model.axes).real-model.target).astype(np.int64)
            energy = np.sum(err.reshape(args.trials, -1)**2, axis=1)
            best = np.minimum(best, energy)
            for trial in np.flatnonzero(active & (energy == 0)):
                seq = seqs[trial]
                if model.kind == 'k4_863' and not np.array_equal(seq.sum(-1), model.y):
                    continue
                if model.profile and not np.array_equal(seq.sum(-1), model.profile):
                    continue
                assert model.integer_check(seq)
                hits.append(dict(trial=int(trial), iterations=it,
                                 seconds=time.perf_counter()-t0,
                                 sequences=seq.reshape(4, -1).astype(int).tolist()))
                active[trial] = False
            if not np.any(active):
                break
            if time.perf_counter()-t0 >= args.seconds:
                stop = 'seconds'
                break
    result = dict(algorithm='binary-product-replicas', kind=args.kind, p=args.p, d=args.d,
                  profile=args.profile, beta=args.beta, trials=args.trials, seed=args.seed,
                  variables=r.nvars, copies=r.ncopies, lag_factors=len(r.lags),
                  iterations=it, seconds=time.perf_counter()-t0, stop=stop,
                  solved=len(hits), best_integer_paf_energy=best.astype(int).tolist(), hits=hits,
                  certificate='direct integer PAF only; no matrix claim')
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'hits'}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--kind', choices=['b0', 'k4_863'], default='b0')
    p.add_argument('--p', type=int, default=61)
    p.add_argument('--d', type=int, default=3)
    p.add_argument('--profile', type=lambda x:list(map(int,x.split(','))))
    p.add_argument('--cert', default='data/h3452_stage1.json')
    p.add_argument('--trials', type=int, default=4)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--beta', type=float, default=.5)
    p.add_argument('--iters', type=int, default=20000)
    p.add_argument('--seconds', type=float, default=240)
    p.add_argument('--check', type=int, default=20)
    p.add_argument('--out')
    run(p.parse_args())
