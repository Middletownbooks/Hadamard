#!/usr/bin/env python3
"""Fixed, coupled spectral-sphere RRR. No imports from the earlier S4 solver.

This is a heuristic, not a complete search. FFTs drive it; direct integer
correlations certify any output. A core is not a verified Hadamard matrix.
Run examples are in report/openai_fresh_answer.md.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np


def prime_orbits(p, d):
    assert p > 2 and all(p % i for i in range(2, math.isqrt(p) + 1))
    assert (p - 1) % d == 0
    for u in range(1, p):
        subgroup = sorted({pow(u, j, p) for j in range(d)})
        if len(subgroup) == d and pow(u, d, p) == 1:
            break
    seen = {0}
    orbits = [[0]]
    for i in range(1, p):
        if i not in seen:
            o = sorted({i * u % p for u in subgroup})
            seen.update(o)
            orbits.append(o)
    return orbits


def exact_paf(x, mods):
    """Direct int64 translation sums; no FFT in certification."""
    x = np.asarray(x, dtype=np.int64).reshape((4,) + tuple(mods))
    out = np.empty(mods, dtype=np.int64)
    axes = tuple(range(1, x.ndim))
    for lag in itertools.product(*(range(m) for m in mods)):
        out[lag] = np.sum(x * np.roll(x, lag, axes), dtype=np.int64)
    return out


class Model:
    def __init__(self, kind, p=61, d=3, profile=None, cert=None):
        self.kind = kind
        self.profile = profile
        if kind == 'b0':
            self.mods = (p,)
            self.orbits = prime_orbits(p, d)
            self.d = d
            self.nonzero = np.array(self.orbits[1:])
            self.labels = np.empty(p, int)
            for j, o in enumerate(self.orbits):
                self.labels[o] = j
            self.target = np.zeros(p, np.int64)
            self.target[0] = 4 * p
            self.options = []
            if profile:
                assert sum(a*a for a in profile) == 4*p
                for a in profile:
                    opts = [(e, ((a-e)//d + len(self.orbits)-1)//2)
                            for e in (-1, 1)
                            if (a-e) % d == 0
                            and ((a-e)//d + len(self.orbits)-1) % 2 == 0
                            and 0 <= ((a-e)//d + len(self.orbits)-1)//2 <= len(self.orbits)-1]
                    if not opts:
                        raise ValueError('unattainable row sum')
                    self.options.append(opts)
        elif kind == 'k4_863':
            self.mods = (2, 2, 5, 43)
            self.orbits = prime_orbits(43, 7)
            self.nonzero = np.array(self.orbits[1:])
            self.labels = np.empty(43, int)
            for j, o in enumerate(self.orbits):
                self.labels[o] = j
            y = np.array(json.loads(Path(cert).read_text())['y'], np.int64).reshape(4, 2, 2, 5)
            self.e = np.where(y % 14 == 1, 1, -1)
            assert np.all((y-self.e) % 14 == 0)
            self.nplus = ((y-self.e)//7 + 6)//2
            assert np.all((0 <= self.nplus) & (self.nplus <= 6))
            self.y = y
            self.target = np.full(self.mods, 4, np.int64)
            self.target[0, 0, :, :] = -12
            self.target[0, 0, 0, 0] = 3440
            assert np.array_equal(exact_paf(y,self.mods[:-1]),self.target.sum(axis=-1))
        else:
            raise ValueError(kind)
        self.axes = tuple(range(-len(self.mods), 0))
        target_fft = np.fft.fftn(self.target)
        assert np.max(np.abs(target_fft.imag)) < 1e-7
        assert np.min(target_fft.real) >= -1e-7
        self.radius = np.sqrt(np.maximum(0., target_fft.real))
        self.n = math.prod(self.mods)

    def pa(self, x):
        scores = x[..., self.nonzero].sum(axis=-1)
        if self.kind == 'b0':
            if self.profile is None:
                vals = np.concatenate((np.where(x[..., :1] >= 0, 1., -1.),
                                       np.where(scores >= 0, 1., -1.)), axis=-1)
            else:
                vals = np.empty(x.shape[:-1] + (len(self.orbits),))
                for b, opts in enumerate(self.options):
                    order = np.argsort(scores[:, b, :], axis=-1, kind='stable')
                    best = np.full(x.shape[0], -np.inf)
                    for e, k in opts:
                        signs = -np.ones_like(scores[:, b, :])
                        if k:
                            np.put_along_axis(signs, order[:, -k:], 1., axis=-1)
                        obj = e*x[:, b, 0] + np.sum(scores[:, b, :]*signs, axis=-1)
                        win = obj > best
                        vals[win, b, 0] = e
                        vals[win, b, 1:] = signs[win]
                        best[win] = obj[win]
        else:
            order = np.argsort(scores, axis=-1, kind='stable')
            ranks = np.argsort(order, axis=-1, kind='stable')
            signs = np.where(ranks >= (6-self.nplus)[None, ..., None], 1., -1.)
            vals = np.concatenate((np.broadcast_to(self.e, x.shape[:-1])[..., None], signs), axis=-1)
        return vals[..., self.labels]

    def pb(self, x):
        f = np.fft.fftn(x, axes=self.axes)
        norms = np.sqrt(np.sum(np.abs(f)**2, axis=1))
        scale = np.divide(self.radius, norms, out=np.zeros_like(norms), where=norms > 1e-15)
        f *= scale[:, None]
        # A zero vector has non-unique projection. A real first-axis tie choice
        # respects conjugate symmetry and the common multiplier invariance.
        zero = norms <= 1e-15
        f[:, 0] = np.where(zero, self.radius, f[:, 0])
        return np.fft.ifftn(f, axes=self.axes).real

    def integer_check(self, a):
        assert np.all(np.abs(a) == 1)
        if self.kind == 'k4_863':
            assert np.array_equal(a.sum(axis=-1), self.y)
        if self.profile:
            assert np.array_equal(a.sum(axis=-1), self.profile)
        return np.array_equal(exact_paf(a, self.mods), self.target)


def run(args):
    model = Model(args.kind, args.p, args.d, args.profile, args.cert)
    rng = np.random.default_rng(args.seed)
    x = rng.normal(size=(args.trials, 4) + model.mods)
    # Begin in the multiplier invariant real subspace.
    x = model.pa(x) + .1*model.pa(rng.normal(size=x.shape))
    active = np.ones(args.trials, bool)
    best = np.full(args.trials, np.inf)
    hits = []
    t0 = time.perf_counter()
    stop = 'iterations'
    for it in range(1, args.iters + 1):
        if args.direction == 'Afirst':
            a = model.pa(x)
            b = model.pb(2*a - x)
            change = b-a
        else:
            b = model.pb(x)
            a = model.pa(2*b - x)
            change = a-b
        gap = np.max(np.abs(a-b).reshape(args.trials, -1), axis=1)
        x += args.beta * change
        if it % args.check == 0 or np.any(active & (gap < 1e-7)):
            f = np.fft.fftn(a, axes=model.axes)
            paf = np.fft.ifftn(np.sum(np.abs(f)**2, axis=1), axes=model.axes).real
            err = np.rint(paf-model.target).astype(np.int64)
            energy = np.sum(err.reshape(args.trials, -1)**2, axis=1)
            best = np.minimum(best, energy)
            for k in np.flatnonzero(active & (energy == 0)):
                assert model.integer_check(a[k])
                hits.append(dict(trial=int(k), iterations=it,
                                 seconds=time.perf_counter()-t0,
                                 rowsums=a[k].reshape(4, -1).sum(axis=1).astype(int).tolist(),
                                 sequences=a[k].reshape(4, -1).astype(int).tolist()))
                active[k] = False
            if not np.any(active):
                break
        if it % 100 == 0 and time.perf_counter()-t0 >= args.seconds:
            stop = 'seconds'
            break
    result = dict(algorithm='fixed-coupled-sphere-RRR', direction=args.direction, kind=args.kind, p=args.p, d=args.d,
                  profile=args.profile, beta=args.beta, seed=args.seed, trials=args.trials,
                  iterations=it, stop=stop, seconds=time.perf_counter()-t0,
                  best_integer_paf_energy=best.astype(int).tolist(),
                  solved=len(hits), hits=hits, certificate='direct integer PAF only; no matrix claim')
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'hits'}), flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--kind', choices=['b0', 'k4_863'], default='b0')
    ap.add_argument('--p', type=int, default=61)
    ap.add_argument('--d', type=int, default=3)
    ap.add_argument('--profile', type=lambda x: list(map(int, x.split(','))))
    ap.add_argument('--cert', default='data/h3452_stage1.json')
    ap.add_argument('--trials', type=int, default=8)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--iters', type=int, default=20000)
    ap.add_argument('--seconds', type=float, default=240)
    ap.add_argument('--beta', type=float, default=.5)
    ap.add_argument('--direction', choices=['Afirst', 'Bfirst'], default='Afirst')
    ap.add_argument('--check', type=int, default=20)
    ap.add_argument('--out')
    run(ap.parse_args())
