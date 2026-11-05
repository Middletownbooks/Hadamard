#!/usr/bin/env python3
"""Independent exact PAF and primitive-character energy checks. No builder imports.

The energy is computed in the group domain using integer centering matrices,
not numerical Fourier transforms or the construction's predicted formula.
"""
import argparse
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np


def paf_check(rec):
    factors = rec['factors']
    mods = ([2, 2] if rec['q'] == 'v4' else [4]) + factors
    coords = np.array(list(itertools.product(*(range(d) for d in mods))), dtype=np.int64)
    rows = np.array(rec['sequences'], dtype=np.int64)
    N = math.prod(mods)
    assert rows.shape == (4, N) and np.isin(rows, (-1, 1)).all()
    paf = []
    for co in coords:
        shifted = (coords + co) % mods
        ix = np.ravel_multi_index(shifted.T, mods)
        paf.append(int((rows * rows[:, ix]).sum()))
    qdim = 2 if rec['q'] == 'v4' else 1
    target = np.where(np.all(coords[:, :qdim] == 0, axis=1), -12, 4)
    target[0] = 4 * N
    bad = np.flatnonzero(np.array(paf) != target)
    return dict(ok=len(bad) == 0, lags=N, bad_lags=len(bad),
                first_bad=None if not len(bad) else dict(lag=coords[bad[0]].tolist(),
                                                       actual=paf[bad[0]], wanted=int(target[bad[0]])))


def primitive_energy(rec):
    factors = rec['factors']
    m = math.prod(factors)
    count = math.prod(d - 1 for d in factors)
    cells = np.array(rec['sequences'], dtype=np.int64).reshape(16, *factors)
    # Apply m P = tensor_j (d_j I - J); P projects onto nonprincipal
    # characters in every factor. No eigenvalues or FFT rounding are used.
    scaled = cells.copy()
    for axis, d in enumerate(factors, 1):
        scaled = d * scaled - scaled.sum(axis=axis, keepdims=True)
    energies = [Fraction(int((a * b).sum()), count) for a, b in zip(cells, scaled)]
    total = sum(energies)
    return dict(cells=[str(e) for e in energies], total=str(total),
                required=16 * m + 12, excess=str(total - (16 * m + 12)),
                bound_16m_plus_16_holds=total >= 16 * m + 16,
                equality=total == 16 * m + 16)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input')
    ap.add_argument('--out', required=True)
    ap.add_argument('--expect', choices=['reject-tensor', 'core'], required=True)
    a = ap.parse_args()
    data = json.loads(Path(a.input).read_text())
    results = []
    for rec in data:
        item = {k: rec[k] for k in ('name', 'q', 'factors', 'm')}
        item['paf'] = paf_check(rec)
        item['primitive_energy'] = primitive_energy(rec)
        item['sequence_sha256_int8'] = hashlib.sha256(np.array(rec['sequences'], np.int8).tobytes()).hexdigest()
        if a.expect == 'reject-tensor':
            assert not item['paf']['ok'], item
            assert item['primitive_energy']['bound_16m_plus_16_holds'], item
        else:
            assert item['paf']['ok'], item
            assert item['primitive_energy']['excess'] == '0', item
        results.append(item)
    out = dict(script='theory/openai_identity_verify.py', expectation=a.expect,
               cases=len(results), all_expected=True, results=results)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(dict(cases=len(results), all_expected=True, output=a.out)))


if __name__ == '__main__':
    main()
