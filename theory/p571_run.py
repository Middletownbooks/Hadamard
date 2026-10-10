r"""Driver for the hierarchical solver (theory/p571_hsolve.py): stage 1 -> stage 2 -> bordered GS completion.

  python3 theory/p571_run.py m gen [gen ...] [--L 4096] [--rounds 1] [--maxy N] [--sub S] [--stop 1] [--tag name]

Stage 1: theta_5-profiled (stage1c) when 5 | q and M maps onto Z_5^*, otherwise plain L1 -> L2 (lists of up to L
distinct X1 classes per block; exhaustive when the class space is smaller).  Stage 2 on every y found (up to maxy),
S-level lists capped at S classes per block (None = all).  Every core is checked exactly (PAF) and completed by
theory/p571_finish.bordered_gs (H H^T = N I); the first Hadamard core is written to data/.
Writes one JSON line per y to stdout and a summary."""
import sys, os, json, time, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p571_hsolve import Inst, stage1, stage1c, stage2, y_to_x_check
from p571_finish import bordered_gs, write_core

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main(m, gens, L=4096, rounds=1, maxy=10 ** 9, sub=None, stop=1, tag='', seed=0, n_prof=10):
    inst = Inst(m, gens); q = inst.q
    t0 = time.time()
    use_c = (q % 5 == 0) and len({t % 5 for t in inst.S.Gq}) == 4
    if use_c: ys = stage1c(inst, L=L, n_profiles=n_prof, seed=seed, log=True)
    else: ys, _ = stage1(inst, L=L, rounds=rounds, seed=seed, log=True)
    if ys:
        Y = np.unique(np.array(ys).reshape(len(ys), -1), axis=0).reshape(-1, 4, 2, inst.R)
    else: Y = np.zeros((0, 4, 2, inst.R), np.int64)
    s1ok = all(y_to_x_check(inst, y) for y in Y[:200])
    print(json.dumps(dict(stage1_distinct_ys=len(Y), S1_exact=bool(s1ok), secs=round(time.time() - t0, 1))), flush=True)
    found = 0; hadamard = 0; t1 = time.time(); n_s = 0
    for i, y in enumerate(Y[:maxy]):
        rec, cores = stage2(inst, y, sub=sub, seed=seed + i, tbits=20)
        n_s += rec['s_solutions'] > 0
        if cores:
            found += len(cores)
            H, r = bordered_gs(cores[0])
            print(json.dumps(dict(y=i, **rec, finish=r)), flush=True)
            if r.get('hadamard'):
                hadamard += 1
                if hadamard == 1:
                    path = os.path.join(ROOT, 'data', f'core_v{inst.v}_b1_p571method{tag}.txt')
                    write_core(path, cores[0], f'v={inst.v} N={4 * inst.v} b=1 m={m} cyclic Z{m}; slice M=<{",".join(map(str, gens))}> '
                               f'(Z6 fixed pointwise); theory/p571_hsolve.py (stage 1 L1/L2, stage 2 S/D); verified by theory/p571_finish.py')
                if hadamard >= stop: break
    print(json.dumps(dict(summary=True, v=inst.v, gens=gens, ys_processed=min(len(Y), maxy, i + 1 if len(Y) else 0),
                          ys_with_S_solutions=n_s, cores=found, hadamard=hadamard, stage2_secs=round(time.time() - t1, 1))), flush=True)

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--L": "4096", "--rounds": "1", "--maxy": str(10 ** 9), "--sub": "0", "--stop": "1",
                             "--tag": "", "--seed": "0", "--np": "10"}
    pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    main(pos[0], pos[1:], L=int(opt['--L']), rounds=int(opt['--rounds']), maxy=int(opt['--maxy']),
         sub=(int(opt['--sub']) or None), stop=int(opt['--stop']), tag=opt['--tag'], seed=int(opt['--seed']), n_prof=int(opt['--np']))
