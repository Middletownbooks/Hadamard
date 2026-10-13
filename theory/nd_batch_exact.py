r"""Restart-safe batch of exact compressed stage-1 decisions (theory/nd_s1_sd.py) [M].

Each job is either cyclic "m mu K" or non-cyclic "fq F p u,e[;u,e]" (Z_F x F_{p^2}, compressed along F_q).  For each
job: the alpha-quadruple count first (fast); the full exhaustive search only if it is at most --max-quads; the record
(alpha_quads, S_solutions, stage1_solutions, complete) is appended (the search stops at 200 solutions:
complete = false then means 'non-empty, not enumerated') to the output jsonl, and jobs already recorded are
skipped on restart.  Solutions (if any) are written to logs/nd/s1sols/<tag>.json.
usage: python3 theory/nd_batch_exact.py JOBFILE OUT.jsonl [--max-quads 2000000]"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nd_s1_sd as S
import nd_fq_split as FS

def problem(job):
    t = job.split()
    if t[0] == "fq":
        F, p = int(t[1]), int(t[2]); gens = [tuple(map(int, g.split(','))) for g in t[3].split(';')]
        C = FS.compressed(F, p, gens)
        UF = sorted({u % F for (u, e) in gens}); MF = [1]
        while True:
            new = sorted({x * u % F for x in MF for u in UF} | set(MF))
            if new == MF: break
            MF = new
        C['MF'] = MF
        P = S.Problem(C['m'], 0, C['q'], explicit=C); P.C_explicit = C
        return P, f"fq_F{F}_p{p}_" + "_".join(f"{u}-{e}" for u, e in gens)
    m, mu, K = map(int, t)
    return S.Problem(m, mu, K), f"cyc_{m}_{mu}_{K}"

def main(jobfile, out, max_quads):
    done = set()
    if os.path.exists(out):
        for l in open(out):
            try: done.add(json.loads(l)['job'])
            except Exception: pass
    os.makedirs("logs/nd/s1sols", exist_ok=True)
    for job in [l.strip() for l in open(jobfile) if l.strip() and not l.startswith('#')]:
        if job in done: continue
        t0 = time.time(); rec = dict(job=job)
        try:
            P, tag = problem(job)
            rec.update(F=P.F, n=P.n, f0=P.f0, orbits=len(P.reps))
            A, Q, info = P.alpha_quads(max_quads=max_quads)
            rec.update(alpha_quads=len(Q), **info)
            if len(Q) <= max_quads:
                o, sols = P.solve(max_sol=200, log=lambda x: None, quads=Q)
                rec.update(S_solutions=o['S_solutions'], stage1_solutions=o['stage1_solutions'],
                           all_verified=o['all_verified'], complete=o['complete'])
                if sols: json.dump(dict(job=job, solutions=sols[:2000]), open(f"logs/nd/s1sols/{tag}.json", "w"))
        except AssertionError as e:
            rec.update(error="not applicable: " + str(e))
        except RuntimeError as e:
            rec.update(error=str(e))
        rec['secs'] = round(time.time() - t0, 1)
        with open(out, "a") as fh: fh.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)

if __name__ == "__main__":
    a = sys.argv[1:]; mq = 2000000
    if "--max-quads" in a: i = a.index("--max-quads"); mq = int(a[i + 1]); a = a[:i] + a[i + 2:]
    main(a[0], a[1], mq)
