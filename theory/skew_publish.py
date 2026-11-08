"""Write pub/skew/<order>/ for a verified certificate (prompt §4 step 4).
  python3 theory/skew_publish.py CERT.json "<structure line>" "<search line>" "<UTC of the hit>"
Bundle: certificate.json (the verify_skew format), sequences.txt, verify.py (= scripts/verify_skew.py), verify_output.txt,
README.md (construction, search, verification, novelty status)."""
import json, os, shutil, subprocess, sys, hashlib
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cert_path, structure, search, utc = sys.argv[1:5]
c = json.load(open(cert_path)); v = c['v']; n = 4 * v
d = os.path.join(ROOT, 'pub', 'skew', str(n)); os.makedirs(d, exist_ok=True)
json.dump(c, open(os.path.join(d, 'certificate.json'), 'w'), indent=1)
with open(os.path.join(d, 'sequences.txt'), 'w') as f:
    f.write(f'# skew Goethals-Seidel quadruple over G = {" x ".join("Z%d" % m for m in c["mods"])} (element index = mixed radix)\n')
    f.write(f'# invariant under the automorphisms (matrices on coordinate vectors): {c["gens"]}\n')
    for k in 'ABCD': f.write(f'{k}={c[k]}\n')
shutil.copy(os.path.join(ROOT, 'scripts', 'verify_skew.py'), os.path.join(d, 'verify.py'))
out = subprocess.run([sys.executable, os.path.join(d, 'verify.py'), os.path.join(d, 'certificate.json')], capture_output=True, text=True)
open(os.path.join(d, 'verify_output.txt'), 'w').write(out.stdout + out.stderr)
lines = [json.loads(l) for l in out.stdout.splitlines() if l.startswith('{')]
ok = out.returncode == 0
rows = lines[0]['rowsums'] if lines else {}
sha = lines[-1].get('sha256_of_packed_rows') if lines else None
readme = f"""# Skew-Hadamard matrix of order {n} = 4 · {v}

**Status:** {'verified' if ok else 'VERIFICATION FAILED'}; {v} is one of the 630 values of the Đoković–Seberry census of undecided skew orders 4v
still open after Kumar's repository at d885022 and Schumer's 919 (`data/skew_open_census.json`).

## Construction
A skew Goethals–Seidel array (Goethals–Seidel 1970; Wallis): four ±1 functions A, B, C, D on an abelian group G of order
{v}, A of skew type (A(0) = 1, A(−g) = −A(g)), with Σ PAF = 4v·δ; blocks are group-developed, X_(g,h) = X(h − g), and R
is the matrix of g ↦ −g:

    H = [[A, BR, CR, DR], [−BR, A, DᵀR, −CᵀR], [−CR, −DᵀR, A, BᵀR], [−DR, CᵀR, −BᵀR, A]]

Then H Hᵀ = {n} I and H + Hᵀ = 2I.

- {structure}
- Row sums: {rows}.

## Search
{search}
Hit logged at {utc} (UTC).

## Verification
`python3 verify.py certificate.json` (shares no code with the search): integer PAF on every lag, skew first block,
invariance under every generator, H Hᵀ = {n} I (float64 BLAS on ±1 entries, exact), H + Hᵀ = 2I.
Output in `verify_output.txt`. SHA-256 of the packed rows of H: `{sha}`.

## Novelty
Checked on 2026-11-08: the census list (Đoković–Seberry, May 2025, via Kumar's `search/data/skew_census_undecided.txt`);
Kumar's `skew/README.md` at d885022 (no order {n}); his NOVELTY.md (Schumer's deposit: 3676, 5436, 6644, 7852; EKLO
arXiv:2610.27175 v1). Not checked: EKLO v2 and Zenodo directly (refused by this session's network policy). Status:
new as far as these sources show; to be confirmed by the owner before any claim is posted.
"""
open(os.path.join(d, 'README.md'), 'w').write(readme)
print(json.dumps(dict(bundle=os.path.relpath(d, ROOT), verified=ok, sha=sha)))
