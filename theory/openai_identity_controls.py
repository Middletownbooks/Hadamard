#!/usr/bin/env python3
"""Validate known small cores and independently border the two m=3 controls.

No new core is constructed here. Source: Claude's independently rebuilt elliptic
controls. The PAF checker and matrix verifier import no candidate-builder code.
"""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]


def main():
    out=ROOT/'logs/openai_identity'
    out.mkdir(parents=True,exist_ok=True)
    old=json.loads((ROOT/'logs/4p3_review/elliptic_cores_for_border.json').read_text())
    selected={}
    for rec in old:
        selected.setdefault((rec['p'],rec['q']),rec)
    controls=[dict(name=r['name'],q=r['q'],m=r['p'],factors=[r['p']],sequences=r['sequences'],
                   source='logs/4p3_review/elliptic_cores_for_border.json')
              for key,r in sorted(selected.items())]
    input_file=out/'positive_small_input.json'
    input_file.write_text(json.dumps(controls,indent=2)+'\n')
    subprocess.run([sys.executable,str(ROOT/'theory/openai_identity_verify.py'),str(input_file),
                    '--expect','core','--out',str(out/'positive_small.json')],check=True)
    border_input=out/'border_input.json'
    border_input.write_text(json.dumps([r for (p,q),r in sorted(selected.items()) if p==3],indent=2)+'\n')
    with tempfile.TemporaryDirectory(prefix='openai_identity_') as scratch:
        subprocess.run([sys.executable,str(ROOT/'theory/4p3_review_border.py'),
                        '--only','elliptic','--elliptic',str(border_input),
                        '--scratch',scratch,'--out',str(out/'border_controls.json')],check=True)
    result=json.loads((out/'border_controls.json').read_text())
    assert result['all_certificates_ok'] and result['all_cores_hadamard']
    assert len(result['cores'])==2
    print(json.dumps(dict(known_small_cores=len(controls),independent_matrices=2,orders=[60,60])))


if __name__=='__main__':
    main()
