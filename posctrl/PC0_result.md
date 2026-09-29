# PC0 (known-solution slice v=47, h=7, p=7, u=29): FAIL (unevaluable)
The blind run could not complete: at p = h = 7 the compression is only the coset sums, so every block's lift list is
the whole symmetric space, and the exact join ran more than 30 min without finishing. Per the registered rules, an
unevaluable bar is a FAIL.

Targeted checks on the known solution (data/rung_c6_v47_sym.txt), via scripts/bge3_posctrl_known47.py:
- its compression is in the exact stage-A block lists for its shift (0, -2, 4, 6), for all four blocks: PASS;
- per-class residue filter and pairwise joint filter at M = 16 and 64: all pass (no false rejection): PASS;
- PAF equations, and assembly with the pipeline's template code into a 188 x 188 matrix: H H^T = 188 I: PASS.
