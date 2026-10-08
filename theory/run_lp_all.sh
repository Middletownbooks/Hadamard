#!/bin/bash
# LP slices at l = 2v-1 for the eleven open v (cyclic Z_l, all multiplier subgroups, r <= RMAX)
R=${1:-30}
for spec in "523:[('C',5),('C',11),('C',19)]" "537:[('C',29),('C',37)]" "571:[('C',7),('C',163)]" "573:[('C',5),('C',229)]" \
            "599:[('E',3,2),('C',7),('C',19)]" "599c:[('C',9),('C',7),('C',19)]" "643:[('C',5),('C',257)]" "647:[('C',3),('C',431)]" \
            "659:[('C',3),('C',439)]" "669:[('C',7),('C',191)]" "719:[('C',3),('C',479)]" "739:[('C',7),('C',211)]"; do
  v=${spec%%:*}; comps=${spec#*:}
  echo "### v=$v comps=$comps" 
  python3 sweep_b1.py "$comps" $R lp
done
