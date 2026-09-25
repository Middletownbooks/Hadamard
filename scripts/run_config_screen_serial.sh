#!/bin/bash
# option-2 configuration screen for v=515 (h=85): calibration u=53 first, then order-8, order-16, order-4 multipliers at p=5
cd /home/user/Hadamard
H=85; P=${P:-5}; NA=${NA:-15}
CFGS=${CFGS:-"53 257 59 161 19 43 121 127 7 11 23 29 31 37 79 173 13 47 89 191 203 217 259 353 361 67 137"}
echo $CFGS | tr ' ' '\n' | xargs -P1 -I{} bash -c "python3 -u scripts/bge3_config_screen.py $H $P {} $NA > logs/config_screen/h${H}_p${P}_u{}.log 2>&1"
touch logs/config_screen/done_p${P}
