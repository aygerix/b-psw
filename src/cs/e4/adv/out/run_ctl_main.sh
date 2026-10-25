#!/bin/bash
# controls for the main sets (1 core; queue A is the other core)
cd ~/b-psw/src/cs/e4/adv
P=${PY:-python3}
timeout 1800 nice -n 10 $P adv.py gen 3 20261025 ctl_main_in.gp mF.out mL.out > gen_main.log 2>&1
printf 'default(parisizemax, 10^9);\nTL=60; IN="ctl_main_in.gp"; OUT="cM.out";\nread("adv_fac.gp");\n' > drv_cM.gp
rm -f cM.out
/usr/bin/time -v -o cM.time timeout 3600 nice -n 10 gp -q drv_cM.gp > cM.log 2>&1
echo CTLDONE > ctl_main_done.flag
