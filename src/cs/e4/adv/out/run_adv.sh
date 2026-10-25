#!/bin/bash
# adversarial E4 re-run: two queues (<= 2 cores), every job under timeout 3600 and nice 10
cd ~/b-psw/src/cs/e4/adv
q() { for j in "$@"; do rm -f $j.out; /usr/bin/time -v -o $j.time timeout 3600 nice -n 10 gp -q drv_$j.gp > $j.log 2>&1; done; }
q mF rA &
q mL rB &
wait
echo ALLDONE > adv_done.flag
