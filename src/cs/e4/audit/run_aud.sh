#!/bin/bash
# Audit rerun of E4 (instrumented fic_aud.gp; same parameters and seeds as the original drivers). Two queues.
cd ~/b-psw/src/cs/e4/audit
q() { for j in "$@"; do /usr/bin/time -v -o $j.time timeout 3600 nice -n 10 gp -q drv_$j.gp > $j.log 2>&1; done; }
q aF apilR &
q aL arep1 arep2 &
wait
echo ALLDONE > done.flag
