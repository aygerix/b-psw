#!/bin/bash
# audC_lin_batch.sh -- runs of audC_lin.py (2 at a time, nice 10, timeout 3000 s each)
cd ~/b-psw/src/wp/mom; P=~/chowla_venv/bin/python3; T="timeout 3000 nice -n 10"
run() { $T $P -W ignore audC_lin.py "$@" > audC_lin_$(echo "$@" | tr ' ' '_').out 2>&1; }
( run rand 3 21 lin 100000; run opt 3 22 lin 60; run rand 4 23 lin 50000; run opt 4 24 lin 40; run opt 5 25 lin 30 ) &
( run opt 2 31 lin 200; run opt 3 32 nolin 60; run genopt 3 33 nolin 40; run genopt 3 34 lin 40; run genopt 2 35 lin 60 ) &
wait
echo DONE > audC_lin_batch.done
