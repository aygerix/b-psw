#!/bin/sh
# run bfx.py after headline worker 1 finished (runs/all_l20_w1.out ends with ALLDONE)
cd ~/b-psw/src/cs/fxchk
while ! tail -1 runs/all_l20_w1.out 2>/dev/null | grep -q ALLDONE; do sleep 15; done
OMP_NUM_THREADS=1 timeout 3600 nice -n 10 $HOME/chowla_venv/bin/python3 bfx.py 8 in/passA_w0.out in/passA_w1.out in/passB_w0.out in/passB_w1.out in/passC_w0.out in/passC_w1.out in/l20_merged.out in/l225_merged.out in/xd_w0.out in/xd_big.out in/l20_D37.out > runs/bfx8.out 2>&1
