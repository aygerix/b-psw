#!/bin/bash
# audC_pc.sh -- positive-control search: general-position configurations with no two-sided single-point test (J > 0).
cd ~/b-psw/src/wp/mom; P=~/chowla_venv/bin/python3; T="timeout 3000 nice -n 10"
run() { $T $P -W ignore audC_lin.py "$@" > audC_pc_$(echo "$@" | tr ' ' '_').out 2>&1; }
( run genopt 4 41 nolin 30; run genopt 6 43 nolin 15; run genopt 8 45 nolin 8 ) &
( run genopt 5 42 lin 20; run genopt 7 44 lin 10; run genopt 9 46 nolin 6 ) &
wait; echo DONE > audC_pc.done
