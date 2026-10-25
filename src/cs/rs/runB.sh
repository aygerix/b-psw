#!/bin/bash
# Experiment B sweep (lattice, real vs controls). Output: expB_*.out (SUMMARY lines)
cd ~/b-psw/src/cs/rs
PY=~/cs_venv/bin/python3; I=../../../data/bpsw_instance_T1500_1e17.json
run() { nice -n 10 timeout 1700 $PY lat.py $I "$@" 2>&1 | grep SUMMARY || echo "TIMEOUT/FAIL $*"; }
( for B in 30 40 50 60; do for c in real unif supp hgt; do run $B +4 10 $B rnd $c 20 enum; done; done ) > expB_scale_enum.out
( for B in 70 80 100 120; do for c in real unif supp hgt; do run $B +4 8 $B rnd $c 20; done; done ) > expB_scale_bkz.out
( for B in 40 60 80; do for c in real hgt unif; do run $B +4 8 $((B+7)) low $c 20; done; done ) > expB_low.out
( for B in 40 60; do for f in 1.5 2 3; do for c in real unif; do run $B $f 6 $((B+11)) rnd $c 20; done; done; done ) > expB_dense.out
echo done > expB.done
