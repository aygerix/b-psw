#!/bin/bash
cd ~/b-psw/src/cs/e3b
PY=$HOME/chowla_venv/bin/python3
for m in 2,1,1,1,1 2,2,2,0,0 0,3,3,0,0 0,2,2,1,1 2,3,3,0,0 2,2,2,1,1 0,3,3,1,1 0,2,2,2,2 2,0,0,3,3; do
  timeout 3000 nice -n 10 $PY strat_lp.py 5 $m 20 1000 $RANDOM --sym --out out/sym_${m//,/}.txt > out/sym_${m//,/}.log 2>&1
done
echo DONE > out/sym_done
