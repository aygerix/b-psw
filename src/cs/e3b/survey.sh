#!/bin/bash
# src/cs/e3b/survey.sh -- E3b numerical survey driver (k = 5): stream A = Sum m = 7, stream B = Sum m = 6 + symmetric.
# Each pattern: strat_lp.py (NB batches x 1000 random starts) -> out/sv_<tag>.txt, then refine.py -> out/sv_<tag>.ref
cd ~/b-psw/src/cs/e3b
PY=$HOME/chowla_venv/bin/python3
run() {  # tag pattern nb [extra]
  timeout 3600 nice -n 10 $PY strat_lp.py 5 $2 $3 1000 $RANDOM --out out/sv_$1.txt $4 > out/sv_$1.log 2>&1
  timeout 3600 nice -n 10 $PY refine.py 5 out/sv_$1.txt 300 > out/sv_$1.ref 2>&1
}
if [ "$1" = A ]; then
  for m in 2,3,2,0,0 2,3,1,1,0 2,2,2,1,0 2,2,1,1,1 1,3,3,0,0 1,3,2,1,0 1,3,1,1,1 1,2,2,2,0 1,2,2,1,1 0,3,3,1,0 0,3,2,2,0 0,3,2,1,1 0,2,2,2,1; do
    run s7_${m//,/} $m 30
  done
fi
if [ "$1" = B ]; then
  for m in 2,3,1,0,0 2,2,2,0,0 2,2,1,1,0 2,1,1,1,1 1,3,2,0,0 1,3,1,1,0 1,2,2,1,0 1,2,1,1,1 0,3,3,0,0 0,3,2,1,0 0,3,1,1,1 0,2,2,2,0 0,2,2,1,1; do
    run s6_${m//,/} $m 30
  done
  for m in 2,2,2,0,0 2,1,1,1,1 0,3,3,0,0 0,2,2,1,1 2,3,3,0,0 2,2,2,1,1 0,3,3,1,1 0,2,2,2,2; do
    run sym_${m//,/} $m 20 --sym
  done
fi
echo DONE $1
