#!/bin/sh
# Batch 3 (single core): heuristic reruns with a binary rebuilt from dso_heur.c, plus the mixed-cell term.
cd ~/b-psw/src/cs/dso/audit
PY=$HOME/chowla_venv/bin/python3
mkdir -p heurrun; cd heurrun; ln -sf ../dso_heur_rebuilt dso_heur
echo "== band (2e3,4e3]"; /usr/bin/time -f "%e s %M KB" timeout 900 nice -n 10 ./dso_heur 4000 /dev/null 16 n 2000 2>&1 | grep -E "TOTAL|BAND|SIGN| s "
echo "== mixed band (2e3,4e3]"; /usr/bin/time -f "%e s" timeout 900 nice -n 10 ../aud_heur_mixed 4000 /dev/null 16 n 2000 stop 2>&1 | grep -E "MIXED| s"
echo "== mixed X=1e4 (primes > 1e4)"; /usr/bin/time -f "%e s" timeout 1800 nice -n 10 ../aud_heur_mixed 10000 /dev/null 16 n 10000 stop 2>&1 | grep -E "MIXED| s"
echo "== heur_cmp X=1e4"; /usr/bin/time -f "%e s %M KB" timeout 1800 nice -n 10 $PY ../../dso_heur_cmp.py ../../dso_X10000.txt 10000 > ../rerun_heur_X10000.out 2>&1; grep -E "TOTAL|SIGN|observed p|expected small|cells by" ../rerun_heur_X10000.out
echo BATCH3-DONE
