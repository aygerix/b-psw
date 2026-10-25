#!/bin/bash
cd ~/b-psw/src/cs/e3b
PY=$HOME/chowla_venv/bin/python3
if [ "$1" = 02 ]; then
  timeout 3600 nice -n 10 $PY node_model.py @out/s7part_02 10 1000 702 --out out/node7_hits.txt >> out/node7.log 2>/dev/null
  echo DONE02 >> out/node7_done
fi
if [ "$1" = 03 ]; then
  while kill -0 510639 2>/dev/null; do sleep 20; done
  echo DONE01 >> out/node7_done
  timeout 3600 nice -n 10 $PY node_model.py @out/s7part_03 10 1000 703 --out out/node7_hits.txt >> out/node7.log 2>/dev/null
  echo DONE03 >> out/node7_done
fi
