#!/bin/bash
cd ~/b-psw/src/cs/e3b
for p in 00 01 02 03; do
  timeout 3600 nice -n 10 $HOME/chowla_venv/bin/python3 node_model.py @out/s7part_$p 10 1000 7$p --out out/node7_hits.txt >> out/node7.log 2>/dev/null
done
echo DONE >> out/node7.log
