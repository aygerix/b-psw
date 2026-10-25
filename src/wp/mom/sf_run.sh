#!/bin/bash
# sf_run.sh LIST OUT : run sf_lm.py (deflated, NFEV 100 per unknown) on each line "MODE DEGS NNEG STARTS SEED" of LIST
P=$HOME/chowla_venv/bin/python3
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
while read -r mode degs nneg starts seed; do
  [ -z "$mode" ] && continue
  timeout 1500 nice -n 10 $P sf_lm.py $mode $degs $nneg $starts $seed defl 100 >> "$2" 2>&1 || echo "TIMEOUT/ERR $mode $degs $nneg" >> "$2"
done < "$1"
echo DONE >> "$2"
