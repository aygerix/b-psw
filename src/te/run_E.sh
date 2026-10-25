#!/bin/sh
# Theorem E production driver: for each D, one anchor3 run (one table build) with one pass per split part
# (src/te/splitparts_1e10.txt; S = 1 included explicitly). Output: runs/E_<X>_D<D>.out
# Usage: sh run_E.sh X THREADS D1 [D2 ...]
X=$1; T=$2; shift 2
mkdir -p runs
for D in "$@"; do
  ARGS=$(awk -v d="$D" '$1 == d { printf " --split %s %s %s", $2, $3, $4 }' splitparts_1e10.txt)
  echo "D=$D passes: $(echo $ARGS | grep -o -- '--split' | wc -l)" > runs/E_${X}_D${D}.out
  /usr/bin/time -f "TIME %e s wall %U s user %M KB" ./anchor3 --D $D --X $X --mode bpsw --threads $T $ARGS >> runs/E_${X}_D${D}.out 2>&1
done
