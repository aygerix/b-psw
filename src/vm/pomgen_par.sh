#!/bin/sh
# Parallel Pomerance-pool generation: sh pomgen_par.sh T log10Xmax W  -> c3/pom_T_X.txt (merged) and per-worker logs.
T=$1; X=$2; W=$3
cd ~/b-psw/src/c3 || exit 1
i=0
while [ $i -lt $W ]; do
  nice -n 2 ./pomgen3 $T $X pom_${T}_${X}_w$i.txt $W $i > pomgen3_${T}_${X}_w$i.out 2>&1 &
  i=$((i+1))
done
wait
sort -n -u pom_${T}_${X}_w*.txt > pom_${T}_${X}.txt
wc -l pom_${T}_${X}.txt > pom_${T}_${X}.done
