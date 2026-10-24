#!/bin/sh
# Locate the Pomerance crossover T*: pools to 1e15 for T = 1100..1400, then the greedy frontier.
cd ~/b-psw/src/c3 || exit 1
for T in 1100 1200 1300 1400; do
  nice -n 15 ./pomgen2 $T 15 pom_${T}_15.txt > pomgen_${T}_15.out 2>&1
  (cd ../pom && nice -n 15 python3 profile2.py ../c3/pom_${T}_15.txt 2000 > prof2_${T}_15.out 2>&1)
done
echo done > tstar.done
