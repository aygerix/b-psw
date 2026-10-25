#!/bin/sh
# FX queue 4: lambda = 0.225 for D = -15 (2,0) on the split sub-instances (Knuth estimate of the whole tree 1.3e7 nodes,
# INC 0.2 or 1.0 alike; each sub-tree is a part of it).  Waits for WAITFILE, then runs the listed sub-instances.
# Usage: sh q4.sh WAITFILE sub1 sub2 ...      (sub = FF, FL, LF_F, LF_L, LL_F_F, LL_F_L, LL_L)
cd ~/b-psw/src/cs/fx
W=$1; shift
while [ ! -f $W ]; do sleep 30; done
mkdir -p l225
for s in "$@"; do
  INC=0.2 PROG=1000000 timeout 3600 nice -n 10 ./bbc_fx split15/m15_$s.txt 0.225 > l225/D-15_2_0_$s.out 2> l225/D-15_2_0_$s.err
  echo "rc $?" >> l225/D-15_2_0_$s.out
done
touch l225/done_$(echo "$@" | tr ' ' '_')
