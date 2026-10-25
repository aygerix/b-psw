#!/bin/sh
# FX queue 3: after queue 1 finishes (l20/q1.done), the D = -11 (2,0) lambda = 0.2 halves (split by row 5).
cd ~/b-psw/src/cs/fx
while [ ! -f l20/q1.done ]; do sleep 30; done
for s in F L; do
  E=$(awk -v n=D-11_2_0_$s '$1 == n {print ($10 <= 9e6) ? 1 : 0}' l20/knuth_q2.log | tail -1)
  if [ "$E" != "1" ]; then echo "D-11_2_0_$s not run (estimate)" >> l20/knuth_q2.log; continue; fi
  INC=0.2832 PROG=1000000 timeout 3600 nice -n 10 ./bbc_fx split11/m11_20_$s.txt 0.2 > l20/D-11_2_0_$s.out 2> l20/D-11_2_0_$s.err
  echo "rc $?" >> l20/D-11_2_0_$s.out
done
echo Q3DONE > l20/q3.done
