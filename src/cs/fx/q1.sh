#!/bin/sh
# FX queue 1 (one core): lambda = 0.2 certificates, run after the D = -7 (2,1) job (pid $1) exits.
# Each job: bbc_fx, progress every 1e6 nodes on stderr, timeout 3600 s.  Node estimates (Knuth, 600 dives) in the comments.
cd ~/b-psw/src/cs/fx
while kill -0 $1 2>/dev/null; do sleep 15; done
job() { NAME=$1; FN=$2; INCV=$3
  if [ -n "$INCV" ]; then export INC=$INCV; else unset INC; fi
  PROG=1000000 timeout 3600 nice -n 10 ./bbc_fx $FN 0.2 > l20/$NAME.out 2> l20/$NAME.err; echo "rc $?" >> l20/$NAME.out; }
job D-7_1_2 ../e1/inst/inst_D-7_1_2.txt ""               # est 3.0e6 (INC 1.0002)
job D-15_2_0_FF split15/m15_FF.txt 6.8722                # est 3.4e5
job D-15_2_0_FL split15/m15_FL.txt 6.8722                # est 3.4e5
job D-15_2_0_LF_F split15/m15_LF_F.txt 6.8722            # est 3.4e6
job D-15_2_0_LF_L split15/m15_LF_L.txt 6.8722            # est 4.8e6
job D-15_2_0_LL_L split15/m15_LL_L.txt 6.8722            # est 7.3e6
echo Q1DONE > l20/q1.done
