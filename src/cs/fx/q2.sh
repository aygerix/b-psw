#!/bin/sh
# FX queue 2 (one core).  (a) Knuth estimates for every remaining lambda = 0.2 job (HARD RULE), logged to l20/knuth_q2.log;
# a job runs only if its estimate is <= 9e6 nodes.  (b) D = 37 (2,0) lambda = 0.2.  (c) local search T0 = 300 vs T0 = 200
# (same settings as the E1 audit: 300 restarts, seed 11, every single-prime climb).  (d) remaining lambda = 0.2 jobs.
cd ~/b-psw/src/cs/fx
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=${PY:-python3}
est() { KNUTH=600 INC=$3 timeout 600 nice -n 10 ./bbc_fx $2 0.2 | sed "s|^|$1 |" >> l20/knuth_q2.log; }
ok() { awk -v n=$1 '$1 == n {print ($10 <= 9e6) ? 1 : 0}' l20/knuth_q2.log | tail -1; }
job() { NAME=$1; FN=$2; INCV=$3
  if [ "$(ok $NAME)" != "1" ]; then echo "$NAME not run (estimate)" >> l20/knuth_q2.log; return; fi
  if [ -n "$INCV" ]; then export INC=$INCV; else unset INC; fi
  PROG=1000000 timeout 3600 nice -n 10 ./bbc_fx $FN 0.2 > l20/$NAME.out 2> l20/$NAME.err; echo "rc $?" >> l20/$NAME.out; unset INC; }
est D37_2_0 instxD/inst_D37_2_0.txt 0
est D-11_2_1 ../e1/inst/inst_D-11_2_1.txt 0
est D-11_2_1i ../e1/inst/inst_D-11_2_1.txt 0.2
est D-11_2_0_F split11/m11_20_F.txt 0.2832
est D-11_2_0_L split11/m11_20_L.txt 0.2832
est D-15_2_0_LL_F_F split15/m15_LL_F_F.txt 6.8722
est D-15_2_0_LL_F_L split15/m15_LL_F_L.txt 6.8722
job D37_2_0 instxD/inst_D37_2_0.txt ""
for spec in "300 -15 2 0 c3_T300_1e10.txt" "300 -7 1 2 c3_T300_1e10.txt" "200 -7 1 2 ../e1/c3_1e10.txt" "300 5 2 0 c3_T300_1e10.txt" "200 5 2 0 ../e1/c3_1e10.txt"; do
  set -- $spec
  T0=$1 DUMP=214,275 timeout 3000 nice -n 10 $PY ls300.py $2 $3 $4 0.2,0.25,0.3 300 11 -1 $5 ../e1/chm200_pool.txt > ls/ls_T$1_D$2_$3_$4.out 2>&1
done
job D-11_2_1 ../e1/inst/inst_D-11_2_1.txt ""
job D-15_2_0_LL_F_F split15/m15_LL_F_F.txt 6.8722
job D-15_2_0_LL_F_L split15/m15_LL_F_L.txt 6.8722
echo Q2DONE > l20/q2.done
