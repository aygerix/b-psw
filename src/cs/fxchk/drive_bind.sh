#!/bin/sh
# Per-D binding instances at lambda 0.2 with T = claimed M_ub (reproduce the certificate value).
cd ~/b-psw/src/cs/fxchk
OUT=runs/bind_l20.out; : > $OUT
while read KEY FN T; do
  R=$(timeout 3600 nice -n 10 ./bb2 $FN 0.2 $T | sed 's/best_exact [^ ]*/best_exact X/')
  echo "$KEY T=$T $R" >> $OUT
done < bind_list.txt
echo ALLDONE >> $OUT
