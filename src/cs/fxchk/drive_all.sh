#!/bin/sh
# Headline check: every certificate instance except D=-15 (2,0), lambda = 0.2, decision threshold T = 1.719
# (0.2 b + 1.719 <= 0.2249 (b + 1) for all b >= 60).  Worker W of 2; instances sorted by pool size, round-robin.
cd ~/b-psw/src/cs/fxchk
W=$1; OUT=runs/all_l20_w$W.out; : > $OUT
( awk -F'\t' '{print $4"\t../e1/"$8"\t"$1" "$2" "$3}' ../e1/out_index_copy.tsv; awk -F'\t' '{print $4"\t../fx/"$8"\t"$1" "$2" "$3}' ../fx/instxD/index.tsv ) \
 | grep -v "	-15 2 0$" | sort -t'	' -k1,1n | awk -v w=$W 'NR % 2 == w' | while IFS='	' read N FN KEY; do
  R=$(timeout 3600 nice -n 10 ./bb2 $FN 0.2 1.719 | sed 's/best_exact [^ ]*/best_exact X/')
  echo "$KEY $N rc=$? $R" >> $OUT
done
echo ALLDONE >> $OUT
