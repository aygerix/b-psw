#!/bin/sh
# E1 production: bbc on the instances in inst/index.tsv, sorted by pool size, lambdas in decreasing order.
# Usage: sh run_bbc.sh BIN LAMLIST NODECAP OUTFILE MINPOOL MAXPOOL WORKER NWORK PERINST_TIMEOUT
#   worker k of n takes the instances with (rank mod n) == k; each instance runs under its own timeout.
# Output: "== D j1 j2 npool nitems nrows bad" header, bbc's lines, "rc <exit code>" (124 = per-instance timeout).
BIN=$1; LAMS=$2; CAP=$3; OUT=$4; MINP=$5; MAXP=$6; W=$7; NW=$8; TL=$9
sort -t "$(printf '\t')" -k4,4n ${INDEX:-inst/index.tsv} | awk -F"\t" -v a=$MINP -v b=$MAXP '$4 >= a && $4 <= b' | \
awk -v w=$W -v n=$NW '(NR - 1) % n == w' | while IFS="$(printf '\t')" read D J1 J2 NP NI NR BAD FN; do
  echo "== $D $J1 $J2 $NP $NI $NR $BAD" >> $OUT
  timeout $TL nice -n 10 $BIN $FN $LAMS $CAP >> $OUT 2>&1
  echo "rc $?" >> $OUT
done
echo "ALLDONE" >> $OUT
