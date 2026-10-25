#!/bin/sh
# FX driver: bbc_fx on the instances of an index file, lambdas in decreasing order, with a Knuth gate (HARD RULE: estimate
# first).  For every instance with pool >= GATE and every lambda, KNUTH=200 dives estimate the node count (INC = 0, an
# overestimate of the DFS tree); if est > MAXN the instance is marked "SKIP" for that lambda and smaller ones.
# Otherwise the exact run gets NODECAP = 3*MAXN and a per-run timeout TL (rc 124 = timeout).
# Usage: sh run_fx.sh INDEX LAMS(comma) OUT WORKER NWORK TL MAXN GATE
IDX=$1; LAMS=$2; OUT=$3; W=$4; NW=$5; TL=$6; MAXN=$7; GATE=$8
BIN=$(dirname $0)/bbc_fx
sort -t "$(printf '\t')" -k4,4n $IDX | awk -v w=$W -v n=$NW '(NR - 1) % n == w' | while IFS="$(printf '\t')" read D J1 J2 NP NI NR BAD FN; do
  echo "== $D $J1 $J2 $NP $NI $NR $BAD" >> $OUT
  for L in $(echo $LAMS | tr ',' ' '); do
    if [ $NP -ge $GATE ]; then
      E=$(KNUTH=200 timeout 600 nice -n 10 $BIN $FN $L | awk '/^knuth/{print $9}')
      echo "knuth $L $E" >> $OUT
      if [ -z "$E" ] || [ $(echo "$E > $MAXN" | awk '{print ($1 > $3) ? 1 : 0}') -eq 1 ]; then echo "SKIP $L est $E" >> $OUT; break; fi
    fi
    timeout $TL nice -n 10 $BIN $FN $L $((3 * MAXN)) >> $OUT 2>&1
    RC=$?; echo "rc $RC" >> $OUT
    if [ $RC -ne 0 ]; then break; fi
    if tail -3 $OUT | grep -q NODECAP; then break; fi
  done
done
echo "ALLDONE" >> $OUT
