#!/bin/sh
# Theorem E production driver: for each D, one anchor3 run (one table build) with one pass per split part
# (src/te/splitparts_1e10.txt; S = 1 included explicitly). Output: runs/E_<X>_D<D>.out
# After each D the run is CHECKED: exit code 0, exactly one PASS line per split part, and the hit total.
# A missing PASS line (OOM kill, die()) must never be read as "no hits": such a D is reported as INCOMPLETE.
# Usage: sh run_E.sh X THREADS D1 [D2 ...]
X=$1; T=$2; shift 2
mkdir -p runs
for D in "$@"; do
  ARGS=$(awk -v d="$D" '$1 == d { printf " --split %s %s %s", $2, $3, $4 }' splitparts_1e10.txt)
  NEXP=$(awk -v d="$D" '$1 == d' splitparts_1e10.txt | wc -l)
  OUT=runs/E_${X}_D${D}.out
  echo "D=$D passes expected: $NEXP" > $OUT
  /usr/bin/time -f "TIME %e s wall %U s user %M KB" ./anchor3 --D $D --X $X --mode bpsw --threads $T --cpair 0.1 $ARGS >> $OUT 2>&1
  RC=$?
  NPASS=$(grep -c "^PASS" $OUT)
  HITS=$(grep "^PASS" $OUT | sed -n 's/.* hits=\([0-9]*\).*/\1/p' | awk '{s += $1} END {print s + 0}')
  NHITLINES=$(grep -c "HIT" $OUT)
  if [ "$RC" -eq 0 ] && [ "$NPASS" -eq "$NEXP" ]; then ST=COMPLETE; else ST=INCOMPLETE; fi
  echo "SUMMARY D=$D X=$X status=$ST rc=$RC passes=$NPASS/$NEXP hits=$HITS hitlines=$NHITLINES" | tee -a $OUT >> runs/E_${X}_summary.txt
done
