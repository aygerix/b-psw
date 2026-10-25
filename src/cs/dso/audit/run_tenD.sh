#!/bin/sh
# Ten-D audit at X = 1e4 (single core, sequential): reproduce the census with the original script, recompute
# orders independently, check completeness with aud_census (Mobius/Mertens + recurrence), run the C exact search.
cd ~/b-psw/src/cs/dso/audit
X=10000
for D in 5 -7 -11 13 -15 17 -19 -23 29 -31; do
  Q=$(( (1 - D) / 4 ))
  tag=$(echo $D | sed 's/-/m/')
  rm -f rep_D${tag}_X$X.txt rows_D${tag}_X$X.txt
  echo "== D = $D (Q = $Q)"
  echo "X=$X; DD=$D; OUT=\"rep_D${tag}_X$X.txt\"; read(\"../dso_census_D.gp\")" | timeout 600 nice -n 10 gp -q -D parisizemax=2000000000 2>&1 | grep -v Warning | tail -2
  orig=../dsoD_${tag}_X$X.txt
  if cmp -s rep_D${tag}_X$X.txt $orig; then echo "census rerun identical to $orig"; else echo "census rerun DIFFERS from $orig"; fi
  echo "IN=\"$orig\"; DD=$D; X=$X; OUTP=\"rows_D${tag}_X$X.txt\"; read(\"aud_orders.gp\")" | timeout 600 nice -n 10 gp -q -D parisizemax=2000000000 2>&1 | tail -2
  cut -d' ' -f1 rows_D${tag}_X$X.txt > prD_${tag}.txt
  timeout 600 nice -n 10 ./aud_census $X $Q prD_${tag}.txt | grep -E "listed|COFACTOR|ERROR"
  timeout 900 nice -n 10 ./aud_search rows_D${tag}_X$X.txt $D 5 0 1 | grep -E "usable|k = |SOL"
done
