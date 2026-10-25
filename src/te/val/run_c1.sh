#!/bin/sh
# C1/C2 driver: engine control modes vs brute-force reference (run from the val directory)
mkdir -p c1 c1w
for D in 5 -7 13; do for X in 2000 6000 20000; do
  for mode in fermat lucas bpsw weak; do for cap in 0 -1 auto 40; do
    out=c1/D${D}_X${X}_${mode}_cap${cap}.out
    if [ "$cap" = auto ]; then capopt=""; else capopt="--cap $cap"; fi
    /usr/bin/time -f "TIME %e s wall %U s user %M KB" timeout 7200 ./anchor3 --D $D --X $X --mode $mode $capopt --threads 4 > $out 2>&1; echo "rc=$?" >> $out
  done; done
  for mode in fermat lucas; do for cap in 0 -1 auto; do
    out=c1w/D${D}_X${X}_${mode}_cap${cap}.out
    if [ "$cap" = auto ]; then capopt=""; else capopt="--cap $cap"; fi
    /usr/bin/time -f "TIME %e s wall %U s user %M KB" timeout 7200 ./anchor3_w1 --D $D --X $X --mode $mode $capopt --threads 4 > $out 2>&1; echo "rc=$?" >> $out
  done; done
done; done
echo ALLDONE > c1/DONE
