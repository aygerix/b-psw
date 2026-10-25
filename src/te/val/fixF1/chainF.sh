#!/bin/sh
# F1 validation chain (run from src/te/valF on the VM); one 4-thread job at a time
T="/usr/bin/time -f"
# 0. self-test of the fixed binary
$T "TIME %e s wall %U s user %M KB" timeout 1800 ./anchor3 --selftest > selftest.out 2>&1; echo "rc=$?" >> selftest.out
# 1. hit-print test: bpsw-mode HIT lines with and without --quiet (sprp2 final-test variants, old and new source)
for v in new old; do for XC in 20000:auto 20000:0 20000:-1 200000:auto 200000:0; do X=${XC%%:*}; cap=${XC#*:}; for q in "" --quiet; do
  if [ "$cap" = auto ]; then capopt=""; else capopt="--cap $cap"; fi
  out=hp/${v}_X${X}_cap${cap}${q}.out
  timeout 1800 ./hp/hp_$v --D 5 --X $X --mode bpsw $capopt $q --threads 4 > $out 2>&1; echo "rc=$?" >> $out
done; done; done
# 2. calibration D = 5, X = 1e7, bpsw (as c4/bpsw_d5_1e7.out), also with --quiet
mkdir -p c4
$T "TIME %e s wall %U s user %M KB" timeout 3600 ./anchor3 --D 5 --X 1e7 --threads 4 > c4/bpsw_d5_1e7.out 2>&1; echo "rc=$?" >> c4/bpsw_d5_1e7.out
$T "TIME %e s wall %U s user %M KB" timeout 3600 ./anchor3 --D 5 --X 1e7 --threads 4 --quiet > c4/bpsw_d5_1e7_quiet.out 2>&1; echo "rc=$?" >> c4/bpsw_d5_1e7_quiet.out
touch stage2.done
# 3. C1/C2 one-torus brute-force equality and bpsw/weak vs references
./run_c1.sh > c1_driver.log 2>&1
# 4. split-part controls (C1 S > 1)
./run_split.sh > split.log 2>&1
./run_split2.sh > split2.log 2>&1
touch stage4.done
# 5. C3 planted controls (default cap, cap 0, cap -1)
mkdir -p c3
./runq.sh q3a.txt
./runq.sh q3c.txt
./runq.sh q8.txt
touch chainF.done
