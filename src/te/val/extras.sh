#!/bin/sh
# extra validation jobs (single process each)
mkdir -p c5
[ -f c5/inert.out ] || timeout 3600 gp -q inert_count.gp > c5/inert.out 2>&1
[ -f c5/pair_1e8.out ] || /usr/bin/time -f "TIME %e s wall %U s user %M KB" timeout 3600 ./harness_pair 90000001 100000000 6 > c5/pair_1e8.out 2>&1
[ -f c5/pair_1e9.out ] || /usr/bin/time -f "TIME %e s wall %U s user %M KB" timeout 3600 ./harness_pair 960000001 1000000000 6 > c5/pair_1e9.out 2>&1
[ -f c5/pair_4e9.out ] || /usr/bin/time -f "TIME %e s wall %U s user %M KB" timeout 3600 ./harness_pair 4260000001 4300000000 6 > c5/pair_4e9.out 2>&1
[ -f c5/pair_1e10.out ] || /usr/bin/time -f "TIME %e s wall %U s user %M KB" timeout 3600 ./harness_pair 9960000001 10000000000 6 > c5/pair_1e10.out 2>&1
touch extras.done
