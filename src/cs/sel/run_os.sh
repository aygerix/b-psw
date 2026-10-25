#!/bin/bash
# Selector line: one-sided index-straggler runs (sequential, one core). Output: os_<mode>_<T>_<lx>.{out,txt}
cd ~/b-psw/src/cs/sel
for job in "1 300 17" "2 300 17" "1 500 17" "2 500 17" "1 1500 15" "2 1500 15" "1 1000 16" "2 1000 16"; do
  set -- $job
  /usr/bin/time -f "%e s" timeout 1800 nice -n 10 ./onesided $1 $2 $3 os_$1_$2_$3.txt > os_$1_$2_$3.out 2>&1
done
echo done > run_os.done
