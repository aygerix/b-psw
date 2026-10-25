#!/bin/sh
# Line OC: full generation, one process, sequential. Outputs in out/full_<fam>.out
cd ~/b-psw/src/cs/oc
for spec in "mer 4000" "wag 5000" "fib 5000" "luc 5000" "gm 4000" "cyc 2000"; do
  set -- $spec
  /usr/bin/time -f "$1 %e s" sh -c "echo 'FAM=\"$1\"; CAP=$2; read(\"oc_gen.gp\")' | timeout 3600 nice -n 10 gp -q -s 800000000 > out/full_$1.out" 2>> out/full_time.txt
done
for spec in "rie" "pro"; do
  /usr/bin/time -f "$spec %e s" sh -c "echo 'FAM=\"$spec\"; CMAX=127; CAP=600; read(\"oc_gen.gp\")' | timeout 3600 nice -n 10 gp -q -s 800000000 > out/full_$spec.out" 2>> out/full_time.txt
done
echo DONE >> out/full_time.txt
