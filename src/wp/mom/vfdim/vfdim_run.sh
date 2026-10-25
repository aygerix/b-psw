#!/bin/bash
# verifier rerun of m9_kcert.gp for the four cases (2 at a time)
cd ~/b-psw/src/wp/mom
ulimit -v 3000000
run() { echo "SOLFILE=\"m9_sol_$1.gp\"; read(\"m9_kcert.gp\")" | timeout 3600 nice -n 10 gp -q -s 400M > vfdim/kcert_$1.out 2>&1; echo "done $1 rc $?" >> vfdim/status; }
run 9_0 & run 9_1 & wait
run 11_0 & run 11_1 & wait
echo ALLDONE >> vfdim/status
