#!/bin/bash
# src/cs/e3b/advchk/rerun_all.sh -- regenerate the saved outputs of the E3b adversarial check (cheap steps; <= 2 cores).
# Stream A (this script, sequential): exact/symbolic checks, recall tests, unions, field tests.
# Stream B (run_B below, launched in parallel): symmetric double-pencil enumeration and k = 3 calibration.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=$HOME/chowla_venv/bin/python3
D=~/b-psw/src/cs/e3b/advchk
cd $D
if [ "$1" = B ]; then
  timeout 1500 nice -n 10 $PY symdp_nodes.py "1|2|3,12|3,1|2|2,12|2,2|3|3,2|2|2" 30 1000 11 --out out/symdp6_rerun.txt > out/symdp6.log 2>/dev/null
  cd ..
  for m in 0,2,2 2,1,1 2,2,1 1,2,2; do
    timeout 600 nice -n 10 $PY strat_lp.py 3 $m 3 1000 7 2>/dev/null | grep -v batch
  done > $D/out/cal3.log
  echo DONE > $D/out/rerun_B_done
  exit 0
fi
timeout 300 nice -n 10 gp -q -s 2G lemmas.gp < /dev/null > out/lemmas.out 2>&1
timeout 120 nice -n 10 $PY genus.py > out/genus.out 2>&1
timeout 300 nice -n 10 gp -q -s 4G vertex_sym.gp < /dev/null > out/vertex_sym.out 2>&1
(cd .. && timeout 600 nice -n 10 $PY advchk/recall_chk.py out/sym_22200.txt 6 > $D/out/recall_chk.out 2>/dev/null)
timeout 600 nice -n 10 $PY sq_selftest.py > out/sq_selftest.out 2>/dev/null
timeout 600 nice -n 10 $PY dpen_rat.py 1.988329 -4.369707 7.059911 0.156827 0.076613j 0.666298 > out/dpen_rat.out 2>&1
timeout 900 nice -n 10 $PY fsp_union.py 1 out/fsptu_1.txt,out/fsptu_1b.txt ../out/fspn_1_s21.txt,../out/fspn_1_s22.txt,../out/fspn_1.txt > out/fsp_union_1.out 2>/dev/null
timeout 900 nice -n 10 $PY fsp_union.py 2 out/fsptu_2.txt,out/fsptu_2b.txt ../out/fspn_2_s21.txt,../out/fspn_2_s22.txt,../out/fspn_2.txt > out/fsp_union_2.out 2>/dev/null
timeout 900 nice -n 10 $PY symdp_rat.py out/symdp6.txt 60 > out/symdp_rat.out 2>/dev/null
timeout 1200 nice -n 10 $PY fsp_rat.py 1 out/fsptu_1.txt 60 > out/fsp_rat_1.out 2>/dev/null
timeout 1200 nice -n 10 $PY fsp_rat.py 2 out/fsptu_2b.txt 60 > out/fsp_rat_2b.out 2>/dev/null
echo DONE > out/rerun_A_done
