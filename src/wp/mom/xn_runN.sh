#!/bin/bash
# xn line: order-N feasibility scan, k = $1, lin = $2, N from $3 to $4, $5 starts each. Output xn_feas_k{k}_l{lin}_N{N}.out
cd ~/b-psw/src/wp/mom
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=~/chowla_venv/bin/python
for N in $(seq $3 $4); do
  timeout 3500 nice -n 10 $PY xn_feas.py $1 $2 $N $5 11 > xn_feas_k$1_l$2_N$N.out 2>&1
done
