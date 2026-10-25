# bbc (the line's certificate code) vs the independent MILP on T0-restricted sub-instances, several lambdas.
# Usage: sh cmp_sub.sh T0 LAMS "D j1 j2" ["D j1 j2" ...]
PY=${PY:-python3}; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
T=$1; L=$2; shift 2
for spec in "$@"; do
  set -- $spec
  echo "#### D=$1 ($2,$3) T0=$T"
  timeout 900 nice -n 10 ./bbc_audit instT$T/inst_D$1_$2_$3.txt $L | grep "^lam"
  timeout 2400 nice -n 10 $PY aud_milp.py $1 $2 $3 $T $L 300 ../c3_1e10.txt ../chm200_pool.txt
done
