# Falsification batch: local search on full instances; compare "best N - lam bits" with the certified M_ub.
PY=${PY:-python3}; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for spec in "$@"; do
  set -- $spec
  timeout 1500 nice -n 10 $PY aud_ls.py $1 $2 $3 0.2,0.25,0.3 300 11 -1 ../c3_1e10.txt ../chm200_pool.txt
done
