# MILP (independent) vs bbc on sub-instances; then MILP on full instances.
PY=${PY:-python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for spec in "-15 2 0 43 0.1,0.15,0.2" "-7 1 2 43 0.1,0.15,0.2,0.25" "17 2 0 43 0.1,0.15,0.2"; do
  set -- $spec
  echo "== MILP D=$1 ($2,$3) T0=$4"; timeout 600 nice -n 10 $PY aud_milp.py $1 $2 $3 $4 $5 120 ../c3_1e10.txt ../chm200_pool.txt
done
