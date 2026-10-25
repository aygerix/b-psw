#!/bin/bash
# ns line job queue: runs the job list given as $1 (one job per line: ENV... -- args to ns_run3.py), sequentially,
# each under timeout 3600 and nice 10, single-threaded BLAS.  Usage: nohup ./ns_queue.sh jobs_A.txt > qA.log 2>&1 &
cd ~/b-psw/src/wp/mom
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=${NS_PY:-python3}   # on the VM set NS_PY to the venv python that has numpy/scipy/mpmath (system python3 lacks numpy)
while IFS= read -r line; do
  [ -z "$line" ] && continue
  case "$line" in \#*) continue;; esac
  name=$(echo "$line" | awk '{print $1}')
  envs=$(echo "$line" | awk -F' -- ' '{print $1}' | cut -d' ' -f2-)
  args=$(echo "$line" | awk -F' -- ' '{print $2}')
  echo "$(date -u +%H:%M:%S) start $name: env [$envs] args [$args]"
  env $envs timeout 3600 nice -n 10 $PY ns_run3.py $args > ns_out/$name.log 2>&1
  echo "$(date -u +%H:%M:%S) done $name (exit $?)"
done < "$1"
