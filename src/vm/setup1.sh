#!/bin/sh
# Session 5 VM setup: per-D copies of the L1 search, then the job queue (at most 7 at once, niced).
cd ~/b-psw/src/lp || exit 1
for D in 41 53 61 73 89 97 101; do
  sed -e "s/DLIST=\[[^]]*\]/DLIST=[$D]/" lp_1e6c.gp > lp_1e6_D$D.gp
done
grep -h "DLIST=" lp_1e6_D*.gp | cut -c1-60
cd ~/b-psw/src/vm || exit 1
nohup xargs -d '\n' -P 7 -I CMD sh -c 'CMD' < jobs1.txt > runq1.log 2>&1 &
echo "queue started, pid $!"
