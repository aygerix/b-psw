#!/bin/sh
# setup for the F1 fix validation (run from src/te on the VM)
set -e
mkdir -p valF/hp
cp anchor3 valF/anchor3_old            # pre-fix production binary (built from the committed source)
cp val/anchor3.c valF/anchor3_old.c     # pre-fix source (same md5 as the committed one)
mv anchor3_new.c anchor3.c
gcc -O2 -fopenmp -Wall -Wno-misleading-indentation -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function -o anchor3 anchor3.c -lgmp -lm 2> valF/build_warnings.txt || { echo BUILD FAILED; cat valF/build_warnings.txt; exit 1; }
md5sum anchor3.c anchor3 valF/anchor3_old valF/anchor3_old.c
cd valF
cp ../anchor3 ../anchor3.c .
for f in run_c1.sh run_split.sh run_split2.sh runq.sh cmp_c1.py cmp_c3.py cmp_split.py mkw1.py ref_split.gp q3a.txt q3c.txt q8.txt \
         plants_D-7_X1e6_capA.txt plants_D-7_X1e6_lcm_capA.txt plants_D5_X1e6_capA.txt plants_D5_X1e6_lcm_capA.txt \
         plants_D5_X1e7_S1.txt plants_D5_X1e7_capA.txt plants_D5_X1e7_lcm_capA.txt; do cp ../val/$f .; done
mkdir -p c6ref; cp ../val/c6ref/*.txt c6ref/
python3 mkw1.py anchor3.c anchor3_w1.c
gcc -O2 -fopenmp -o anchor3_w1 anchor3_w1.c -lgmp -lm
# hit-print test variants: final test in bpsw mode replaced by strong Fermat base 2 (so bpsw-mode HIT lines occur)
for v in new old; do
  if [ $v = new ]; then s=anchor3.c; else s=anchor3_old.c; fi
  python3 - "$s" hp/hp_$v.c <<'PY'
import sys
src = open(sys.argv[1]).read()
a = 'bpsw_mpz(&ts->b, ts->n)) record_hit("HIT"'
c = src.count(a)
if c != 3: sys.exit("count %d" % c)
open(sys.argv[2], 'w').write(src.replace(a, 'sprp2_mpz(&ts->b, ts->n)) record_hit("HIT"'))
print("hp ok", sys.argv[2])
PY
  gcc -O2 -fopenmp -o hp/hp_$v hp/hp_$v.c -lgmp -lm
done
ls -la
