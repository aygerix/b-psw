#!/bin/sh
# F1 hit-print test, second form (run from src/te/valF after chainF.done). Variant hp2: bpsw mode with the one-torus
# modulus M = ord_p(2) and final test strong Fermat base 2, so that bpsw-mode HIT lines (branch A and C) occur.
while [ ! -f chainF.done ]; do sleep 30; done
mkdir -p hp2
for v in new old; do
  if [ $v = new ]; then s=anchor3.c; else s=anchor3_old.c; fi
  python3 - "$s" hp2/hp2_$v.c <<'PY'
import sys
src = open(sys.argv[1]).read()
for a, b, k in (('bpsw_mpz(&ts->b, ts->n)) record_hit("HIT"', 'sprp2_mpz(&ts->b, ts->n)) record_hit("HIT"', 3),
                ('  return (u128)o * (z / gcd64(o, z)); }\n', '  return o; }\n', 1)):
    c = src.count(a)
    if c != k: sys.exit("count %d for %r" % (c, a))
    src = src.replace(a, b)
open(sys.argv[2], 'w').write(src)
print("hp2 ok", sys.argv[2])
PY
  gcc -O2 -fopenmp -o hp2/hp2_$v hp2/hp2_$v.c -lgmp -lm
done
for v in new old; do for XC in 20000:auto 20000:0 20000:-1 100000:auto 100000:0; do X=${XC%%:*}; cap=${XC#*:}; for q in "" --quiet; do
  if [ "$cap" = auto ]; then capopt=""; else capopt="--cap $cap"; fi
  out=hp2/${v}_X${X}_cap${cap}${q}.out
  timeout 1800 ./hp2/hp2_$v --D 5 --X $X --mode bpsw $capopt $q --threads 4 > $out 2>&1; echo "rc=$?" >> $out
done; done; done
touch chainF2.done
