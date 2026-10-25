"""Export a DSO census (D = 5) as E1 columns.
Usage: python3 dso_export.py ROWS_FILE X OUT
Columns: p ord z kron v2ord v2z R1=factor(odd(ord)) R2=factor(odd(z)) usable
  kron = (5/p) (-1 inert, +1 split, 0 for p = 5); usable = 1 iff p can divide a strong BPSW psp with D = 5 at all:
  inert, or split with gcd(odd ord, odd z) = 1, and 5 divides neither ord nor z (Lemma 12(b)).
In E1 terms: the Fermat-side requirement of p is R1 (each q^e | odd(ord) needs q^e | n - 1), the Lucas-side
requirement is R2 (q^e | n + 1), and p fixes the 2-adic class (v2ord, v2z).
"""
import sys, re, math
from sympy import factorint
fn, X, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
rows = sorted(tuple(map(int, t)) for t in re.findall(r'\[(\d+), (\d+), (\d+), (-?\d+)\]', open(fn).read()))
v2 = lambda x: (x & -x).bit_length() - 1
odd = lambda x: x >> v2(x)
fs = lambda n: '1' if n == 1 else '*'.join(f"{q}^{e}" if e > 1 else f"{q}" for q, e in sorted(factorint(n).items()))
nu = 0
with open(out, 'w') as f:
    f.write(f"# E2 doubly-small-order census, D = 5: ALL primes p with ord_p(2) <= {X} and z(p) <= {X} "
            f"(complete: gcd of prod Phi_d(2) and prod Phi_e(alpha,beta), src/cs/dso/dso_census.gp; "
            f"independently re-derived by direct sieve for all p <= 1e9 by src/cs/dso/dso_direct.c, identical; max p here {rows[-1][0]}).\n")
    f.write("# p ord z kron v2ord v2z R1=odd(ord) R2=odd(z) usable\n")
    for p, o, z, k in rows:
        us = int(bool(k != 0 and (k == -1 or math.gcd(odd(o), odd(z)) == 1) and o % 5 != 0 and z % 5 != 0))
        nu += us
        f.write(f"{p} {o} {z} {k} {v2(o)} {v2(z)} {fs(odd(o))} {fs(odd(z))} {us}\n")
print(f"wrote {len(rows)} rows ({nu} usable) to {out}")
