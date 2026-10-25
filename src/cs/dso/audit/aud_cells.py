"""Cell census on the audit's recomputed rows ("p o z s"): histogram of (ord, z) cell sizes; every cell of size >= 2
classified (sign pattern, shared odd prime). Also tests every odd-inert-count subset of every multi-prime cell with
bpsw.py.  usage: python3 aud_cells.py SRC rows1 [rows2 ...]"""
import sys, math, itertools
from collections import defaultdict, Counter
sys.path.insert(0, sys.argv[1]); from bpsw import bpsw_pass
v2 = lambda x: (x & -x).bit_length() - 1
odd = lambda x: x >> v2(x)
for fn in sys.argv[2:]:
    R = [tuple(map(int, l.split())) for l in open(fn)]
    C = defaultdict(list)
    for p, o, z, s in R: C[(o, z)].append((p, s))
    h = Counter(len(v) for v in C.values())
    multi = {c: v for c, v in C.items() if len(v) >= 2}
    kinds = Counter((tuple(sorted(s for _, s in v)), math.gcd(odd(c[0]), odd(c[1])) > 1) for c, v in multi.items())
    tested = passes = 0
    for c, v in multi.items():
        for r in range(2, len(v) + 1):
            for S in itertools.combinations(v, r):
                tested += 1; passes += bpsw_pass(math.prod(p for p, _ in S))
    print(f"{fn}: rows {len(R)}, cell sizes {sorted(h.items())}, multi-cells by (signs, shared odd prime): {dict(kinds)}; "
          f"subsets tested {tested}, bpsw passes {passes}")
