"""Group structure of the BPSW subset-product problem for annealed CHM pools.
G = (Z/A)^* x (Z/B)^*, A = prod q^e1, B = prod q^e2. Decompose into l-primary cyclic components,
(Z/q^e)^* = C_(q-1) x C_(q^(e-1)). Report for each l: number of components and bits, and the split
'structured' (l <= L0) versus 'cyclic' (larger l). The odd target is 0 in every odd l-component."""
import sys, math, json
from collections import defaultdict
from sympy import primerange, factorint
from pool import build, evaluate
from frontier import anneal_l

z = int(sys.argv[1]); rows = build(z, sys.argv[2]); L0 = int(sys.argv[3]) if len(sys.argv) > 3 else 13
bycls = defaultdict(list)
for r in rows: bycls[r['cls']].append(r)
qs = [q for q in primerange(3, z + 1) if q != 5]
for c in sorted(bycls, key=lambda c: -len(bycls[c]))[:2]:
    for lam in (0.1, 0.35):
        b = max((anneal_l(bycls[c], qs, lam, 3000, s) for s in range(3)), key=lambda t: t[0])
        N, lg, e1, e2, pool = evaluate(bycls[c], b[1])
        comp = defaultdict(list)   # l -> list of exponents f (component C_(l^f))
        for q, e in list(e1.items()) + list(e2.items()):
            for l, f in factorint(q - 1).items(): comp[l].append(f)
            if e >= 2: comp[q].append(e - 1)
        bits = {l: sum(f for f in fs) * math.log2(l) for l, fs in comp.items()}
        tot = sum(bits.values())
        two = bits.get(2, 0.0)
        struct = sum(v for l, v in bits.items() if 2 < l <= L0)
        cyc = sum(v for l, v in bits.items() if l > L0)
        print(f"z={z} class {c} lambda={lam}: N={N}, log2|G|={lg:.1f} (sum of components {tot:.1f})")
        print(f"   2-part {two:.1f} bits in {len(comp.get(2,[]))} comps (max 2^{max(comp.get(2,[0]))}); "
              f"odd l<={L0}: {struct:.1f} bits; odd l>{L0}: {cyc:.1f} bits")
        print("   per l: " + ", ".join(f"{l}:{len(fs)}x{max(fs)}" for l, fs in sorted(comp.items())))
