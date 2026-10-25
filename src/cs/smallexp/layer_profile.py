#!/usr/bin/env python3
"""Layer profile of a Pomerance-type BPSW instance and the best hybrid (LA layers + k-tree) element count.

Group G = (Z/L1)^* x (Z/L2)^* x Z/2 (parity), (Z/q^e)^* = Z/(q-1) x Z/q^(e-1).
Layer (l, k) has rank rho_{l,k} = #{cyclic factors with l^k | order}.  Killing a layer by the polynomial
methods costs a multiplicative factor f_l(rho) in elements (f_2 = rho+1, f_3 = (rho+1)^2, f_5 = (rho+1)^5)
and removes rho*log2(l) bits.  The rest is left to Wagner's k-tree with lists of 2^ell subset sums:
k = 2^(ceil(B'/ell) - 1) lists of ell super-elements each.  We report the minimal base-element count.
Also: GF(l)-rank of the pool's layer-1 images for l = 2, 3 (is the generated group full in the shallow layers?)
Usage: python3 layer_profile.py instance.json
"""
import sys, json, math
from collections import defaultdict
import numpy as np
from sympy import factorint, primitive_root

F = {2: lambda r: r + 1, 3: lambda r: (r + 1) ** 2, 5: lambda r: (r + 1) ** 5}

def cyclic_orders(fac):
    out = []
    for q, e in fac.items():
        q, e = int(q), int(e)
        out.append(q - 1)
        if e > 1:
            out.append(q ** (e - 1))
    return out

def layers(orders):
    rho = defaultdict(int)
    for m in orders:
        for l, a in factorint(m).items():
            for k in range(1, a + 1):
                rho[(l, k)] += 1
    return dict(rho)

def best_hybrid(rho, N, ells=range(8, 81)):
    B = sum(r * math.log2(l) for (l, k), r in rho.items())
    # candidate killed sets: prefixes per prime in {2,3,5}
    depth = {l: max(k for (ll, k) in rho if ll == l) for l in (2, 3, 5) if any(ll == l for ll, _ in rho)}
    best = None
    import itertools
    rng = [range(0, min(depth.get(l, 0), 4) + 1) for l in (2, 3, 5)]
    for d2, d3, d5 in itertools.product(*rng):
        logF, killed = 0.0, 0.0
        for l, d in ((2, d2), (3, d3), (5, d5)):
            for k in range(1, d + 1):
                r = rho[(l, k)]
                logF += math.log2(F[l](r)); killed += r * math.log2(l)
        Bp = B - killed
        for ell in ells:
            levels = math.ceil(Bp / ell)
            logN = (levels - 1) + math.log2(ell) + logF
            if best is None or logN < best[0]:
                best = (logN, (d2, d3, d5), ell, Bp)
    return B, best

def layer1_rank(primes, rows_sign, l):
    """GF(l) rank of the layer-1 images (l-th power characters) of the pool, rows with l | q - 1."""
    cols = []
    rows = [(q, s) for q, s in rows_sign if (q - 1) % l == 0]
    gens = {q: primitive_root(q) for q, _ in rows}
    tab = {}
    for q, _ in rows:
        z = pow(gens[q], (q - 1) // l, q)
        tab[q] = {pow(z, j, q): j for j in range(l)}
    M = np.zeros((len(primes), len(rows)), dtype=np.int64)
    for i, p in enumerate(primes):
        for j, (q, s) in enumerate(rows):
            x = (s * p) % q
            M[i, j] = tab[q][pow(x, (q - 1) // l, q)]
    import zsolve
    A = M.T.copy() % l
    piv = zsolve._rref(A, l)
    return len(rows), len(piv), M

def main():
    d = json.load(open(sys.argv[1]))
    primes = [int(t[0]) for t in d['primes']]
    N = len(primes)
    orders = cyclic_orders(d['L1_factors']) + cyclic_orders(d['L2_factors']) + [2]
    rho = layers(orders)
    B, best = best_hybrid(rho, N)
    print(f"N = {N}, log2|G| = {B:.1f}  (instance file says {d['log2_G']})")
    byl = defaultdict(float)
    for (l, k), r in rho.items():
        byl[l] += r * math.log2(l)
    print("bits by prime l (top 8):", ", ".join(f"{l}:{b:.0f}" for l, b in sorted(byl.items(), key=lambda t: -t[1])[:8]))
    for l in (2, 3, 5, 7):
        print(f"  l={l}: ranks by depth", [rho.get((l, k), 0) for k in range(1, 12) if rho.get((l, k), 0)])
    cheap1 = sum(rho.get((l, 1), 0) * math.log2(l) for l in (2, 3, 5))
    print(f"depth-1 bits at l in {{2,3,5}}: {cheap1:.0f} of {B:.0f} ({100*cheap1/B:.1f}%)")
    logN, ds, ell, Bp = best
    print(f"best hybrid: kill depths (2,3,5) = {ds}, ell = {ell}, residual bits {Bp:.0f}, "
          f"log2(base elements needed) = {logN:.1f} vs log2 N = {math.log2(N):.1f}: deficit {logN - math.log2(N):.1f} bits")
    B0, best0 = best_hybrid({k: v for k, v in rho.items()}, N, ells=[ell])
    # pure k-tree at the same ell and LA-only, for comparison
    lev = math.ceil(B / ell); print(f"pure k-tree at ell={ell}: log2 N needed = {lev - 1 + math.log2(ell):.1f}")
    for ell2 in (30, 40, 50, 64):
        lev = math.ceil(B / ell2); print(f"   pure k-tree ell={ell2}: log2 N = {lev - 1 + math.log2(ell2):.1f}")
    # shallow-layer generated group
    rs = [(int(q), 1) for q in d['L1_factors']] + [(int(q), -1) for q in d['L2_factors']]
    for l in (2, 3):
        nr, rk, _ = layer1_rank(primes, rs, l)
        print(f"layer-1 images, l={l}: {nr} rows with l | q-1, GF({l}) rank of the {N} pool images = {rk}")

if __name__ == "__main__":
    main()
