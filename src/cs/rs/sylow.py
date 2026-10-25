"""Sylow decomposition of the BPSW subset-product instance, plus the component matrix used by the experiments.

Input: data/bpsw_instance_T1500_1e17.json (format: src/pom/instance.py). Rows are the prime powers q^e || L1 L2,
side F (q^e | L1, target n == 1) or L (q^e | L2, target n == -1). For row (q, e) with phi = (q-1) q^(e-1) and
g a primitive root mod q^e, the column of p is log_g(p) mod phi; by CRT Z/phi = sum over l^v || phi of Z/l^v.
A component is (l, v, q, e, side); its entry for p is log_g(p) mod l^v, its target is 0 (F) or phi/2 mod l^v
(= 0 for odd l, = 2^(v-1) for l = 2).

Prints, per prime l: number of components (rank of G_l), bits, exponent, and the rank of the pool columns in G/lG.
Also: fraction of entries that vanish (p == +-1 mod q forces the odd part of row q to 0) per l.
Writes an .npz with the component metadata, the integer matrix (components x primes), targets, primes and a
'forced' mask (entry forced by q | p^2 - 1).
Usage: python3 sylow.py instance.json out.npz
"""
import sys, json, math
import numpy as np

def factor(m):
    f, q = {}, 2
    while q * q <= m:
        while m % q == 0: f[q] = f.get(q, 0) + 1; m //= q
        q += 1
    if m > 1: f[m] = f.get(m, 0) + 1
    return f

def prim_root(q, e):
    m, ph = q ** e, (q - 1) * q ** (e - 1)
    fs = list(factor(ph))
    for g in range(2, m):
        if g % q and all(pow(g, ph // r, m) != 1 for r in fs): return g
    raise ValueError

def rank_mod_l(M, l):
    """Rank over F_l of an integer matrix (rows x cols), entries reduced mod l. Gaussian elimination in numpy."""
    A = (M % l).astype(np.int64)
    r, (nr, nc) = 0, A.shape
    for c in range(nc):
        piv = np.nonzero(A[r:, c])[0]
        if len(piv) == 0: continue
        i = r + piv[0]
        A[[r, i]] = A[[i, r]]
        A[r] = A[r] * pow(int(A[r, c]), -1, l) % l
        nz = np.nonzero(A[:, c])[0]
        nz = nz[nz != r]
        if len(nz): A[nz] = (A[nz] - np.outer(A[nz, c], A[r])) % l
        r += 1
        if r == nr: break
    return r

def main():
    inst = json.load(open(sys.argv[1]))
    ps = [int(t[0]) for t in inst["primes"]]
    N = len(ps)
    rows = [(int(q), e, 'F') for q, e in inst["L1_factors"].items()] + \
           [(int(q), e, 'L') for q, e in inst["L2_factors"].items()]
    rows.sort()
    meta, mat, tgt, forced = [], [], [], []
    for q, e, side in rows:
        m, ph = q ** e, (q - 1) * q ** (e - 1)
        g = prim_root(q, e)
        lg = np.zeros(m, dtype=np.int64); x = 1
        for k in range(ph): lg[x] = k; x = x * g % m
        col = lg[np.array([p % m for p in ps], dtype=np.int64)]
        frc = np.array([(p % q) in (1, q - 1) for p in ps])
        t = 0 if side == 'F' else ph // 2
        for l, v in factor(ph).items():
            mod = l ** v
            meta.append((l, v, q, e, side)); mat.append(col % mod); tgt.append(t % mod); forced.append(frc)
    mat = np.array(mat, dtype=np.int64); tgt = np.array(tgt, dtype=np.int64)
    np.savez_compressed(sys.argv[2], meta=np.array(meta, dtype=object), mat=mat, tgt=tgt,
                        primes=np.array([str(p) for p in ps]), forced=np.array(forced))
    # consistency: product of a random subset matches its residue sums
    rng = np.random.default_rng(1)
    S = rng.choice(N, 51, replace=False); n = math.prod(ps[i] for i in S)
    for (q, e, side) in rows[:40]:
        m, ph = q ** e, (q - 1) * q ** (e - 1); g = prim_root(q, e)
        comps = [i for i, mt in enumerate(meta) if mt[2] == q]
        # recombine log mod phi from the components by CRT
        res, mod = 0, 1
        for i in comps:
            l, v = meta[i][0], meta[i][1]; mm = l ** v; a = int(mat[i, S].sum()) % mm
            res = res + mod * ((a - res) * pow(mod, -1, mm) % mm); mod *= mm
        assert mod == ph and pow(g, res, m) == n % m
    print(f"N = {N} primes, {len(rows)} rows, {len(meta)} Sylow components; CRT consistency check passed")
    tot = sum(math.log2(l ** v) for l, v, *_ in meta)
    print(f"log2|G| = {tot:.4f} (instance file: {inst['log2_G']})")
    by = {}
    for i, (l, v, q, e, side) in enumerate(meta): by.setdefault(l, []).append((i, v, q, side))
    print(" l   rank  bits    exp   v-distribution        rank(G/lG) of columns  zero-frac  forced-frac  1/l")
    cum = 0.0
    out = []
    for l in sorted(by):
        lst = by[l]; bits = sum(v for _, v, _, _ in lst) * math.log2(l); cum += bits
        vs = {}
        for _, v, _, _ in lst: vs[v] = vs.get(v, 0) + 1
        idx = [i for i, *_ in lst]
        rk = rank_mod_l(mat[idx], l) if l < 50 or len(idx) > 1 else int(np.any(mat[idx] % l))
        zf = float(np.mean(mat[idx] % l == 0)) if l > 2 else float('nan')
        ff = float(np.mean(np.array(forced)[idx]))
        out.append((l, len(lst), bits, max(v for _, v, _, _ in lst)))
        print(f"{l:4d} {len(lst):5d} {bits:7.1f} {l}^{max(v for _, v, _, _ in lst):<3d} {str(dict(sorted(vs.items()))):22s} "
              f"{rk:5d}/{len(lst):<5d}        {zf:8.4f}   {ff:8.4f}  {1/l:.4f}")
    b2 = sum(b for l, r, b, x in out if l == 2); b3 = sum(b for l, r, b, x in out if l == 3)
    r2 = sum(r for l, r, b, x in out if l == 2); r3 = sum(r for l, r, b, x in out if l == 3)
    el2 = r2; el3 = r3 * math.log2(3)
    small = sum(b for l, r, b, x in out if l <= 7); big = tot - small
    print(f"\n2-part: {b2:.1f} bits, rank {r2} (elementary quotient G/2G: {el2} bits); 3-part: {b3:.1f} bits, rank {r3} "
          f"(G/3G: {el3:.1f} bits)")
    print(f"l <= 3: {b2 + b3:.1f} bits ({100*(b2+b3)/tot:.1f}%); l <= 7: {small:.1f} ({100*small/tot:.1f}%); "
          f"l > 7: {big:.1f} bits ({100*big/tot:.1f}%)")
    print(f"exponent-2-only part (G/2G) {el2} bits = {100*el2/tot:.1f}% ; exponent-3 quotient {el3:.1f} = {100*el3/tot:.1f}%")
    for L0 in (3, 5, 7, 13, 31, 101):
        hb = sum(b for l, r, b, x in out if l > L0)
        print(f"   bits in l > {L0:3d}: {hb:7.1f}   (rank-1 l's above {L0}: "
              f"{sum(1 for l, r, b, x in out if l > L0 and r == 1)} of {sum(1 for l, r, b, x in out if l > L0)})")

if __name__ == "__main__":
    main()
