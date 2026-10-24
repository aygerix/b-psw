"""Discrete-log form of the BPSW subset-product instance, plus a generation check. Stdlib only.

Input: a subpool file (one prime per line, e.g. mincut.py output). Orders are recomputed with profile2.load.
A = lcm odd(ord_p 2) (F side), B = lcm odd(z(p)) (L side), G = (Z/A)^* x (Z/B)^* = prod_{q^e} Z/phi(q^e).
For each odd prime power q^e || AB, with g_q the least primitive root mod q^e:
  c[q^e][i] = log_{g_q}(p_i) mod phi(q^e),   target t = 0 (q^e | A) or phi(q^e)/2 (q^e | B).
A subset x in {0,1}^N gives a BPSW psp iff sum_i c[q^e][i] x_i == t (mod phi(q^e)) for every q^e, and |x| >= 3.
|x| is then odd automatically when 3 | B (every pool prime is == -1 mod 3).
Generation: <p_i> = G iff for every prime l | |G| the columns span G/lG = F_l^{r_l}, with r_l = #{q^e : l | phi(q^e)}.
Usage: python3 dlog_system.py subpool.txt out.txt [Tmax]
Output format (text):
  line 1: N R
  line 2: the N primes
  then for each of the R prime powers: "q e phi g side t" followed by a line of N residues.
"""
import sys, math, random
from profile2 import primes_upto, load

def prim_root(q, e):
    m = q ** e; ph = (q - 1) * q ** (e - 1)
    fs = [r for r in primes_upto(q) if (q - 1) % r == 0] + ([q] if e > 1 else [])
    for g in range(2, m):
        if g % q and all(pow(g, ph // r, m) != 1 for r in fs): return g
    raise ValueError

def rank_mod(rows, l):
    """Rank over F_l of a list of integer row vectors (entries reduced mod l)."""
    if l == 2:
        basis = {}                                   # leading bit -> row bitmask
        r = 0
        for row in rows:
            v = 0
            for j, a in enumerate(row):
                if a & 1: v |= 1 << j
            while v:
                h = v.bit_length() - 1
                if h in basis: v ^= basis[h]
                else: basis[h] = v; r += 1; break
        return r
    piv = []                                         # list of (pivot col, normalized row)
    r = 0
    for row in rows:
        v = [a % l for a in row]
        for c, pr in piv:
            if v[c]:
                f = v[c]; v = [(a - f * b) % l for a, b in zip(v, pr)]
        nz = next((j for j, a in enumerate(v) if a), None)
        if nz is None: continue
        inv = pow(v[nz], -1, l); v = [a * inv % l for a in v]
        piv.append((nz, v)); r += 1
    return r

def main():
    P = primes_upto(int(sys.argv[3]) if len(sys.argv) > 3 else 5000)
    rows = load(sys.argv[1], P)                      # (p, {q: e}) with q == 1 mod 4 (F) or 3 mod 4 (L)
    ps = [p for p, _ in rows]
    em = {}
    for _, m in rows:
        for q, e in m.items(): em[q] = max(em.get(q, 0), e)
    assert all(q not in em for q in ps), "unit condition violated"
    mods = sorted(em.items())
    out = [f"{len(ps)} {len(mods)}", " ".join(map(str, ps))]
    C = []
    for q, e in mods:
        m = q ** e; ph = (q - 1) * q ** (e - 1); g = prim_root(q, e)
        lg = {}; x = 1
        for k in range(ph): lg[x] = k; x = x * g % m
        side = 'F' if q % 4 == 1 else 'L'
        t = 0 if side == 'F' else ph // 2
        row = [lg[p % m] for p in ps]
        C.append((q, e, ph, g, side, t, row))
        out.append(f"{q} {e} {ph} {g} {side} {t}"); out.append(" ".join(map(str, row)))
    open(sys.argv[2], "w").write("\n".join(out) + "\n")
    logG = sum(math.log2(ph) for _, _, ph, *_ in C)
    print(f"system: N = {len(ps)} unknowns, R = {len(C)} prime-power rows, log2|G| = {logG:.4f}, "
          f"max phi = {max(ph for _, _, ph, *_ in C)}")
    # consistency check on random subsets: the residue sum must match log of the actual product
    random.seed(7)
    for _ in range(5):
        S = random.sample(range(len(ps)), 101)
        n = math.prod(ps[i] for i in S)
        for q, e, ph, g, side, t, row in C:
            assert pow(g, sum(row[i] for i in S) % ph, q ** e) == n % q ** e
    print("consistency: 5 random subsets match the products exactly")
    # generation certificate: ranks of G/lG projections
    order = math.prod(ph for _, _, ph, *_ in C)
    ells = [l for l in primes_upto(max(ph for _, _, ph, *_ in C)) if order % l == 0]
    bad = []
    for l in ells:
        sub = [row for q, e, ph, g, side, t, row in C if ph % l == 0]
        rk = rank_mod(sub, l)
        if rk < len(sub): bad.append((l, rk, len(sub)))
    print(f"generation: {len(ells)} primes l | |G|; projections to G/lG of full rank for all but {len(bad)}: {bad[:10]}")

if __name__ == "__main__":
    main()
