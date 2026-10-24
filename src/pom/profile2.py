"""Pomerance pool (D = 5, class (1,2)): true orders and greedy frontier. Faster profile.py, stdlib only.

Same model as profile.py. For pool primes p (p - 1 = 2A, p + 1 = 4B, A and B T-smooth under the mod-4 rule)
take o = ord_p 2 and z = z(p). The subset-product group of a sub-pool S is G(S) = prod_q (Z/q^e_q)^*, where
e_q is the largest exponent of q in odd(o) or odd(z) over S. Then log2 E(#solutions) ~ |S| - 1 - log2|G(S)|.
Greedy move: drop every prime attaining the top exponent of one modulus q, choosing q to maximise
|S| - log2|G(S)|. Gains are computed incrementally (exponent counts per modulus).
Usage: python3 profile2.py poolfile [Tmax]
"""
import sys, math

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

def factor_smooth(x, P):
    f = {}
    for q in P:
        if q * q > x: break
        if x % q == 0:
            e = 0
            while x % q == 0: x //= q; e += 1
            f[q] = e
    if x > 1: f[x] = f.get(x, 0) + 1          # x is prime here
    return f

def v2(x): return (x & -x).bit_length() - 1

def order_mod(a, p, fac):
    o = p - 1
    for q, e in fac.items():
        for _ in range(e):
            if pow(a, o // q, p) == 1: o //= q
            else: break
    return o

def fib_pair(n, p):
    # (F_n, F_{n+1}) mod p, iterative fast doubling
    a, b = 0, 1
    for bit in bin(n)[2:]:
        c = a * ((2 * b - a) % p) % p; d = (a * a + b * b) % p
        a, b = (d, (c + d) % p) if bit == '1' else (c, d)
    return a, b

def zrank(p, fac):
    z = p + 1
    for q, e in fac.items():
        for _ in range(e):
            if fib_pair(z // q, p)[0] == 0: z //= q
            else: break
    return z

def bits(q, e): return 0.0 if e == 0 else math.log2(q - 1) + (e - 1) * math.log2(q)

def load(fn, P):
    rows = []
    for p in sorted(int(l) for l in open(fn) if l.strip()):
        if p < 5: continue
        fm, fp = factor_smooth(p - 1, P), factor_smooth(p + 1, P)
        o, z = order_mod(2, p, fm), zrank(p, fp)
        if (v2(o), v2(z)) != (1, 2): continue
        F = factor_smooth(o >> v2(o), P); L = factor_smooth(z >> v2(z), P)
        assert all(q % 4 == 1 for q in F) and all(q % 4 == 3 for q in L), p
        rows.append((p, {**F, **L}))
    return rows

def main():
    P = primes_upto(int(sys.argv[2]) if len(sys.argv) > 2 else 5000)
    rows = load(sys.argv[1], P)
    alive = [True] * len(rows)
    # cnt[q][e] = number of alive rows with exponent e (>= 1) of q
    cnt = {}
    for _, m in rows:
        for q, e in m.items():
            cnt.setdefault(q, {}); cnt[q][e] = cnt[q].get(e, 0) + 1
    def emax(q): return max(cnt[q]) if cnt[q] else 0
    at = {}          # (q, e) -> list of row indices with exponent e of q
    for i, (_, m) in enumerate(rows):
        for q, e in m.items(): at.setdefault((q, e), []).append(i)
    N = len(rows); B = sum(bits(q, emax(q)) for q in cnt)
    print(f"pool {N} primes in class (1,2), true-order log2|G| = {B:.1f}, moduli {len(cnt)}, "
          f"N/log2|G| = {N / B:.3f}, N - 1 - log2|G| = {N - 1 - B:.1f}", flush=True)
    curve = [(B, N)]; best = (N - 1 - B, N, B)
    while N > 0:
        bestc = None
        for q in cnt:
            e = emax(q)
            if e == 0: continue
            drop = [i for i in at[(q, e)] if alive[i]]
            # how many rows at each (q', e'_max) does the drop remove?
            hit = {}
            for i in drop:
                for q2, e2 in rows[i][1].items():
                    if e2 == emax(q2): hit[q2] = hit.get(q2, 0) + 1
            dB = 0.0
            for q2, h in hit.items():
                e2 = emax(q2)
                if h == cnt[q2][e2]:                     # top level of q2 emptied: find next level alive
                    lower = [ee for ee in cnt[q2] if ee < e2 and cnt[q2][ee] > sum(
                        1 for i in drop if rows[i][1].get(q2) == ee)]
                    dB += bits(q2, e2) - bits(q2, max(lower) if lower else 0)
            gain = dB - len(drop)                        # change in (N - B) is -len(drop) + dB
            if bestc is None or gain > bestc[0]: bestc = (gain, q, e, drop)
        if bestc is None: break
        _, q, e, drop = bestc
        for i in drop:
            alive[i] = False
            for q2, e2 in rows[i][1].items():
                cnt[q2][e2] -= 1
                if cnt[q2][e2] == 0: del cnt[q2][e2]
        N -= len(drop); B = sum(bits(q2, emax(q2)) for q2 in cnt)
        curve.append((B, N))
        if N - 1 - B > best[0]: best = (N - 1 - B, N, B)
    pts = [c for c in curve if c[0] >= 40]
    mr = max(pts, key=lambda c: c[1] / c[0]) if pts else (1, 0)
    step = max(1, len(curve) // 16)
    print("frontier (log2|G|, N):", " ".join(f"({b:.0f},{n})" for b, n in curve[::step]))
    print(f"max N/log2|G| on frontier (bits >= 40): {mr[1] / mr[0]:.3f} at ({mr[0]:.0f} bits, N = {mr[1]})")
    print(f"best N - 1 - log2|G| = {best[0]:.1f} at (N, bits) = ({best[1]}, {best[2]:.1f})")

if __name__ == "__main__":
    main()
