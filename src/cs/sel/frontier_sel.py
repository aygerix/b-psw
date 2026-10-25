"""Pool economics for the selector line: greedy frontier (profile2.py logic) for the Pomerance pool alone and for the
pool plus index stragglers, and the k-tree gap at each frontier point.

Rows: inert p, class (1,2), with o1 = odd(ord_p 2) supported on F = {q == 1 mod 4, q != 5, q <= T} and
o2 = odd(z(p)) supported on L = {q == 3 mod 4, q <= T}. G(S) = prod_q (Z/q^e_q)^*, e_q = max exponent over S.

k-tree requirement (Wagner, 2^m lists, list size 2^(b/(m+1)) built from disjoint blocks of primes):
  N >= 2^m * b/(m+1), time ~ 2^(m + b/(m+1)).  N_req(b, C) = min over m with m + b/(m+1) <= C of 2^m b/(m+1).
gap(C) = min over frontier points of N_req(b, C)/N: the factor by which the usable-prime count must grow at fixed
group bits for a 2^C attack.  (Existence alone needs N >= b + 1.)
Usage: python3 frontier_sel.py T poolfile [stragglerfile ...]   (straggler files: first column p; duplicates ignored)
"""
import sys, math
from sympy import factorint, n_order

def fibmod(k, p):
    a, b = 0, 1
    for bit in bin(k)[2:]:
        c = a * ((2 * b - a) % p) % p; d = (a * a + b * b) % p
        a, b = (d, (c + d) % p) if bit == '1' else (c, d)
    return a

def zrank(p):
    z = p + 1
    for q, e in factorint(p + 1).items():
        for _ in range(e):
            if fibmod(z // q, p) == 0: z //= q
            else: break
    return z

def v2(x): return (x & -x).bit_length() - 1

def bits(q, e): return 0.0 if e == 0 else math.log2(q - 1) + (e - 1) * math.log2(q)

def row(p, T):
    if p < 7 or p % 5 not in (2, 3): return None
    o, z = n_order(2, p), zrank(p)
    if (v2(o), v2(z)) != (1, 2): return None
    F = factorint(o >> v2(o)); L = factorint(z >> v2(z))
    if not all(q % 4 == 1 and q != 5 and q <= T for q in F): return None
    if not all(q % 4 == 3 and q <= T for q in L): return None
    return {**F, **L}

def frontier(rows):
    alive = [True] * len(rows)
    cnt = {}
    for m in rows:
        for q, e in m.items(): cnt.setdefault(q, {}); cnt[q][e] = cnt[q].get(e, 0) + 1
    def emax(q): return max(cnt[q]) if cnt[q] else 0
    at = {}
    for i, m in enumerate(rows):
        for q, e in m.items(): at.setdefault((q, e), []).append(i)
    N = len(rows); B = sum(bits(q, emax(q)) for q in cnt)
    curve = [(B, N)]
    while N > 0:
        bestc = None
        for q in cnt:
            e = emax(q)
            if e == 0: continue
            drop = [i for i in at[(q, e)] if alive[i]]
            hit = {}
            for i in drop:
                for q2, e2 in rows[i].items():
                    if e2 == emax(q2): hit[q2] = hit.get(q2, 0) + 1
            dB = 0.0
            for q2, h in hit.items():
                e2 = emax(q2)
                if h == cnt[q2][e2]:
                    lower = [ee for ee in cnt[q2] if ee < e2 and cnt[q2][ee] > sum(
                        1 for i in drop if rows[i].get(q2) == ee)]
                    dB += bits(q2, e2) - bits(q2, max(lower) if lower else 0)
            gain = dB - len(drop)
            if bestc is None or gain > bestc[0]: bestc = (gain, q, e, drop)
        if bestc is None: break
        _, q, e, drop = bestc
        for i in drop:
            alive[i] = False
            for q2, e2 in rows[i].items():
                cnt[q2][e2] -= 1
                if cnt[q2][e2] == 0: del cnt[q2][e2]
        N -= len(drop); B = sum(bits(q2, emax(q2)) for q2 in cnt)
        curve.append((B, N))
    return curve

def nreq(b, C):
    best = None
    for m in range(0, 64):
        t = m + b / (m + 1)
        if t <= C:
            v = 2 ** m * b / (m + 1)
            if best is None or v < best: best = v
    return best

def report(name, curve):
    N0, B0 = curve[0][1], curve[0][0]
    print(f"{name}: N = {N0}, log2|G| = {B0:.1f}, N/B = {N0 / B0:.3f}")
    pts = [c for c in curve if c[1] >= 2 and c[0] >= 20]
    step = max(1, len(pts) // 14)
    print("  frontier (bits, N):", " ".join(f"({b:.0f},{n})" for b, n in pts[::step]))
    for C in (60, 70, 80, 100):
        g = min(((nreq(b, C) / n, b, n) for b, n in pts if nreq(b, C) is not None), default=None)
        if g: print(f"  C = {C:3d}: min gap N_req/N = {g[0]:.1f} at (bits {g[1]:.0f}, N {g[2]}), N_req = {g[0] * g[2]:.0f}")
    ex = min(((b + 1) / n, b, n) for b, n in pts)
    print(f"  existence: min (b+1)/N = {ex[0]:.2f} at (bits {ex[1]:.0f}, N {ex[2]})")

def main():
    T = int(sys.argv[1])
    pool = sorted({int(l.split()[0]) for l in open(sys.argv[2]) if l.strip()})
    extra = set()
    for fn in sys.argv[3:]:
        for l in open(fn):
            if l.strip(): extra.add(int(l.split()[0]))
    extra -= set(pool)
    rp = [r for r in (row(p, T) for p in pool) if r is not None]
    rx = [r for r in (row(p, T) for p in sorted(extra)) if r is not None]
    print(f"T = {T}: pool rows {len(rp)} (of {len(pool)}), straggler rows {len(rx)} (of {len(extra)} extra primes)")
    report("pool", frontier(rp))
    if rx: report("pool + stragglers", frontier(rp + rx))

if __name__ == "__main__":
    main()
