"""W5: do the arithmetic columns carry structure beyond the known support pattern? Stdlib only.

Input: a dlog_system.py output file. Coordinates are the prime powers q^e | AB, and c[r][i] = log of p_i mod q^e.
Matched random model (per prime i and coordinate q^e, independently):
  - side F (q == 1 mod 4): v = v_q(p_i - 1). If v >= 1, draw x uniform among units with x == 1 mod q^min(v,e) and,
    when v < e, x != 1 mod q^(v+1). If v = 0, draw x uniform among units with x != +-1 mod q.
  - side L: the same with p_i + 1 and -1.
  This keeps every forced residue (+-1 on the support, never -1 on the F side or +1 on the L side, exact q-adic
  depth) and randomises everything else.
Statistics on the projection to the first k coordinates (coordinates ordered by increasing phi), over all columns:
  P2 = #pairs with equal projection (signed weight-2 relations x_i - x_j),
  Z2 = #pairs with projections summing to 0 (signed: x_i + x_j == 0),
  T2 = #pairs summing to the target projection (unsigned weight-2 hits),
  T3 = #triples summing to the target projection (unsigned weight-3 hits), counted when C(N,3)/|H_k| is moderate.
Usage: python3 structtest.py system.txt [reps] [seed]
"""
import sys, math, random
from collections import Counter

def read_system(fn):
    f = open(fn)
    N, R = map(int, f.readline().split())
    ps = list(map(int, f.readline().split()))
    coords = []
    for _ in range(R):
        q, e, ph, g, side, t = f.readline().split()
        row = list(map(int, f.readline().split()))
        coords.append((int(q), int(e), int(ph), int(g), side, int(t), row))
    return ps, coords

def vq(x, q):
    v = 0
    while x % q == 0: x //= q; v += 1
    return v

def random_columns(ps, coords, rng):
    out = []
    for q, e, ph, g, side, t, row in coords:
        m = q ** e
        lg = {}; x = 1
        for k in range(ph): lg[x] = k; x = x * g % m
        s = 1 if side == 'F' else -1
        newrow = []
        for p in ps:
            v = vq(p - s, q)
            while True:
                if v >= 1:
                    a = min(v, e)
                    # x == s mod q^a: x = s + q^a * k
                    x = (s + q ** a * rng.randrange(q ** (e - a))) % m
                    if v < e and (x - s) % q ** (a + 1) == 0: continue
                else:
                    x = rng.randrange(1, m)
                    if x % q == 0 or (x - 1) % q == 0 or (x + 1) % q == 0: continue
                break
            newrow.append(lg[x])
        out.append((q, e, ph, g, side, t, newrow))
    return out

def stats(coords, ks, N):
    res = []
    order = sorted(range(len(coords)), key=lambda r: coords[r][2])
    for k in ks:
        sel = [coords[r] for r in order[:k]]
        bits = sum(math.log2(c[2]) for c in sel)
        proj = [tuple(c[6][i] for c in sel) for i in range(N)]
        mods = [c[2] for c in sel]; tgt = tuple(c[5] for c in sel)
        cnt = Counter(proj)
        P2 = sum(v * (v - 1) // 2 for v in cnt.values())
        neg = lambda h: tuple((-a) % m for a, m in zip(h, mods))
        sub = lambda u, h: tuple((a - b) % m for a, b, m in zip(u, h, mods))
        Z2 = sum(cnt[h] * cnt.get(neg(h), 0) for h in cnt)
        Z2 = (Z2 - sum(cnt[h] for h in cnt if neg(h) == h)) // 2
        T2 = sum(cnt[h] * cnt.get(sub(tgt, h), 0) for h in cnt)
        T2 = (T2 - sum(cnt[h] for h in cnt if sub(tgt, h) == h)) // 2
        T3 = None
        if math.comb(N, 3) / 2 ** bits < 5e4 and bits >= 20:
            # count triples i<j<l with h_i + h_j + h_l = tgt: for each pair, look up the third index
            T3 = 0
            idx = {}
            for i, h in enumerate(proj): idx.setdefault(h, []).append(i)
            for i in range(N):
                hi = proj[i]
                for j in range(i + 1, N):
                    need = tuple((a - b - c) % m for a, b, c, m in zip(tgt, hi, proj[j], mods))
                    for l in idx.get(need, ()):
                        if l > j: T3 += 1
        res.append((k, bits, P2, Z2, T2, T3))
    return res

def main():
    ps, coords = read_system(sys.argv[1])
    reps = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    N = len(ps)
    order = sorted(range(len(coords)), key=lambda r: coords[r][2])
    cum = 0.0; ks = []
    targets = [12, 16, 20, 24, 28, 32, 36, 40]
    for k, r in enumerate(order, 1):
        cum += math.log2(coords[r][2])
        if targets and cum >= targets[0]:
            ks.append(k)
            while targets and cum >= targets[0]: targets.pop(0)
    real = stats(coords, ks, N)
    rand = [stats(random_columns(ps, coords, rng), ks, N) for _ in range(reps)]
    print(f"N = {N}; projections onto the k smallest coordinates; real vs matched random ({reps} replicas, mean +- sd)")
    for t, (k, bits, P2, Z2, T2, T3) in enumerate(real):
        def ms(j):
            xs = [r[t][j] for r in rand if r[t][j] is not None]
            if not xs: return "   -"
            mu = sum(xs) / len(xs); sd = (sum((x - mu) ** 2 for x in xs) / max(1, len(xs) - 1)) ** 0.5
            return f"{mu:9.1f} +- {sd:7.1f}"
        print(f"k={k:3d} bits={bits:6.1f} | P2 {P2:8d} vs {ms(2)} | Z2 {Z2:8d} vs {ms(3)} | T2 {T2:7d} vs {ms(4)} | "
              f"T3 {T3 if T3 is not None else '-'} vs {ms(5)}", flush=True)

if __name__ == "__main__":
    main()
