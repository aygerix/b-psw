"""l-primary profile of the instance group G = prod_{q^e} Z/phi(q^e) from a dlog_system.py output file.
For each prime l | |G|: rank r_l (number of coordinates with l | phi), bits = sum v_l(phi) log2 l, max exponent."""
import sys, math
from profile2 import primes_upto
f = open(sys.argv[1]); N, R = map(int, f.readline().split()); f.readline()
phis = []
for _ in range(R):
    q, e, ph, g, side, t = f.readline().split(); f.readline(); phis.append(int(ph))
prof = {}
for ph in phis:
    x = ph
    for l in primes_upto(int(math.isqrt(ph)) + 1) + [None]:
        if l is None:
            if x > 1: prof.setdefault(x, []).append(1)
            break
        if x % l == 0:
            v = 0
            while x % l == 0: x //= l; v += 1
            prof.setdefault(l, []).append(v)
tot = sum(math.log2(ph) for ph in phis)
print(f"N = {N}, R = {R}, log2|G| = {tot:.2f}")
rows = sorted(((l, len(vs), sum(vs) * math.log2(l), max(vs)) for l, vs in prof.items()), key=lambda r: -r[2])
for lo, hi in [(2, 2), (3, 7), (11, 50), (51, 200), (201, 10**9)]:
    b = sum(r[2] for r in rows if lo <= r[0] <= hi); n = sum(1 for r in rows if lo <= r[0] <= hi)
    print(f"  l in [{lo},{hi}]: {n} primes l, {b:.1f} bits ({100 * b / tot:.1f}%)")
print("  top l by bits (l: rank r_l, bits, max exponent):", "; ".join(f"{l}: {r}, {b:.0f}, {m}" for l, r, b, m in rows[:14]))
print("  large l with rank >= 2:", [(l, r) for l, r, b, m in rows if l > 50 and r >= 2][:20])
