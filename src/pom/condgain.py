"""W4: gain from conditional moduli inside a fixed subpool S*. For uniform random S (each prime w.p. 1/2),
Delta(S) = b(S*) - b(S) >= 0 bits are saved because S needs only its own levels. Reports log2 E[2^Delta]
(the factor by which the heuristic solution count within S* exceeds 2^(|S*|-1-b(S*))) and the Delta distribution.
Usage: python3 condgain.py subpool.txt [samples] [Tmax]"""
import sys, math, random
from profile2 import primes_upto, load, bits
P = primes_upto(int(sys.argv[3]) if len(sys.argv) > 3 else 2000)
rows = load(sys.argv[1], P)
def b_of(sel):
    em = {}
    for m in sel:
        for q, e in m.items(): em[q] = max(em.get(q, 0), e)
    return sum(bits(q, e) for q, e in em.items())
ms = [m for _, m in rows]; B0 = b_of(ms)
rng = random.Random(5); K = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
ds = []
for _ in range(K):
    ds.append(B0 - b_of([m for m in ms if rng.random() < 0.5]))
mx = max(ds); lg = mx + math.log2(sum(2 ** (d - mx) for d in ds) / K)
ds.sort()
print(f"S*: N = {len(ms)}, b = {B0:.3f}; random halves: Delta mean {sum(ds)/K:.2f}, median {ds[K//2]:.2f}, "
      f"max {mx:.2f}; log2 E[2^Delta] ~ {lg:.2f} (sampled, K = {K})")
# levels needed by few primes: they drive the gain
cnt = {}
for m in ms:
    for q, e in m.items():
        for ee in range(1, e + 1): cnt[(q, ee)] = cnt.get((q, ee), 0) + 1
em = {}
for m in ms:
    for q, e in m.items(): em[q] = max(em.get(q, 0), e)
top = [(q, em[q], cnt[(q, em[q])], bits(q, em[q]) - bits(q, em[q] - 1)) for q in em]
rare = sorted(t for t in top if t[2] <= 3)
print(f"top levels needed by <= 3 primes: {len(rare)}; exact independent-level estimate of the gain: "
      f"{sum(math.log2(1 + (2 ** w - 1) / 2 ** c) for _, _, c, w in top):.2f} bits")
