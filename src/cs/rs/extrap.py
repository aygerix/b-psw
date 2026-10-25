"""Extrapolation of Experiments A/B to the full instance (pure arithmetic; inputs from sylow_T1500.out, expA_*.out).

Measured cost law (expA, real and random columns alike, success >= 0.8):
  GEN   needs N_gen(B)  = b * 2^(ceil((B + 2)/b) - 1) primes (k-tree, leaf lists of 2^b subsets);
  SPLIT needs N_split   = (r + 4) * N_gen(B - r) primes, r = rank of the elementary 2-subgroup handled by GF(2).
Full instance: B = 2525.41 bits, G[2] has rank 237, so the split works in G/G[2] (2288.41 bits). Pool N = 4451.
"""
import math
B, r, N = 2525.4088, 237, 4451
def ngen(B, b): return b * 2 ** (math.ceil((B + 2) / b) - 1)
print(" b   log2 N_gen   log2 N_split(r=237)   best r' in [0,237]   capacity at N=4451 (bits)   deficit (doublings)")
for b in (20, 30, 40, 50, 60):
    g = math.log2(ngen(B, b)); s = math.log2((r + 4) * ngen(B - r, b))
    best = min(range(0, r + 1), key=lambda k: (k + 4 if k else 1) * ngen(B - k, b))
    h = int(math.floor(math.log2(N / b))); cap = (h + 1) * b - 2
    print(f"{b:3d}   {g:9.1f}    {s:14.1f}          {best:5d}                {cap:6d}                  {g - math.log2(N):6.1f}")
print("\nSplit break-even: r = b log2(r + 4)  ->  b* =", round(r / math.log2(r + 4), 1))
print("Oracle bound (all Sylow parts with l <= L0 free): remaining bits 1773.5 (L0=3), 1471.6 (7), 1000.4 (31), 611.1 (101);")
print("   k-tree capacity at N = 4451 is <= %d bits even at b = 60" % ((int(math.floor(math.log2(N / 60))) + 1) * 60 - 2))
# exchange rate along Xmax at T = 1500 (min-cut instances, data/bpsw_instance_T1500_1e15.json -> _1e17.json)
dB = 2525.4088 - 2375.8483; dl = math.log2(4451 / 3296)
print(f"\nExchange rate rho = dB / dlog2 N along Xmax (1e15 -> 1e17, T = 1500): {dB:.1f} / {dl:.3f} = {dB/dl:.0f} bits per doubling")
print("A pool direction lowers the k-tree deficit iff rho < b (rho = 345 vs b <= 60).")
# lattice: enumeration time per trial (expB): d ~ 34: < 0.05 s, d ~ 44: 0.4-0.8 s, d ~ 54: 50 s
print(f"Enumeration growth d 44 -> 54: log2(50/0.6)/10 = {math.log2(50/0.6)/10:.2f} bits per dimension; "
      f"sieving floor 2^(0.292 d) at d = 2527: 2^{0.292*2527:.0f}")
