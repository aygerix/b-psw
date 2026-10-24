"""(I) at small scale: can a designed prime set S beat an interval for counting twin S-units?

Input: a list of x with x and x + 1 both y-smooth (data/sn_200.txt, CHM closure, ~99% complete).
For S a subset of the primes <= y, N(S) = #{x in list : every prime of x(x+1) lies in S}.
Compare intervals S = {primes <= y'} with S chosen greedily (backward elimination from all primes <= y:
repeatedly drop the prime whose removal loses the fewest pairs), at equal |S| and at equal sum log2 q.
Usage: python3 sunit_design.py sn_file y
"""
import sys, math

def primes_upto(n):
    return [p for p in range(2, n + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]

def main():
    y = int(sys.argv[2]); P = primes_upto(y); idx = {p: i for i, p in enumerate(P)}
    masks = []
    for line in open(sys.argv[1]):
        if not line.strip(): continue
        x = int(line); m = 0
        for v in (x, x + 1):
            for p in P:
                if v % p == 0:
                    m |= 1 << idx[p]
                    while v % p == 0: v //= p
            assert v == 1, x
        masks.append(m)
    from collections import Counter
    cnt = Counter(masks)
    print(f"{len(masks)} pairs, {len(cnt)} distinct supports, {len(P)} primes")
    def N(Smask): return sum(c for m, c in cnt.items() if m & ~Smask == 0)
    # intervals
    print("intervals: s, max prime, sum log2 q, N")
    full = (1 << len(P)) - 1
    inter = {}
    for s in range(1, len(P) + 1):
        Sm = (1 << s) - 1; b = sum(math.log2(p) for p in P[:s]); n = N(Sm); inter[s] = (b, n)
        if s % 3 == 0 or s == len(P): print(f"  {s:3d} {P[s-1]:4d} {b:8.1f} {n:8d}")
    # greedy backward elimination
    print("greedy elimination: s, dropped prime, sum log2 q, N, N / N(interval of same s)")
    Sm = full
    for s in range(len(P), 1, -1):
        best = None
        for i in range(len(P)):
            if not (Sm >> i) & 1: continue
            n = N(Sm & ~(1 << i))
            if best is None or n > best[0]: best = (n, i)
        n, i = best; Sm &= ~(1 << i)
        b = sum(math.log2(P[j]) for j in range(len(P)) if (Sm >> j) & 1)
        if (s - 1) % 3 == 0 or s - 1 <= 6:
            print(f"  {s-1:3d} {P[i]:4d} {b:8.1f} {n:8d}  ratio {n / max(1, inter[s-1][1]):.3f}  "
                  f"S = {[P[j] for j in range(len(P)) if (Sm >> j) & 1] if s - 1 <= 12 else ''}")

if __name__ == "__main__":
    main()
