# estimated branch-A candidate count of a --cap -1 plant run: every override prime p is an anchor in EVERY pass,
# count ~ p^2 / (2 lcm(M'(p), W)), W = lcm(A, B) of the pass (plus the background of non-planted primes, ignored)
import sys, math
def lcm(a, b): return a // math.gcd(a, b) * b
ov = {}; passes = set()
for l in open(sys.argv[1]):
    if l.startswith('O'): _, p, M = l.split(); ov[int(p)] = int(M)
    elif l.startswith('T'):
        t = l.split(); passes.add((int(t[4]), int(t[5]), int(t[6])))
passes.add((1, 1, 1))
tot = 0
for (S, A, B) in passes:
    W = lcm(A, B); s = sum(p * p // (2 * lcm(M, W)) for p, M in ov.items()); tot += s
print(sys.argv[1], "passes", len(passes), "overrides", len(ov), "est. candidates %.3g" % tot)
