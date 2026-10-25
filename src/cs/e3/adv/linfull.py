# Linear sets p_d = 1 + d m (d in S): member d is NOT dead iff N_d = d^(k-1) - prod_{e != d}(d - 2e) = 0 (Theorem U).
# Integral Korselt with m = lcm(S) t for all t  <=>  lcm(S) | sum(S).  Search S subset of [1, 40], k = 3..5.
import itertools, math
from math import prod
hits = []
for k in (3, 4, 5):
    cnt = 0
    for S in itertools.combinations(range(1, 41), k):
        if math.gcd(*S) != 1: continue
        L = math.lcm(*S)
        if sum(S) % L: continue
        cnt += 1
        full = [d for d in S if d ** (k - 1) == prod(d - 2 * e for e in S if e != d)]
        if full: hits.append((k, S, full))
    print("k", k, "integral-Korselt primitive sets", cnt)
print("sets with a non-dead member:", hits[:20], len(hits))
# also: over Q (no integrality), examples
ex = [(S, [d for d in S if d ** 2 == prod(d - 2 * e for e in S if e != d)]) for S in itertools.combinations(range(1, 13), 3)]
print("k=3 over Q, S in [1,12], with a non-dead member:", [e for e in ex if e[1]][:8])
