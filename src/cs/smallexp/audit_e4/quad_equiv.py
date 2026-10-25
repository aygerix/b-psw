#!/usr/bin/env python3
"""Sanity: for x_i = a_i + 2 b_i in (Z/4)^r, sum c_i x_i == 0 mod 4  <=>  L(c) = sum c_i a_i == 0 mod 2 and
Q_k(c) = sum_{i<j} c_i c_j a_ik a_jk + sum_i c_i b_ik == 0 mod 2.  And: Q is additive on unions of disjoint
blocks each with L = 0 (cross terms s_k(B) s_k(B') vanish mod 2) -- i.e. the layered algorithm is the
KPG-type linearisation of this MQ system on a subspace."""
import random
rng = random.Random(3)
bad = 0; tests = 0; add_bad = 0
for _ in range(3000):
    r = rng.randrange(1, 8); N = rng.randrange(2, 14)
    X = [[rng.randrange(4) for _ in range(r)] for _ in range(N)]
    A = [[v % 2 for v in x] for x in X]; Bv = [[v // 2 for v in x] for x in X]
    def L(c): return [sum(c[i] * A[i][k] for i in range(N)) % 2 for k in range(r)]
    def Q(c): return [(sum(c[i] * c[j] * A[i][k] * A[j][k] for i in range(N) for j in range(i + 1, N))
                       + sum(c[i] * Bv[i][k] for i in range(N))) % 2 for k in range(r)]
    c = [rng.randrange(2) for _ in range(N)]
    z = all(sum(c[i] * X[i][k] for i in range(N)) % 4 == 0 for k in range(r))
    q = not any(L(c)) and not any(Q(c))
    tests += 1; bad += (z != q)
    # additivity: two disjoint blocks with L = 0
    c1 = [rng.randrange(2) for _ in range(N)]; c2 = [0 if c1[i] else rng.randrange(2) for i in range(N)]
    if not any(L(c1)) and not any(L(c2)):
        cu = [c1[i] | c2[i] for i in range(N)]
        add_bad += any((Q(cu)[k] - Q(c1)[k] - Q(c2)[k]) % 2 for k in range(r))
print(f"equivalence mismatches {bad}/{tests}; additivity violations {add_bad}")
