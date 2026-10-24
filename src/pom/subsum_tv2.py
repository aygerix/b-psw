"""W5b, replicated: subset-sum distribution of real blocks vs many matched-random replicas (2-coordinate projections).

For each trial: a random pair of coordinates (q1^e1, q2^e2) with |H| = phi1 * phi2 in [2^14, 2^16], and a random
block of `blk` consecutive pool primes. The exact distribution of sum_{i in S} c_i over all 2^blk subsets S is
computed on H = Z/phi1 x Z/phi2 (2-D cyclic shifts). Statistics: TV distance from uniform and the normalised
collision mass C = |H| * sum_h P(h)^2 (1 for uniform). Each real value is compared with `reps` independent
matched-random replicas (structtest.random_columns), giving a z-score. Across trials, a systematic positive z on
C would mean extra mid-weight relations among the real columns.
Usage: python3 subsum_tv2.py system.txt [trials] [blk] [reps] [seed]
"""
import sys, math, random
from structtest import read_system, random_columns

def dist2(colA, colB, idx, m1, m2):
    D = [[0] * m2 for _ in range(m1)]; D[0][0] = 1
    for i in idx:
        s1, s2 = colA[i] % m1, colB[i] % m2
        rows = D[-s1:] + D[:-s1] if s1 else D
        new = []
        for r_old, r_sh in zip(D, rows):
            sh = r_sh[-s2:] + r_sh[:-s2] if s2 else r_sh
            new.append([a + b for a, b in zip(r_old, sh)])
        D = new
    return D

def stats(D, blk, size):
    tot = 2 ** blk
    tv = 0.5 * sum(abs(c / tot - 1 / size) for row in D for c in row)
    col = size * sum(c * c for row in D for c in row) / tot ** 2
    return tv, col

def main():
    ps, coords = read_system(sys.argv[1])
    trials = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    blk = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    reps = int(sys.argv[4]) if len(sys.argv) > 4 else 20
    rng = random.Random(int(sys.argv[5]) if len(sys.argv) > 5 else 5)
    N = len(ps)
    rcs = [random_columns(ps, coords, rng) for _ in range(reps)]
    print(f"N = {N}, block {blk}, {reps} matched-random replicas per trial; z = (real - mean) / sd", flush=True)
    zs_tv, zs_col = [], []
    for t in range(trials):
        while True:
            r1, r2 = rng.sample(range(len(coords)), 2)
            m1, m2 = coords[r1][2], coords[r2][2]
            if 2 ** 14 <= m1 * m2 <= 2 ** 16: break
        start = rng.randrange(0, N - blk); idx = range(start, start + blk)
        size = m1 * m2
        tr, cr = stats(dist2(coords[r1][6], coords[r2][6], idx, m1, m2), blk, size)
        sims = [stats(dist2(rc[r1][6], rc[r2][6], idx, m1, m2), blk, size) for rc in rcs]
        def z(x, xs):
            mu = sum(xs) / len(xs); sd = (sum((y - mu) ** 2 for y in xs) / (len(xs) - 1)) ** 0.5
            return (x - mu) / sd if sd > 0 else 0.0, mu, sd
        ztv, mtv, stv = z(tr, [s[0] for s in sims]); zc, mc, sc = z(cr, [s[1] for s in sims])
        zs_tv.append(ztv); zs_col.append(zc)
        q1 = coords[r1][0] ** coords[r1][1]; q2 = coords[r2][0] ** coords[r2][1]
        print(f"  H = ({q1}, {q2}) |H| = {size}: TV {tr:.4f} vs {mtv:.4f} +- {stv:.4f} (z {ztv:+.2f}); "
              f"C {cr:.4f} vs {mc:.4f} +- {sc:.4f} (z {zc:+.2f})", flush=True)
    m = lambda xs: sum(xs) / len(xs)
    print(f"mean z: TV {m(zs_tv):+.3f}, C {m(zs_col):+.3f} over {trials} trials "
          f"(sd of a mean of {trials} independent N(0,1) is {1 / trials ** 0.5:.3f})")

if __name__ == "__main__":
    main()
