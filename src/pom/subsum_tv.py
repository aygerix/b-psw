"""W5b: equidistribution of subset sums, real columns vs matched random columns (exact convolution). Stdlib only.

Choose a projection H = prod Z/phi(q^e) over a random set of coordinates with |H| in [2^14, 2^16]. Take a block of
`blk` consecutive pool primes (real columns) and compute exactly the distribution of sum_{i in S} c_i over all 2^blk
subsets S, projected to H. Report the total-variation distance from uniform, and the same for blocks of
matched-random columns (structtest.random_columns: exact q-adic depth on the support, forbidden +-1 off it).
If the arithmetic columns had exploitable mid-weight structure, their TV distance would be systematically larger.
Usage: python3 subsum_tv.py system.txt [trials] [blk] [seed]
"""
import sys, math, random
from structtest import read_system, random_columns

def tv_block(cols, idx, sel, mods):
    size = math.prod(mods)
    rad = [1] * len(mods)
    for j in range(len(mods) - 2, -1, -1): rad[j] = rad[j + 1] * mods[j + 1]
    dist = [0] * size; dist[0] = 1
    for i in idx:
        shift = [cols[r][6][i] % m for r, m in zip(sel, mods)]
        new = dist[:]
        # add the vector `shift` to every index (mixed radix with per-coordinate wraparound)
        for x in range(size):
            c = dist[x]
            if not c: continue
            y = 0; rem = x
            for j, m in enumerate(mods):
                d = rem // rad[j]; rem -= d * rad[j]
                y += ((d + shift[j]) % m) * rad[j]
            new[y] += c
        dist = new
    tot = 2 ** len(idx)
    return 0.5 * sum(abs(c / tot - 1 / size) for c in dist), size

def main():
    ps, coords = read_system(sys.argv[1])
    trials = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    blk = int(sys.argv[3]) if len(sys.argv) > 3 else 24
    rng = random.Random(int(sys.argv[4]) if len(sys.argv) > 4 else 3)
    N = len(ps)
    rc = random_columns(ps, coords, rng)
    print(f"N = {N}, block size {blk}: TV distance of the subset-sum distribution from uniform on H")
    for t in range(trials):
        while True:
            sel = rng.sample(range(len(coords)), rng.randint(2, 4))
            mods = [coords[r][2] for r in sel]
            if 2 ** 14 <= math.prod(mods) <= 2 ** 16: break
        start = rng.randrange(0, N - blk)
        idx = list(range(start, start + blk))
        tr, size = tv_block(coords, idx, sel, mods)
        tq, _ = tv_block(rc, idx, sel, mods)
        qs = [coords[r][0] ** coords[r][1] for r in sel]
        print(f"  H = {qs} (|H| = {size}, {math.log2(size):.1f} bits), primes {start}..{start + blk - 1}: "
              f"TV real {tr:.4f}, matched random {tq:.4f}", flush=True)

if __name__ == "__main__":
    main()
