#!/usr/bin/env python3
"""Random-instance tests of zsolve: correctness (exact sums), elements consumed vs the deterministic bound.
Usage: python3 test_zsolve.py [quick|full] [seed]"""
import sys, time, random
import zsolve as z

def run(moduli, trials, rng, label, want_target=True):
    G = z.Group(moduli)
    B = z.bound(G)
    Nmax = int(min(B, 3_000_000))
    used_zs, used_t, sizes, t0 = [], [], [], time.time()
    for _ in range(trials):
        X = [tuple(rng.randrange(m) for m in G.m) for _ in range(Nmax)]
        S = z.Solver(G, X)
        I = S.zero_sum()
        assert I is not None and len(I) > 0 and len(set(I)) == len(I)
        tot = G.zero()
        for i in I: tot = G.add(tot, X[i])
        assert tot == G.zero(), "zero-sum check failed"
        used_zs.append(S.pulled); sizes.append(len(I))
        if want_target:
            u = tuple(rng.randrange(m) for m in G.m)
            X2 = X + [tuple(rng.randrange(m) for m in G.m) for _ in range(min(Nmax, 3_000_000) + 64)]
            S2 = z.Solver(G, X2)
            J = S2.target(u)
            assert J is not None and len(set(J)) == len(J)
            tot = u
            for i in J: tot = G.add(tot, X2[i])
            assert tot == G.zero(), "target check failed"
            used_t.append(S2.pulled)
    lay = [(l, k, len(a)) for l, k, a in G.layers()]
    mean = lambda v: sum(v) / len(v) if v else float('nan')
    print(f"{label:28s} layers {lay}  bound {B:>10}  pulled(zs) mean {mean(used_zs):9.1f} max {max(used_zs):7d}"
          f"  pulled(target) mean {mean(used_t):9.1f}  |I| mean {mean(sizes):7.1f}  {time.time()-t0:6.1f}s", flush=True)
    return mean(used_zs)

def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'quick'
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    cases = [
        ([2] * 10, 20), ([2] * 50, 20), ([2] * 200, 5),
        ([3] * 5, 20), ([3] * 20, 10), ([3] * 60, 3),
        ([4] * 5, 20), ([4] * 20, 5),
        ([5] * 2, 20), ([5] * 4, 10), ([5] * 6, 3),
        ([6] * 6, 10), ([9] * 4, 5), ([8] * 4, 10), ([12] * 4, 5), ([10] * 3, 5),
        ([2, 4, 8, 3, 9], 20), ([4, 12, 36], 10),
    ]
    if mode == 'full':
        cases += [([2] * 400, 2), ([3] * 100, 2), ([5] * 8, 2), ([5] * 10, 1), ([4] * 60, 2), ([6] * 20, 2),
                  ([25] * 2, 3), ([20] * 3, 2), ([30] * 3, 2)]
    for mods, tr in cases:
        lab = "x".join(f"Z{m}" for m in sorted(set(mods))) + f" ^{len(mods)}" if len(set(mods)) == 1 else str(mods)
        run(mods, tr, rng, lab)

if __name__ == "__main__":
    main()
