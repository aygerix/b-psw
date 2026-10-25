#!/usr/bin/env python3
"""Independent audit of the exponent-2 claim (does not import zsolve for parts A-C).
A: exhaustive Davenport constant of F_2^r for r <= 4 (max length of a zero-sum-free sequence).
B: independent bitmask GF(2) elimination on r+1 random / adversarial vectors, r up to 2600.
C: the standard basis (r vectors) is zero-sum free (optimality witness), checked exhaustively for r <= 16.
D: re-run zsolve on F_2^r instances (random + adversarial) and check exact sums and pulled counts."""
import sys, random, itertools, time

def zero_sum_free(seq):
    n = len(seq)
    # subset XOR via DP over reachable sums with nonempty flag
    reach = set()
    for x in seq:
        new = {x}
        for s in reach:
            new.add(s ^ x)
        if 0 in new:
            return False
        reach |= new
    return True

def partA(rmax):
    for r in range(1, rmax + 1):
        t0 = time.time()
        # longest zero-sum-free sequence: check every sequence of length r+1 has a zero sum (multisets suffice)
        ok_r1 = all(not zero_sum_free(s) for s in itertools.combinations_with_replacement(range(1 << r), r + 1))
        exists_r = any(zero_sum_free(s) for s in itertools.combinations_with_replacement(range(1 << r), r))
        print(f"A r={r}: every (r+1)-multiset has nonempty zero-sum: {ok_r1}; some r-multiset zero-sum free: {exists_r}"
              f"  => D = {r+1 if ok_r1 and exists_r else '??'}  ({time.time()-t0:.1f}s)", flush=True)

def dependency(vecs):
    """Independent elimination: vecs are ints (bitmasks). Returns nonempty index set with XOR 0, or None."""
    basis = {}   # pivot bit -> (vector, combination as int bitmask over indices)
    for i, v in enumerate(vecs):
        comb = 1 << i
        while v:
            p = v.bit_length() - 1
            if p in basis:
                bv, bc = basis[p]
                v ^= bv; comb ^= bc
            else:
                basis[p] = (v, comb); break
        if v == 0:
            return [j for j in range(i + 1) if comb >> j & 1]
    return None

def partB(rng):
    for r in [1, 2, 10, 50, 200, 400, 1000, 2600]:
        trials = 200 if r <= 400 else 5
        worst = 0
        for t in range(trials):
            mode = t % 4
            if mode == 0:   vecs = [rng.getrandbits(r) for _ in range(r + 1)]
            elif mode == 1: vecs = [1 << (i % r) for i in range(r)] + [rng.getrandbits(r) | 1]          # basis + 1
            elif mode == 2: vecs = [rng.getrandbits(r) | (1 << (r - 1)) for _ in range(r + 1)]       # affine hyperplane
            else:           vecs = [(1 << (i + 1)) - 1 for i in range(r)] + [(1 << r) - 1]             # nested, repeat
            I = dependency(vecs)
            assert I is not None and len(I) > 0
            x = 0
            for i in I: x ^= vecs[i]
            assert x == 0
            worst = max(worst, max(I) + 1)
        print(f"B r={r}: {trials} trials, dependency always found among r+1 vectors, exact XOR = 0; max index used {worst}", flush=True)

def partC():
    for r in range(1, 17):
        basis = [1 << i for i in range(r)]
        assert zero_sum_free(basis)
    print("C: standard basis of F_2^r is zero-sum free for r = 1..16 (so D >= r+1)", flush=True)

def partD(rng):
    sys.path.insert(0, '..')
    import zsolve as z
    for r, trials in [(10, 20), (50, 20), (200, 5), (400, 2), (1000, 1)]:
        G = z.Group([2] * r)
        pulled, tpulled = [], []
        for t in range(trials):
            if t % 2 == 0:
                X = [tuple(rng.randrange(2) for _ in range(r)) for _ in range(r + 1)]
            else:  # adversarial: basis vectors repeated + all-ones, no zero element
                X = [tuple(1 if j == (i % r) else 0 for j in range(r)) for i in range(r)] + [tuple([1] * r)]
            S = z.Solver(G, X)
            I = S.zero_sum()
            assert I and len(set(I)) == len(I)
            tot = G.zero()
            for i in I: tot = G.add(tot, X[i])
            assert tot == G.zero()
            pulled.append(S.pulled)
            # target: u in span guaranteed by appending a basis
            Y = X + [tuple(1 if j == i else 0 for j in range(r)) for i in range(r)]
            u = tuple(rng.randrange(2) for _ in range(r))
            S2 = z.Solver(G, Y)
            J = S2.target(u)
            assert J is not None and len(set(J)) == len(J)
            tot = u
            for i in J: tot = G.add(tot, Y[i])
            assert tot == G.zero()
            tpulled.append(S2.pulled)
        print(f"D zsolve r={r}: zero-sum pulled {sorted(set(pulled))} (claim r+1={r+1}); target pulled range "
              f"[{min(tpulled)}, {max(tpulled)}]; all exact sums verified", flush=True)
    # adversarial: target with u outside span must fail (no false positive)
    r = 20
    G = z.Group([2] * r)
    X = [tuple([0] + [rng.randrange(2) for _ in range(r - 1)]) for _ in range(200)]
    u = tuple([1] + [0] * (r - 1))
    S = z.Solver(G, X)
    print(f"D target outside span (first coord 0 in all x, u_0 = 1): result {S.target(u)} (expected None)", flush=True)

if __name__ == "__main__":
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 12345)
    partA(4); partC(); partB(rng); partD(rng)
