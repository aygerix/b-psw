import sys, math, random, json
from collections import defaultdict
from sympy import primerange
from pool import build, evaluate

def anneal_l(rows, qs, lam, iters, seed, use_orders=True):
    rnd = random.Random(seed)
    assign = {q: rnd.choice('AB') for q in qs}
    N, lg, *_ = evaluate(rows, assign, use_orders); cur = N - lam * lg
    best = (cur, dict(assign), N, lg)
    for it in range(iters):
        T = 2.0 * (1 - it / iters) + 0.02
        q = rnd.choice(qs); old = assign[q]
        assign[q] = rnd.choice([c for c in 'ABX' if c != old])
        N, lg, *_ = evaluate(rows, assign, use_orders); val = N - lam * lg
        if val >= cur or rnd.random() < math.exp((val - cur) / T):
            cur = val
            if val > best[0]: best = (val, dict(assign), N, lg)
        else: assign[q] = old
    return best

if __name__ == "__main__":
    z = int(sys.argv[1]); rows = build(z, sys.argv[2]); iters = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
    bycls = defaultdict(list)
    for r in rows: bycls[r['cls']].append(r)
    qs = [q for q in primerange(3, z + 1) if q != 5]
    res = {}
    for c in sorted(bycls, key=lambda c: -len(bycls[c]))[:3]:
        pts = []
        for lam in (0.0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75):
            b = max((anneal_l(bycls[c], qs, lam, iters, s) for s in range(3)), key=lambda t: t[0])
            pts.append((round(b[3], 1), b[2], round(b[2] - b[3], 1)))
        print(f"z={z} class {c} (|class|={len(bycls[c])}) frontier (log2|G|, N, N-log2|G|):", sorted(set(pts)))
        res[str(c)] = sorted(set(pts))
    json.dump(res, open(f"../data/frontier_z{z}.json", "w"))
