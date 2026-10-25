"""Independent brute-force check of bbc.c on small instances (no min cut, no branch-and-bound): enumerate EVERY
configuration (each row: off, or side F / L with cap 1..E), count usable primes, and report
M(lambda) = max N - lambda bits, plus the exact frontier N_max(b) at integer b.  Compare with bbc's M_ub.
Usage: python3 bf_check.py instfile lam1,lam2,...      (instance format of bbc.c / mkinst.py)"""
import sys, itertools
import numpy as np

def main():
    fn = sys.argv[1]; lams = [float(x) for x in sys.argv[2].split(',')]
    f = open(fn); nP, nI, nR = map(int, f.readline().split())
    items = []
    for _ in range(nI):
        r, s, e, par, c = f.readline().split(); items.append((int(r), int(s), int(e), int(par), float(c)))
    reqs = []
    for _ in range(nP):
        a = f.readline().split(); k = int(a[1]); reqs.append([items[int(x)][:3] for x in a[2:2 + k]])
    E = {}
    for r, s, e, par, c in items: E[(r, s)] = max(E.get((r, s), 0), e)
    cost = {(r, s, e): c for r, s, e, par, c in items}
    opts = []
    for r in range(nR):
        o = [(None, 0, 0.0)]
        for s in (0, 1):
            acc = 0.0
            for e in range(1, E.get((r, s), 0) + 1):
                acc += cost[(r, s, e)]; o.append((s, e, acc))
        opts.append(o)
    # requirement matrix: for each prime and row, needed (side, a) or none
    need = np.full((nP, nR, 2), 0, dtype=np.int64)    # need[p, r, s] = exponent needed on side s (0 = none)
    for i, rq in enumerate(reqs):
        for r, s, a in rq: need[i, r, s] = a
    ncfg = 1
    for o in opts: ncfg *= len(o)
    print(f"{fn}: primes {nP}, rows {nR}, items {nI}, configurations {ncfg}", flush=True)
    best = {lam: (-1e300, None) for lam in lams}
    front = {}
    for cfg in itertools.product(*opts):
        cap = np.zeros((nR, 2), dtype=np.int64)
        bits = 0.0
        for r, (s, e, c) in enumerate(cfg):
            if s is not None: cap[r, s] = e; bits += c
        ok = np.all(need <= cap[None, :, :], axis=(1, 2))
        N = int(ok.sum())
        for lam in lams:
            v = N - lam * bits
            if v > best[lam][0]: best[lam] = (v, (N, bits))
        b0 = int(np.ceil(bits - 1e-9))
        if N > front.get(b0, -1): front[b0] = N
    for lam in lams:
        v, (N, b) = best[lam]
        print(f"lam {lam:.6f} BRUTE M = {v:.6f} at (bits {b:.3f}, N {N})")
    run = -1; fr = {}
    for b in range(0, max(front) + 1):
        run = max(run, front.get(b, -1)); fr[b] = run
    print("exact frontier N_max(b):", " ".join(f"{b}:{fr[b]}" for b in sorted(fr) if b % 10 == 0 or b == max(fr)))

if __name__ == "__main__":
    main()
