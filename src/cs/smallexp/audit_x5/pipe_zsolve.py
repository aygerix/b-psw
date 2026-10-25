#!/usr/bin/env python3
"""Force zsolve through its full l = 5 pipeline: zero-valued QR blocks (the 'lucky' shortcut that ends most
small-r runs with |I| ~ 5) are discarded, so every output is a genuine ternary -> double -> double -> LA product.
Checks exact 0/1 zero-sums and targets (numpy, independent of zsolve's Group.add).
Usage: python3 pipe_zsolve.py seed"""
import sys, time
sys.path.insert(0, '..')
import numpy as np
import zsolve as z

DISC = [0]
_orig = z.QRRel._block
def _block_nozero(self):
    while True:
        blk = _orig(self)
        if blk[0].any():
            return blk
        DISC[0] += 1
z.QRRel._block = _block_nozero

def chk(M, X, I, u=None):
    s = np.zeros(len(M), dtype=np.int64) if u is None else np.array(u, dtype=np.int64)
    if I: s = s + X[np.array(I)].sum(axis=0)
    assert len(set(I)) == len(I) and not (s % M).any()

def run(mods, trials, g, with_target=True):
    G = z.Group(mods); B = z.bound(G); M = np.array(mods, dtype=np.int64)
    N = int(min(3 * B, 3_000_000))
    for _ in range(trials):
        DISC[0] = 0; t0 = time.time()
        X = g.integers(0, 10**9, size=(N, len(mods))) % M
        XL = [tuple(int(a) for a in row) for row in X]
        S = z.Solver(G, XL); I = S.zero_sum()
        assert I is not None, f"pool exhausted {mods}"
        chk(M, X, I)
        msg = f"{str(mods):20s} bound={B:9d} pulled={S.pulled:9d} ratio={S.pulled/B:.3f} |I|={len(I):7d} zero-blocks discarded={DISC[0]}"
        if with_target:
            DISC[0] = 0
            u = tuple(int(a) for a in g.integers(0, 10**9, size=len(mods)) % M)
            S2 = z.Solver(G, XL); J = S2.target(u)
            assert J is not None, "target: pool exhausted"
            chk(M, X, J, u)
            msg += f" | target pulled={S2.pulled} |J|={len(J)}"
        print(msg + f" [{time.time()-t0:.0f}s]", flush=True)

def main():
    g = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
    for r in range(3, 13):
        run([5] * r, 2 if r <= 10 else 1, g, with_target=(r <= 9))
    for mods in ([25, 25], [125], [50, 10], [75, 15], [20, 30, 45], [60, 30], [150], [125, 5], [625]):
        run(mods, 1, g)

if __name__ == "__main__":
    main()
