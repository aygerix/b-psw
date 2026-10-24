"""Exact subpool optimisation for the Pomerance pool by minimum cut (maximum-weight closure).

The fixed-side model is the one in profile2.py: each prime p has odd orders a_p | odd(p-1) and b_p | odd(p+1), and
b(S) = sum over levels (q, e) needed by S of w(q, e), with w(q,1) = log2(q-1) and w(q,e) = log2 q for e >= 2.
For a fixed lambda, max_S |S| - lambda b(S) is a closure problem:
  source -> prime (cap 1); prime -> level (q, e_p(q)) (cap inf); level (q,e) -> level (q,e-1) (cap inf);
  level -> sink (cap lambda w(q,e)).
The source side of a minimum cut is an optimal S. Capacities are scaled by SCALE and rounded, so the optimum is
exact for the rounded weights; b(S) of the selected set is then recomputed exactly.
Unit condition (advisor review, C4): primes dividing the union modulus are removed first (iterated).
Usage: python3 mincut.py poolfile [Tmax] [outfile]
Prints the lambda = 1 optimum (max |S| - b(S)) and the max-ratio point (Dinkelbach iteration on lambda).
"""
import sys, math
from collections import deque
from profile2 import primes_upto, load, bits

SCALE = 10 ** 6
INF = 1 << 60

class Dinic:
    def __init__(self, n):
        self.n = n; self.g = [[] for _ in range(n)]
    def add(self, u, v, c):
        self.g[u].append([v, c, len(self.g[v])]); self.g[v].append([u, 0, len(self.g[u]) - 1])
    def bfs(self, s, t):
        self.lv = [-1] * self.n; self.lv[s] = 0; q = deque([s])
        while q:
            u = q.popleft()
            for v, c, _ in self.g[u]:
                if c > 0 and self.lv[v] < 0: self.lv[v] = self.lv[u] + 1; q.append(v)
        return self.lv[t] >= 0
    def dfs(self, u, t, f):
        if u == t: return f
        while self.it[u] < len(self.g[u]):
            e = self.g[u][self.it[u]]; v, c, r = e
            if c > 0 and self.lv[v] == self.lv[u] + 1:
                d = self.dfs(v, t, min(f, c))
                if d > 0:
                    e[1] -= d; self.g[v][r][1] += d; return d
            self.it[u] += 1
        return 0
    def maxflow(self, s, t):
        flow = 0
        while self.bfs(s, t):
            self.it = [0] * self.n
            while True:
                f = self.dfs(s, t, INF)
                if f == 0: break
                flow += f
        return flow
    def source_side(self, s):
        seen = [False] * self.n; seen[s] = True; q = deque([s])
        while q:
            u = q.popleft()
            for v, c, _ in self.g[u]:
                if c > 0 and not seen[v]: seen[v] = True; q.append(v)
        return seen

def clean_units(rows):
    while True:
        mods = {q for _, m in rows for q in m}
        keep = [r for r in rows if r[0] not in mods]
        if len(keep) == len(rows): return rows
        rows = keep

def exact_b(sel):
    em = {}
    for _, m in sel:
        for q, e in m.items(): em[q] = max(em.get(q, 0), e)
    return sum(bits(q, e) for q, e in em.items()), em

def solve(rows, lam):
    levels = {}
    for _, m in rows:
        for q, e in m.items():
            for ee in range(1, e + 1): levels.setdefault((q, ee), len(levels))
    n = 2 + len(rows) + len(levels); s, t = 0, 1
    D = Dinic(n)
    for i, (_, m) in enumerate(rows):
        D.add(s, 2 + i, SCALE)
        for q, e in m.items(): D.add(2 + i, 2 + len(rows) + levels[(q, e)], INF)
    for (q, e), j in levels.items():
        D.add(2 + len(rows) + j, t, round(lam * SCALE * bits(q, 1) if e == 1 else lam * SCALE * math.log2(q)))
        if e >= 2: D.add(2 + len(rows) + j, 2 + len(rows) + levels[(q, e - 1)], INF)
    D.maxflow(s, t)
    side = D.source_side(s)
    return [rows[i] for i in range(len(rows)) if side[2 + i]]

def main():
    P = primes_upto(int(sys.argv[2]) if len(sys.argv) > 2 else 5000)
    rows0 = [(p, m) for p, m in load(sys.argv[1], P)]
    rows = clean_units(rows0)
    b0, _ = exact_b(rows)
    print(f"pool {len(rows0)} -> {len(rows)} after unit cleaning (removed {sorted(set(r[0] for r in rows0) - set(r[0] for r in rows))}); "
          f"full b = {b0:.3f}, ratio {len(rows) / b0:.4f}", flush=True)
    S1 = solve(rows, 1.0); b1, _ = exact_b(S1)
    print(f"lambda = 1: N = {len(S1)}, b = {b1:.4f}, N - 1 - b = {len(S1) - 1 - b1:.2f}, ratio {len(S1) / max(b1, 1e-9):.4f}", flush=True)
    # Dinkelbach for the max ratio N/b (restricted to b >= 40 bits as in profile.py, by starting from S1 or the pool)
    lam = len(rows) / b0; best = None
    for it in range(60):
        S = solve(rows, lam); b, _ = exact_b(S)
        if not S or b < 40: break
        r = len(S) / b
        print(f"  Dinkelbach it {it}: lambda {lam:.6f} -> N = {len(S)}, b = {b:.4f}, ratio {r:.6f}", flush=True)
        if best is None or r > best[0] + 1e-12: best = (r, S, b)
        if abs(r - lam) < 1e-9: break
        lam = r
    if best:
        r, S, b = best
        print(f"max ratio: N = {len(S)}, b = {b:.4f}, ratio {r:.6f}")
    if len(sys.argv) > 3:
        open(sys.argv[3], "w").write("".join(f"{p}\n" for p, _ in sorted(S1)))
        print(f"wrote lambda = 1 subpool ({len(S1)} primes) to {sys.argv[3]}")

if __name__ == "__main__":
    main()
