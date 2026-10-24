"""C3: exact exhaustion of two-sided-smooth pools.

Pool P_y = { p prime : p^2 - 1 is y-smooth }  (from a CHM list of x with x, x+1 y-smooth; p = 2x + 1).
For a Selfridge D, every BPSW psp n with D* = D and ALL prime factors in P_y is a product of
 - inert p (eps = (D/p) = -1), or
 - admissible split p (Lemma 17: gcd(odd ord, odd z) = 1, min(v2 ord, v2 z) <= 1),
all in ONE 2-adic class (Lemma 1, F1), with the odd parts of the orders of different primes on
disjoint sides (F2), an odd number of inert primes (Lemma 17b), p not dividing D*Q*, and, when D = l*
is a prime discriminant, l not dividing ord (and not z if D > 0) (Lemma 12b).
We enumerate EVERY such subset S (independent sets of the side-conflict graph inside a class),
test the exact congruences of Lemma 1, and run the verbatim bpsw on survivors.
Usage: python3 exact.py y pairs_file D1 [D2 ...] [--cap N]
"""
import sys, os, math, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from sympy import isprime, factorint
from collections import defaultdict
from bpsw import bpsw_pass, jacobi, strong_fermat2, strong_lucas, selfridge_D

def v2(x): return (x & -x).bit_length() - 1
def odd(x): return x >> v2(x)

def order_mod(a, p, fac):
    o = p - 1
    for q, e in fac.items():
        for _ in range(e):
            if pow(a, o // q, p) == 1: o //= q
            else: break
    return o

def lucas_U(k, P, Q, p):
    """U_k(P,Q) mod p by the 2x2 matrix power."""
    def mul(A, B):
        return [[(A[0][0]*B[0][0] + A[0][1]*B[1][0]) % p, (A[0][0]*B[0][1] + A[0][1]*B[1][1]) % p],
                [(A[1][0]*B[0][0] + A[1][1]*B[1][0]) % p, (A[1][0]*B[0][1] + A[1][1]*B[1][1]) % p]]
    R = [[1, 0], [0, 1]]; M = [[P % p, (-Q) % p], [1, 0]]
    while k:
        if k & 1: R = mul(R, M)
        M = mul(M, M); k >>= 1
    return R[1][0]

def zrank(p, D, Q, eps, fac):
    o = p - eps
    for q, e in fac.items():
        for _ in range(e):
            if lucas_U(o // q, 1, Q, p) == 0: o //= q
            else: break
    return o

def prime_disc_l(D):
    l = abs(D)
    return l if isprime(l) else None

MODE = 'both'
MAXK = 0

def build(xs, D):
    Q = (1 - D) // 4
    l = prime_disc_l(D)
    rows = []
    for x in xs:
        p = 2 * x + 1
        if p < 3 or not isprime(p) or D % p == 0 or Q % p == 0:
            continue
        eps = jacobi(D % p, p)
        if eps == 0:
            continue
        fm = factorint(p - 1)
        o = order_mod(2, p, fm)
        fz = fm if eps == 1 else factorint(p + 1)
        z = zrank(p, D, Q, eps, fz)
        if l is not None and MODE != 'fermat':
            if o % l == 0 and MODE == 'both': continue    # Lemma 12(b)
            if D > 0 and z % l == 0: continue
        if MODE == 'fermat': z = 1
        if MODE == 'lucas': o = 1
        R1, R2 = factorint(odd(o)), factorint(odd(z))
        if eps == 1 and MODE == 'both':                    # Lemma 17(a): admissible split primes only
            if set(R1) & set(R2) or min(v2(o), v2(z)) > 1:
                continue
        rows.append(dict(p=p, eps=eps, o=o, z=z, cls=(v2(o), v2(z)), S1=frozenset(R1), S2=frozenset(R2),
                         O1=odd(o), O2=odd(z)))
    return rows

def run_class(R, cap):
    """Enumerate independent sets of the conflict graph; test Lemma 1 congruences exactly."""
    k = len(R)
    conf = [0] * k
    for i in range(k):
        for j in range(k):
            if i != j and ((R[i]['S1'] & R[j]['S2']) or (R[i]['S2'] & R[j]['S1'])):
                conf[i] |= 1 << j
    j1, j2 = R[0]['cls']
    cnt = [0, 0]; hits = []
    # DFS over indices in order; state: chosen mask, forbidden mask, product, lcm1, lcm2, inert parity
    stack = [(0, 0, 0, 1, 1, 1, 0)]
    while stack:
        start, chosen, forb, n, L1, L2, par = stack.pop()
        if chosen and (chosen & (chosen - 1)) and (par == 1 or MODE == 'fermat'):
            cnt[0] += 1
            if (n - 1) % L1 == 0 and (n + 1) % L2 == 0 and v2(n - 1) >= j1 and v2(n + 1) >= j2:
                cnt[1] += 1
                ok = bpsw_pass(n) if MODE == 'both' else (strong_fermat2(n) if MODE == 'fermat' else
                     (strong_lucas(n) and selfridge_D(n) == DCUR))
                hits.append((n, [R[i]['p'] for i in range(k) if chosen >> i & 1], ok))
            if cap and cnt[0] >= cap:
                return cnt, hits, False
        if MAXK and bin(chosen).count('1') >= MAXK:
            continue
        for i in range(start, k):
            if forb >> i & 1:
                continue
            r = R[i]
            stack.append((i + 1, chosen | (1 << i), forb | conf[i], n * r['p'], math.lcm(L1, r['O1']),
                          math.lcm(L2, r['O2']), par ^ (1 if r['eps'] == -1 else 0)))
    return cnt, hits, True

if __name__ == "__main__":
    args = [a for a in sys.argv[1:]]
    cap = 0
    if '--maxk' in args:
        i = args.index('--maxk'); MAXK = int(args[i + 1]); del args[i:i + 2]
    if '--mode' in args:
        i = args.index('--mode'); MODE = args[i + 1]; del args[i:i + 2]
    if '--cap' in args:
        i = args.index('--cap'); cap = int(args[i + 1]); del args[i:i + 2]
    y, f = int(args[0]), args[1]
    Ds = [int(a) for a in args[2:]]
    xs = [int(l) for l in open(f) if l.strip()]
    for D in Ds:
        DCUR = D
        t0 = time.time()
        rows = build(xs, D)
        by = defaultdict(list)
        for r in rows: by[r['cls']].append(r)
        tot_sets = tot_cong = 0; complete = True; allhits = []
        for c, R in sorted(by.items(), key=lambda t: -len(t[1])):
            cnt, hits, done = run_class(R, cap)
            tot_sets += cnt[0]; tot_cong += cnt[1]; complete &= done; allhits += hits
        nsplit = sum(1 for r in rows if r['eps'] == 1)
        print(f"y={y} D={D}: pool {len(rows)} ({nsplit} admissible split), classes {len(by)}, "
              f"max class {max(len(v) for v in by.values()) if by else 0}, consistent odd subsets {tot_sets}, "
              f"pass congruences {tot_cong}, complete {complete}, {time.time() - t0:.1f}s", flush=True)
        for h in allhits:
            print("  CONGRUENCE HIT", h[0], h[1], "exact test", h[2], flush=True)
