"""E2 analysis of a doubly-small-order census (rows [p, ord_p 2, z_D(p), (D/p)] from dso_census.gp).

Usage: python3 dso_analyze.py ROWS_FILE D X [KMAX]

1. Cells (ord, z): size histogram, every cell of size >= 2 listed; every subset of size >= 2 tested
   with src/bpsw.py.
2. Usability: inert p (always: ord | p-1, z | p+1 so gcd(ord, z) | 2), split p iff gcd(odd ord, odd z) = 1.
3. Pairwise-closed graph (task's sufficient condition, inert primes, one 2-adic class):
   p ~ q iff odd(o_q) | p-1, odd(o_p) | q-1, odd(z_q) | p+1, odd(z_p) | q+1. Maximal cliques listed; every odd
   sub-clique of size >= 3 is a strong BPSW psp (Lemma, cs_rank entry 9) and is tested.
4. Exact search, all subsets of size k <= KMAX of usable primes in one 2-adic class (v2 ord, v2 z),
   odd number of inert members. Exact criterion (squarefree n, D selected because (D/n) = -1 and D is the
   first Selfridge value with symbol -1 -- checked by bpsw.py at the end):
       n == 1 mod odd(o_q) and n == -1 mod odd(z_q) for every q in S
   (the 2-parts are automatic: n == 1 mod 2^a, n == -1 mod 2^b from the class and the parity).
   Anchor = member with the largest M_q = odd(o_q) odd(z_q); necessary condition at the anchor
   prod(S - anchor) == t_anchor mod M_anchor, then full verification, then bpsw.py.
"""
import sys, re, os, itertools, math, time
from collections import Counter, defaultdict
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from bpsw import bpsw_pass, selfridge_D, strong_fermat2
import networkx as nx

fn, D, X = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
KMAX = int(sys.argv[4]) if len(sys.argv) > 4 else 5
# DSO_CONTROL=1: planted positive control. Drop the Lucas side (z := 1), all primes usable, no parity/selection;
# the exact search then looks for strong base-2 pseudoprimes, which exist; every hit is checked with strong_fermat2
# and the search is compared with brute force on the small classes.
CONTROL = os.environ.get('DSO_CONTROL') == '1'
FOUND = set()
rows = [tuple(map(int, t)) for t in re.findall(r'\[(\d+), (\d+), (\d+), (-?\d+)\]', open(fn).read())]
v2 = lambda x: (x & -x).bit_length() - 1
odd = lambda x: x >> v2(x)
print(f"# D = {D}, X = {X}, census rows = {len(rows)}, max p = {max(r[0] for r in rows)}")
print(f"# inert {sum(r[3] == -1 for r in rows)}, split {sum(r[3] == 1 for r in rows)}, ramified {sum(r[3] == 0 for r in rows)}")
print(f"# p > X: {sum(r[0] > X for r in rows)}")

ntested = 0; nhits = 0
def test(S, why):
    global ntested, nhits
    n = math.prod(S)
    ntested += 1
    if CONTROL:
        if not strong_fermat2(n):
            print(f"CONTROL FAILURE: {S} satisfies the criterion but is not a strong base-2 psp")
        else:
            nhits += 1
        return True
    if bpsw_pass(n):
        nhits += 1
        print(f"HIT {why}: n = {n} = {' * '.join(map(str, S))}  selfridge D = {selfridge_D(n)}")
        return True
    return False

# 1. cells
cells = defaultdict(list)
for p, o, z, k in rows:
    cells[(o, z)].append((p, k))
h = Counter(len(v) for v in cells.values())
print("cell-size histogram:", sorted(h.items()))
for c, v in sorted(cells.items()):
    if len(v) >= 2 and not CONTROL:
        o, z = c
        sh = math.gcd(odd(o), odd(z))
        print(f"  cell (ord, z) = {c}: {v}  gcd(odd ord, odd z) = {sh}  {'USABLE' if sh == 1 or all(k == -1 for _, k in v) else 'unusable (shared odd prime)'}")
        ps = [p for p, _ in v]
        for r in range(2, len(ps) + 1):
            for S in itertools.combinations(ps, r):
                test(list(S), f"cell {c}")

# 2. usable primes, 2-adic classes
usable = []
for p, o, z, k in rows:
    if k == 0 or p == 2:
        continue
    if CONTROL:
        usable.append(dict(p=p, o=o, z=1, s=-1, a=v2(o), b=0, oo=odd(o), oz=1)); continue
    if k == 1 and math.gcd(odd(o), odd(z)) != 1:
        continue
    if D == 5 and (o % 5 == 0 or z % 5 == 0):
        continue     # Lemma 12(b): impossible in a BPSW psp with D = 5 (n == +-1 mod 5 is excluded)
    usable.append(dict(p=p, o=o, z=z, s=k, a=v2(o), b=v2(z), oo=odd(o), oz=odd(z)))
cls = defaultdict(list)
for u in usable:
    cls[(u['a'], u['b'])].append(u)
print(f"usable primes: {len(usable)} (inert {sum(u['s'] == -1 for u in usable)}, split {sum(u['s'] == 1 for u in usable)})")
print("2-adic classes (v2 ord, v2 z): #usable [#inert]:",
      sorted(((c, len(v), sum(u['s'] == -1 for u in v)) for c, v in cls.items()), key=lambda t: -t[1])[:20])

# 3. pairwise-closed graph on inert primes
tot_edges = 0; maxcl = 1
for c, mem in sorted(cls.items()):
    ine = [] if CONTROL else [u for u in mem if u['s'] == -1]
    G = nx.Graph()
    G.add_nodes_from(u['p'] for u in ine)
    for u, w in itertools.combinations(ine, 2):
        if (u['p'] - 1) % w['oo'] == 0 and (w['p'] - 1) % u['oo'] == 0 and (u['p'] + 1) % w['oz'] == 0 and (w['p'] + 1) % u['oz'] == 0:
            G.add_edge(u['p'], w['p'])
    tot_edges += G.number_of_edges()
    for cl in nx.find_cliques(G):
        if len(cl) >= 2:
            print(f"  closed clique class {c}: {sorted(cl)}")
        maxcl = max(maxcl, len(cl))
        for r in range(3, len(cl) + 1, 2):
            for S in itertools.combinations(sorted(cl), r):
                test(list(S), f"closed set class {c}")
print(f"pairwise-closed graph: {tot_edges} edges, max clique {maxcl}")

# 4. exact small-subset search
def crt_t(oo, oz, s):
    # t == 1 mod oo, t == -s mod oz  (oo, oz coprime)
    M = oo * oz
    t = (1 * oz * pow(oz, -1, oo) + (-s) * oo * pow(oo, -1, oz)) % M if M > 1 else 0
    return M, t

# Selfridge values preceding D: D is selected only if (D'/n) != -1 for each of them (squares never give -1)
EARLIER = []
_d = 5
while _d != D:
    if not (_d > 0 and math.isqrt(_d) ** 2 == _d):
        EARLIER.append(_d)
    _d = -_d - 2 if _d > 0 else -_d + 2
    assert abs(_d) <= abs(D) + 2
def jac(a, n):
    a %= n; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: r = -r
        a %= n
    return r if n == 1 else 0
def selected(n):
    return all(jac(dp, n) != -1 for dp in EARLIER) and jac(D, n) == -1
def full_ok(S):
    n = math.prod(u['p'] for u in S)
    if CONTROL:
        return all((n - 1) % u['oo'] == 0 for u in S)
    if sum(u['s'] == -1 for u in S) % 2 == 0:
        return False
    return all((n - 1) % u['oo'] == 0 and (n + 1) % u['oz'] == 0 for u in S) and selected(n)

t0 = time.time()
cnt_anchor = Counter(); cnt_full = Counter(); subsets_considered = Counter()
for c, mem in sorted(cls.items()):
    for u in mem:
        u['M'], u['t'] = crt_t(u['oo'], u['oz'], u['s'])
    mem.sort(key=lambda u: (u['M'], u['p']))
    m = len(mem)
    P = np.array([u['p'] for u in mem], dtype=object)
    for i in range(m):
        A = mem[i]; M = A['M']; t = A['t']
        prev = mem[:i]
        if not prev:
            continue
        if M == 1:
            # no constraint at the anchor; whole subset space below is checked by full_ok (tiny classes only)
            for k in range(2, KMAX + 1):
                for rest in itertools.combinations(prev, k - 1):
                    S = list(rest) + [A]
                    subsets_considered[k] += 1
                    cnt_anchor[k] += 1
                    if full_ok(S):
                        cnt_full[k] += 1; FOUND.add(tuple(sorted(u['p'] for u in S))); test(sorted(u['p'] for u in S), f"exact k={k} class {c}")
            continue
        res = np.array([u['p'] % M for u in prev], dtype=np.int64)
        okinv = np.array([math.gcd(u['p'], M) == 1 for u in prev])
        inv = np.array([pow(u['p'], -1, M) if math.gcd(u['p'], M) == 1 else 0 for u in prev], dtype=np.int64)
        bydres = defaultdict(list)
        for j, r in enumerate(res):
            if okinv[j]:
                bydres[int(r)].append(j)
        def check(idx, k):
            S = [prev[j] for j in idx] + [A]
            cnt_anchor[k] += 1
            if full_ok(S):
                cnt_full[k] += 1
                FOUND.add(tuple(sorted(u['p'] for u in S)))
                test(sorted(u['p'] for u in S), f"exact k={k} class {c}")
        # k = 2
        if KMAX >= 2:
            subsets_considered[2] += i
            for j in bydres.get(t % M, []):
                check([j], 2)
        # k = 3
        if KMAX >= 3:
            subsets_considered[3] += i * (i - 1) // 2
            for j in range(i):
                if not okinv[j]:
                    continue
                need = (t * int(inv[j])) % M
                for l in bydres.get(need, []):
                    if l > j:
                        check([j, l], 3)
        # k = 4
        if KMAX >= 4 and i >= 3:
            subsets_considered[4] += math.comb(i, 3)
            resset = np.array(sorted(bydres.keys()), dtype=np.int64)
            for j in range(i):
                if not okinv[j]:
                    continue
                base = (t * int(inv[j])) % M
                ls = np.arange(j + 1, i)
                ls = ls[okinv[j + 1:i]]
                if len(ls) == 0:
                    continue
                need = (base * inv[ls]) % M if M < 3037000499 else np.array([(base * int(inv[l])) % M for l in ls], dtype=np.int64)
                hit = np.isin(need, resset)
                for l, nd in zip(ls[hit], need[hit]):
                    for q in bydres[int(nd)]:
                        if q > l:
                            check([j, int(l), q], 4)
        # k = 5: meet in the middle on pair products
        if KMAX >= 5 and i >= 4:
            subsets_considered[5] += math.comb(i, 4)
            good = np.nonzero(okinv)[0]
            if len(good) >= 4:
                jj, ll = np.triu_indices(len(good), 1)
                jj = good[jj]; ll = good[ll]
                if M < 3037000499:
                    pr = (res[jj] * res[ll]) % M
                    pinv = (inv[jj] * inv[ll]) % M
                    need = (t * pinv) % M
                else:
                    pr = np.array([(int(res[a]) * int(res[b])) % M for a, b in zip(jj, ll)], dtype=np.int64)
                    need = np.array([(t * int(inv[a]) * int(inv[b])) % M for a, b in zip(jj, ll)], dtype=np.int64)
                seen5 = set()
                order = np.argsort(pr, kind='stable')
                prs = pr[order]
                lo = np.searchsorted(prs, need, 'left'); hi = np.searchsorted(prs, need, 'right')
                for x in np.nonzero(hi > lo)[0]:
                    a, b = int(jj[x]), int(ll[x])
                    for y in order[lo[x]:hi[x]]:
                        cc, d = int(jj[y]), int(ll[y])
                        # pairs {a,b}, {cc,d} disjoint; each 4-set is checked once
                        key = tuple(sorted((a, b, cc, d)))
                        if len(set(key)) == 4 and key not in seen5:
                            seen5.add(key)
                            check(list(key), 5)
print(f"exact search (k <= {KMAX}) done in {time.time() - t0:.1f} s")
for k in range(2, KMAX + 1):
    print(f"  k = {k}: subsets in scope {subsets_considered[k]:.4g}, anchor survivors {cnt_anchor[k]}, full solutions {cnt_full[k]}")
print(f"bpsw.py tests: {ntested}, passes: {nhits}" + ("  [CONTROL: passes = strong base-2 psp confirmed]" if CONTROL else ""))
if CONTROL:
    # brute force on small classes
    agree = 0; disagree = 0; nb = 0
    for c, mem in cls.items():
        if len(mem) > 70:
            continue
        for k in range(2, min(KMAX, 4) + 1):
            for S in itertools.combinations(mem, k):
                if full_ok(list(S)):
                    nb += 1
                    if tuple(sorted(u['p'] for u in S)) in FOUND: agree += 1
                    else: disagree += 1; print("CONTROL MISS:", sorted(u['p'] for u in S))
    print(f"CONTROL brute force (classes <= 70 members, k <= 4): {nb} solutions, found by search {agree}, missed {disagree}")
