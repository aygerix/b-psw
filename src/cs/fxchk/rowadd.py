"""Adversarial FX check, part 3 claims (new code):
 (1) T0=300 instance for D=-15 (2,0) contains the whole T0=200 pool (superset claim);
 (2) the 112-prime T0=300 configuration and the 60-prime T0=200 configuration lie in the respective pools;
 (3) at the 60-prime configuration, how many extra T0=300 pool primes does ONE new row q in (200,300] (either side, any cap) add;
     and the best TWO rows."""
import sys, os, re, math, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from instlib import read_inst, item_q
B = os.path.expanduser('~/b-psw/src/cs/')
def load(files, D, c):
    pool = {}
    for fn in files:
        for l in open(fn):
            a = l.split()
            if len(a) == 7 and int(a[0]) == D and (int(a[3]), int(a[4])) == c: pool[int(a[1])] = (int(a[5]), int(a[6]))
    return pool
p200 = load([B + 'e1/c3_1e10.txt', B + 'e1/chm200_pool.txt'], -15, (2, 0))
p300 = load([B + 'fx/c3_T300_1e10.txt', B + 'e1/chm200_pool.txt'], -15, (2, 0))
nP, nI, nR, items, prim = read_inst(B + 'fx/inst300/inst_D-15_2_0.txt')
s300i = {p for p, _ in prim}
print(f"T0=200 pool {len(p200)}; T0=300 census+CHM200 pool {len(p300)}; inst300 file primes {nP}, rows {nR}; "
      f"T0=200 pool minus inst300 = {len(set(p200) - s300i)}; inst300 minus (census300+CHM) = {len(s300i - set(p300))}")
def parse_gp(fn):
    s = open(fn).read(); m = re.search(r'vcfg\((.*)\);', s, re.S); args = eval('[' + m.group(1) + ']')
    return args  # D, j1, j2, E1, E2, L
for fn, pool in ((B + 'fxchk/cfg/v_T200_60.gp', p200), (B + 'fxchk/cfg/v_T300_l25.gp', p300)):
    D, j1, j2, E1, E2, L = parse_gp(fn)
    print(f"{os.path.basename(fn)}: {len(L)} primes, in pool: {sum(1 for p in L if p in pool)}")
def fac(x):
    f = {}; q = 3
    while x > 1:
        while x % q == 0: f[q] = f.get(q, 0) + 1; x //= q
        q += 2
    return f
D, j1, j2, E1, E2, L = parse_gp(B + 'fxchk/cfg/v_T200_60.gp')
cap = [{q: e for q, e in E1}, {q: e for q, e in E2}]
req = {p: (fac(o1), fac(o2)) for p, (o1, o2) in p300.items()}
def usable(cap):
    out = []
    for p, (f1, f2) in req.items():
        if all(q in cap[0] and cap[0][q] >= a and q not in cap[1] for q, a in f1.items()) and \
           all(q in cap[1] and cap[1][q] >= a and q not in cap[0] for q, a in f2.items()): out.append(p)
    return out
base = usable(cap); print("usable at the 60-prime caps in the T0=300 pool:", len(base))
newrows = [q for q in range(201, 300, 2) if all(q % d for d in range(3, int(q ** .5) + 1, 2))]
single = []
for q in newrows:
    for s in (0, 1):
        for e in (1, 2):
            c2 = [dict(cap[0]), dict(cap[1])]; c2[s][q] = e
            add = len(usable(c2)) - len(base)
            bits = math.log2(q - 1) + (e - 1) * math.log2(q)
            single.append((add, q, 'FL'[s], e, round(bits, 2)))
single.sort(reverse=True)
print("best single new rows (added primes, q, side, cap, bits):", single[:8])
print("max added by one new row at cap 1:", max(a for a, q, s, e, b in single if e == 1), " histogram:",
      {k: sum(1 for a, *_ in single if a == k and _[2] == 1) for k in range(0, 6)})
# best pair of new rows (cap 1)
best = (0,)
cand = [(q, s) for a, q, s, e, b in single if e == 1 and a >= 1]
for (q1, s1), (q2, s2) in itertools.combinations(cand, 2):
    if q1 == q2: continue
    c2 = [dict(cap[0]), dict(cap[1])]; c2['FL'.index(s1)][q1] = 1; c2['FL'.index(s2)][q2] = 1
    add = len(usable(c2)) - len(base)
    if add > best[0]: best = (add, q1, s1, q2, s2, round(math.log2(q1 - 1) + math.log2(q2 - 1), 2))
print("best pair of new rows (cap 1):", best)
