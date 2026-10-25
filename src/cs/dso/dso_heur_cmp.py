"""Compare the DSO census with the heuristic of dso_heur.c (random part p > X only).
Usage: python3 dso_heur_cmp.py ROWS_FILE X [K]   (runs ./dso_heur)
Reports: observed vs expected primes p > X by sign and by category (shared odd prime, ratio class (d/g, e/g)),
expected 2-cells (large-large from the C sums, small-large from E of every cell holding a prime <= X),
and every occupied cell whose observed number of large primes is improbable (P < 1e-3) under its E.
"""
import sys, re, math, subprocess
from collections import Counter, defaultdict
fn, X = sys.argv[1], int(sys.argv[2]); K = sys.argv[3] if len(sys.argv) > 3 else '16'
rows = [tuple(map(int, t)) for t in re.findall(r'\[(\d+), (\d+), (\d+), (-?\d+)\]', open(fn).read())]
rows = [r for r in rows if r[3] != 0]
v2 = lambda x: (x & -x).bit_length() - 1
odd = lambda x: x >> v2(x)
cells = defaultdict(lambda: [0, 0])     # (o, z, s) -> [n_small, n_large]
for p, o, z, s in rows:
    cells[(o, z, s)][p > X] += 1
cf = f'/tmp/dso_cells_{X}.txt'
with open(cf, 'w') as f:
    for (o, z, s) in cells:
        f.write(f"{o} {z} {s}\n")
out = subprocess.run(['./dso_heur', str(X), cf, K], capture_output=True, text=True).stdout
E = {}
for line in out.splitlines():
    if line.startswith('CELL'):
        t = line.split()
        E[(int(t[1]), int(t[2]), int(t[3]))] = (float(t[4][2:]), float(t[5][2:]))
    elif not line.startswith('LBIN'):
        print(line)
def cat(o, z, s):
    g = math.gcd(o, z); sh = int(math.gcd(odd(o), odd(z)) > 1)
    rc = f"{o // g}:{z // g}" if (o // g <= 4 and z // g <= 4) else 'generic'
    return ('inert' if s == -1 else 'split', sh, rc)
obs = Counter(cat(o, z, s) for p, o, z, s in rows if p > X)
print("OBSERVED p > X by category:")
for k, v in sorted(obs.items()):
    print(f"  OBS {k[0]} shared={k[1]} ratio={k[2]}  n={v}")
print("observed p > X: inert", sum(1 for r in rows if r[0] > X and r[3] == -1), " split", sum(1 for r in rows if r[0] > X and r[3] == 1))
sl = sum(E[c][0] for c, (ns, nl) in cells.items() if ns >= 1)
sl_in = sum(E[c][0] for c, (ns, nl) in cells.items() if ns >= 1 and c[2] == -1)
print(f"expected small-large 2-cells: {sl:.3f} (inert {sl_in:.4f})")
h = Counter((ns, nl) for (ns, nl) in cells.values())
print("observed cells by (n_small, n_large):", sorted(h.items()))
print("cells with >= 2 primes and their E (expected large primes in the cell):")
for c, (ns, nl) in sorted(cells.items()):
    if ns + nl >= 2:
        e, q = E[c]
        # P(at least nl large) for Bernoulli sum ~ Poisson(e) when terms are small
        pt = 1 - sum(math.exp(-e) * e ** j / math.factorial(j) for j in range(nl))
        print(f"  cell {c} small={ns} large={nl} E={e:.4g}  P(>= {nl} large) ~ {pt:.3g}")
# calibration: observed large primes by the E of their cell, against sum of E over all cells in that E-bin
eb = defaultdict(float)
for line in out.splitlines():
    if line.startswith('EBIN'):
        t = line.split(); eb[(t[1], int(t[2].split('=')[1]))] = float(t[3].split('=')[1])
ob = Counter()
for c, (ns, nl) in cells.items():
    if nl:
        ob[('inert' if c[2] == -1 else 'split', math.floor(math.log10(E[c][0])) if E[c][0] > 0 else -99)] += nl
print("calibration by cell E (log10 bin): expected = sum of E over ALL cells in the bin, observed = large primes in cells of the bin")
for k in sorted(set(eb) | set(ob)):
    print(f"  EB {k[0]} log10E={k[1]:4d}  expected {eb.get(k, 0):9.3f}  observed {ob.get(k, 0)}")
