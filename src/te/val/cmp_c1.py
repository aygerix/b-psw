# compare engine hit sets with ref_out lists; usage: python3 cmp_c1.py REFDIR
import sys, os, re, glob
ref = sys.argv[1]
def refset(D, X, pred):
    fn = os.path.join(ref, "D%s_X%s_%s.txt" % (D, X, pred))
    if not os.path.exists(fn): return None
    s = set()
    for line in open(fn):
        if line.startswith('#') or not line.strip(): continue
        n, a, b, c = map(int, line.split()); assert n == a*b*c; s.add((n, a, b, c))
    return s
hitre = re.compile(r'^(FERMAT|LUCAS|HIT|WEAK|W1ONE) n = (\d+) = (\d+) \* (\d+) \* (\d+)\s+\[branch ([AC])\]')
def engset(fn):
    s = set(); txt = open(fn).read(); br = {'A': 0, 'C': 0}; bad = []
    for line in txt.splitlines():
        if line.startswith(('FERMAT', 'LUCAS', 'HIT', 'WEAK', 'W1ONE')):
            m = hitre.match(line)
            if not m: bad.append(line); continue
            n, a, b, c = map(int, m.group(2, 3, 4, 5)); assert n == a*b*c and a < b < c
            if (n, a, b, c) in s: bad.append("DUP " + line)
            s.add((n, a, b, c)); br[m.group(6)] += 1
    ok = ('rc=0' in txt) and ('PASS ' in txt)
    pl = [l for l in txt.splitlines() if l.startswith('PASS ')]
    return s, ok, br, bad, (pl[0] if pl else '')
allok = True
rows = []
for sub, predmap in (('c1', {'fermat': 'fermat', 'lucas': 'lucas', 'bpsw': 'bpsw', 'weak': 'weak'}), ('c1w', {'fermat': 'weakord', 'lucas': 'weakz'})):
    for fn in sorted(glob.glob(os.path.join(sub, '*.out'))):
        b = os.path.basename(fn)[:-4]
        m = re.match(r'D(-?\d+)_X(\d+)_(\w+)_cap(.+)$', b); D, X, mode, cap = m.groups()
        r = refset(D, X, predmap[mode])
        e, ok, br, bad, pl = engset(fn)
        if r is None: status = 'NOREF'
        else: status = 'PASS' if (ok and e == r and not bad) else 'FAIL'
        if status == 'FAIL': allok = False
        pm = re.search(r'Aprimes=(\d+) Acands=(\d+) Cprimes=(\d+) Cpairs=(\d+) Cinner=(\d+)', pl)
        rows.append((sub, D, X, mode, cap, status, len(e), None if r is None else len(r), br['A'], br['C'], len(bad), ok, pm.groups() if pm else None))
        if status == 'FAIL':
            print("FAIL detail", fn, "missing:", sorted(r - e)[:5], "extra:", sorted(e - r)[:5], "bad:", bad[:3], "ok:", ok)
for row in rows: print(*row)
print("ALL C1/C2 OK" if allok else "SOME FAILURES")
