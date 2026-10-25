# compare c6 engine outputs with c6ref brute-force lists
import glob, os, re, sys
sys.path.insert(0, '.')
hitre = re.compile(r'^(FERMAT|LUCAS|W1ONE) n = (\d+) = (\d+) \* (\d+) \* (\d+) \* (\d+)\s+\[branch ([AC])\]')
def refset(fn):
    s = set()
    for l in open(fn):
        if l.startswith('#') or not l.strip(): continue
        s.add(tuple(map(int, l.split())))
    return s
allok = True
for fn in sorted(glob.glob('c6/*.out')):
    b = os.path.basename(fn)[:-4]
    D, X, S, mode, w1, cap = re.match(r'D(-?\d+)_X(\d+)_S(\d+)_(fermat|lucas)(_w1)?_cap(.+)$', b).groups()
    pred = {('fermat', None): 'fermat', ('lucas', None): 'lucas', ('fermat', '_w1'): 'ford', ('lucas', '_w1'): 'fz'}[(mode, w1)]
    r = refset('c6ref/D%s_X%s_S%s_%s.txt' % (D, X, S, pred))
    txt = open(fn).read(); e = set(); br = {'A': 0, 'C': 0}; bad = 0
    for l in txt.splitlines():
        if l.startswith(('FERMAT', 'LUCAS', 'W1ONE')):
            m = hitre.match(l)
            if not m: bad += 1; continue
            n, a, bb, c, s = map(int, m.group(2, 3, 4, 5, 6))
            if n != a*bb*c*s or s != int(S): bad += 1
            e.add((n, a, bb, c)); br[m.group(7)] += 1
    ok = 'rc=0' in txt and 'PASS ' in txt and e == r and bad == 0
    allok &= ok
    print(b, pred, 'engine', len(e), 'ref', len(r), 'A', br['A'], 'C', br['C'], 'PASS' if ok else 'FAIL')
    if not ok: print('   missing', sorted(r - e)[:4], 'extra', sorted(e - r)[:4], 'bad', bad)
print('ALL SPLIT OK' if allok else 'SPLIT FAILURES')
