# F1 validation analysis (run in src/te/valF on the VM). Prints a pass/fail table.
import glob, os, re, sys, subprocess, collections
rows = []
def row(group, name, ok, info=''): rows.append((group, name, 'PASS' if ok else 'FAIL', info))
def txt(fn): return open(fn).read() if os.path.exists(fn) else ''
# 0. self-test
t = txt('selftest.out'); m = re.search(r'SELFTEST (PASSED|FAILED): (\d+) failures, ([\d.]+)s', t)
row('selftest', 'anchor3 --selftest', bool(m) and m.group(1) == 'PASSED' and 'rc=0' in t, m.group(0) if m else 'no summary')
# 1. hit-print test
hitl = re.compile(r'^HIT n = (\d+) = ')
res = {}
for fn in sorted(glob.glob('hp2/*.out')):
    b = os.path.basename(fn)[:-4]; t = txt(fn)
    v, X, cap, q = re.match(r'(new|old)_X(\d+)_cap([^-]*?|-1)(--quiet)?$', b).groups()
    N = sum(int(x) for x in re.findall(r'^PASS .* hits=(\d+) ', t, re.M))
    lines = [l for l in t.splitlines() if l.startswith('HIT ')]
    ns = sorted(int(hitl.match(l).group(1)) for l in lines if hitl.match(l))
    res[(v, X, cap, bool(q))] = (N, len(lines), ns, 'rc=0' in t)
for (v, X, cap, q), (N, nl, ns, rc) in sorted(res.items()):
    name = 'hp_%s X=%s cap=%s %s' % (v, X, cap, 'quiet' if q else 'verbose')
    if v == 'new' or not q: ok = rc and nl == N and N > 0
    else: ok = rc and nl == 0 and N > 0          # old source + --quiet: the reviewed defect (lines lost)
    ref = res.get((v, X, cap, False))
    same = ref is not None and ref[2] == ns
    exp = 'HIT lines == hits' if (v == 'new' or not q) else 'expected defect: 0 HIT lines'
    row('F1 hit-print', name, ok and (same or (v == 'old' and q)), 'hits=%d HITlines=%d (%s)%s' % (N, nl, exp, '' if (v == 'old' and q) else (' same n-set as verbose' if same else ' n-set DIFFERS')))
for X in ('20000', '100000'):
    for cap in ('auto', '0', '-1'):
        a, b = res.get(('new', X, cap, True)), res.get(('old', X, cap, False))
        if a and b: row('F1 hit-print', 'new quiet vs old verbose X=%s cap=%s' % (X, cap), a[2] == b[2], '%d vs %d n' % (len(a[2]), len(b[2])))
# generic identity with the pre-fix outputs (timing fields stripped, line multisets)
strip = re.compile(r'\s(wall|time|cpu|threadsec_A|threadsec_C|us_per_Acand|us_per_Cpair)=\S+')
def norm(fn):
    out = []
    for l in txt(fn).splitlines():
        if l.startswith(('TIME', 'progress', 'Command ')): continue
        out.append(strip.sub('', l))
    return collections.Counter(out)
def ident(sub):
    fs = sorted(glob.glob(sub + '/*.out')); same = 0; diff = []
    for fn in fs:
        old = os.path.join('../val', fn)
        if not os.path.exists(old): diff.append(fn + ' (no old)'); continue
        if norm(fn) == norm(old): same += 1
        else: diff.append(fn)
    return len(fs), same, diff
# 2. calibration
t = txt('c4/bpsw_d5_1e7.out'); tq = txt('c4/bpsw_d5_1e7_quiet.out')
p = re.search(r'^PASS .*$', t, re.M); pq = re.search(r'^PASS .*$', tq, re.M)
okc = bool(p) and 'rc=0' in t and ' hits=0 ' in p.group(0)
row('calibration', 'bpsw D=5 X=1e7', okc and norm('c4/bpsw_d5_1e7.out') == norm('../val/c4/bpsw_d5_1e7.out'), (strip.sub('', p.group(0)) if p else 'no PASS') + ' | identical to pre-fix run')
row('calibration', 'bpsw D=5 X=1e7 --quiet', bool(pq) and 'rc=0' in tq and strip.sub('', pq.group(0)) == strip.sub('', p.group(0) if p else ''), 'PASS line equal to verbose run')
# 3. C1/C2
if os.path.exists('c1/DONE'):
    o = subprocess.run([sys.executable, 'cmp_c1.py', '../ref_out'], capture_output=True, text=True).stdout
    st = collections.Counter(); 
    for l in o.splitlines():
        f = l.split()
        if len(f) >= 6 and f[0] in ('c1', 'c1w'): st[(f[0], f[3] in ('bpsw', 'weak'), f[5])] += 1
    c1p = sum(v for (s, bw, r), v in st.items() if s == 'c1' and not bw and r == 'PASS'); c1t = sum(v for (s, bw, r), v in st.items() if s == 'c1' and not bw)
    c2p = sum(v for (s, bw, r), v in st.items() if s == 'c1' and bw and r == 'PASS'); c2t = sum(v for (s, bw, r), v in st.items() if s == 'c1' and bw)
    wp = sum(v for (s, bw, r), v in st.items() if s == 'c1w' and r == 'PASS'); wt = sum(v for (s, bw, r), v in st.items() if s == 'c1w')
    row('C1', 'one-torus fermat/lucas vs brute force (D 5,-7,13; X 2000/6000/20000; caps 0,-1,auto,40)', c1p == c1t == 72, '%d/%d' % (c1p, c1t))
    row('C1', 'W1-final-test variant vs weakord/weakz (caps 0,-1,auto)', wp == wt == 54, '%d/%d' % (wp, wt))
    row('C2', 'bpsw/weak vs reference', c2p == c2t == 72, '%d/%d' % (c2p, c2t))
    for sub in ('c1', 'c1w'):
        n, s, d = ident(sub); row('identity', sub + ' outputs vs pre-fix binary', n > 0 and s == n, '%d/%d identical %s' % (s, n, d[:3]))
if os.path.exists('stage4.done'):
    o = subprocess.run([sys.executable, 'cmp_split.py'], capture_output=True, text=True).stdout
    np_ = sum(1 for l in o.splitlines() if l.endswith(' PASS')); nt = sum(1 for l in o.splitlines() if l.endswith(' PASS') or l.endswith(' FAIL'))
    row('C1', 'split parts S > 1 vs split brute force (9 S, fermat/lucas, plain + W1 variant, caps 0,-1,auto)', np_ == nt == 108, '%d/%d %s' % (np_, nt, 'ALL SPLIT OK' in o))
    n, s, d = ident('c6'); row('identity', 'c6 outputs vs pre-fix binary', n > 0 and s == n, '%d/%d identical %s' % (s, n, d[:3]))
# 4. C3
pf = {}
for q in ('q3a.txt', 'q3c.txt', 'q8.txt'):
    for l in open(q):
        f = l.split()
        if not f: continue
        pf[f[0]] = f[f.index('--plants') + 1]
tot = collections.Counter()
for out, plants in pf.items():
    if not os.path.exists(out) or 'rc=' not in txt(out): continue
    o = subprocess.run([sys.executable, 'cmp_c3.py', out, plants], capture_output=True, text=True).stdout
    m = re.search(r'plants=(\d+) reported=(\d+) passes=(\d+)/(\S+) found=(\d+) missed=(\d+) invalid=(\d+) dup=(\d+) (PASS|FAIL)', o)
    cap = re.search(r'--cap (\S+)', open('q3a.txt').read() + open('q3c.txt').read() + open('q8.txt').read())
    capv = re.search(r'_cap(-1|0|auto)\.out$', out).group(1)
    ok = bool(m) and m.group(9) == 'PASS'
    if m: tot[capv + '_found'] += int(m.group(5)); tot[capv + '_tot'] += int(m.group(1))
    row('C3', out, ok, ('found %s/%s missed %s invalid %s dup %s passes %s/%s' % (m.group(5), m.group(1), m.group(6), m.group(7), m.group(8), m.group(3), m.group(4))) if m else o.strip()[:200])
for capv in ('auto', '0', '-1'):
    if tot[capv + '_tot']: row('C3', 'total cap=' + capv, tot[capv + '_found'] == tot[capv + '_tot'], '%d/%d found' % (tot[capv + '_found'], tot[capv + '_tot']))
if os.path.isdir('c3'):
    n, s, d = ident('c3'); row('identity', 'c3 outputs vs pre-fix binary', n > 0 and s == n, '%d/%d identical %s' % (s, n, d[:3]))
for r in rows: print(' | '.join(r))
print('SUMMARY', collections.Counter(r[2] for r in rows))
