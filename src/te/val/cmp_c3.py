# summarise plant runs: per class and per cap, found / missed / invalid. usage: python3 cmp_c3.py OUTFILE PLANTFILE
import sys, re, collections
outf, pf = sys.argv[1], sys.argv[2]
cls = []; cur = None
for line in open(pf):
    if line.startswith('# plant '): cur = re.search(r'class=(\S+)', line).group(1)
    elif line.startswith('T'): cls.append(cur); cur = None
txt = open(outf).read()
st = collections.defaultdict(lambda: [0, 0, 0, 0])   # found, missed, invalid, dup
seen = set(); bad = []
for line in txt.splitlines():
    m = re.match(r'PLANT (\d+) (\d+) (\d+) (\d+) S=(\S+) A=(\S+) B=(\S+) valid=(\d) found=(\d+)', line)
    if not m: continue
    k = int(m.group(1)); seen.add(k); c = cls[k]; v = int(m.group(8)); fd = int(m.group(9))
    if not v: st[c][2] += 1; bad.append(line)
    elif fd == 0: st[c][1] += 1; bad.append(line)
    else:
        st[c][0] += 1
        if fd > 1: st[c][3] += 1
summ = re.findall(r'PLANTS pass=\d+ valid=(\d+) found=(\d+) MISSED=(\d+) duplicates=(\d+) invalid=(\d+)', txt)
npass = len(re.findall(r'^PASS ', txt, re.M)); nhdr = re.search(r'(\d+) triples, (\d+) passes', txt)
ok = ('rc=0' in txt) and len(seen) == len(cls) and not bad and nhdr and npass == int(nhdr.group(2))
tot = [sum(st[c][i] for c in st) for i in range(4)]
print("%s plants=%d reported=%d passes=%d/%s found=%d missed=%d invalid=%d dup=%d %s" % (outf, len(cls), len(seen), npass, nhdr.group(2) if nhdr else '?', tot[0], tot[1], tot[2], tot[3], 'PASS' if ok else 'FAIL'))
print("   per class: " + "  ".join("%s %d/%d" % (c, st[c][0], st[c][0] + st[c][1] + st[c][2]) for c in sorted(st)))
for l in bad[:10]: print("   BAD:", l)
t = re.search(r'TIME (\S+) s wall (\S+) s user (\d+) KB', txt)
if t: print("   time: wall %s s, user %s s, maxrss %s KB" % t.groups())
