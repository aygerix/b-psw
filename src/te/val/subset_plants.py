# keep the plant blocks (comment + O lines + T line) whose branch-A candidate count candA <= maxcand
# usage: python3 subset_plants.py IN OUT MAXCAND [ONLY_S1]
import sys, re
inp, out, mx = sys.argv[1], sys.argv[2], float(sys.argv[3]); only_s1 = len(sys.argv) > 4 and sys.argv[4] == '1'
hdr, blocks, cur = [], [], None
for line in open(inp):
    if line.startswith('# plant '):
        cur = [line]; blocks.append(cur); continue
    if cur is None: hdr.append(line); continue
    cur.append(line)
keep = []
for b in blocks:
    c = int(re.search(r'candA=[^=]*=(\d+)', b[0]).group(1))
    t = [l for l in b if l.startswith('T')]; assert len(t) == 1 and b[-1].startswith('T')
    s1 = t[0].split()[4:7] == ['1', '1', '1']
    if c <= mx and (s1 or not only_s1): keep.append(b)
with open(out, 'w') as f:
    for l in hdr: f.write(l)
    f.write('# subset: candA <= %g%s, %d of %d plants\n' % (mx, ' S=1 only' if only_s1 else '', len(keep), len(blocks)))
    for b in keep:
        for l in b: f.write(l)
print(out, len(keep), "of", len(blocks), "candA sum %.3g" % sum(int(re.search(r'candA=[^=]*=(\d+)', b[0]).group(1)) for b in keep))
