# choose plants for a --cap -1 (all branch A) control so that the total branch-A work stays within BUDGET.
# In plant mode every override prime is an anchor in every pass: count(p, pass) ~ p^2 / (2 lcm(M'(p), W_pass)).
# Greedy: add plants in increasing order of their own worst anchor count, keep if the exact total stays <= BUDGET.
# usage: python3 subset_capA.py IN OUT BUDGET
import sys, re, math
def lcm(a, b): return a // math.gcd(a, b) * b
inp, out, budget = sys.argv[1], sys.argv[2], float(sys.argv[3])
hdr, blocks, cur = [], [], None
for line in open(inp):
    if line.startswith('# plant '): cur = [line]; blocks.append(cur); continue
    if cur is None: hdr.append(line); continue
    cur.append(line)
info = []
for b in blocks:
    ov = {int(l.split()[1]): int(l.split()[2]) for l in b if l.startswith('O')}
    t = [l.split() for l in b if l.startswith('T')][0]
    info.append((b, ov, (int(t[4]), int(t[5]), int(t[6]))))
def total(sel):
    ov = {}; passes = {(1, 1, 1)}
    for b, o, sab in sel: ov.update(o); passes.add(sab)
    return sum(p * p // (2 * lcm(M, lcm(A, B))) for (S, A, B) in passes for p, M in ov.items())
order = sorted(info, key=lambda x: max(p * p // (2 * M) for p, M in x[1].items()))
sel = []
for x in order:
    if total(sel + [x]) <= budget: sel.append(x)
sel.sort(key=lambda x: blocks.index(x[0]))
with open(out, 'w') as f:
    for l in hdr: f.write(l)
    f.write('# cap -1 subset: estimated total branch-A count <= %g (%.3g), %d of %d plants\n' % (budget, total(sel), len(sel), len(blocks)))
    for b, o, s in sel:
        for l in b: f.write(l)
cls = {}
for b, o, s in sel:
    c = re.search(r'class=(\S+)', b[0]).group(1); cls[c] = cls.get(c, 0) + 1
print(out, len(sel), "of", len(blocks), "est %.3g" % total(sel), sorted(cls.items()))
