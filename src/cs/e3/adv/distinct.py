import sys, math
for fn in sys.argv[1:]:
    seen = {}
    for l in open(fn):
        r = l.split()
        n = int(r[0])
        if n in seen: continue
        seen[n] = r
    inert = [r for r in seen.values() if len(r) > 11]
    s2 = sum(int(r[8]) for r in seen.values()); sl = sum(int(r[9]) for r in seen.values()); kor = sum(int(r[10]) for r in seen.values())
    exph = sum(1.0 / int(r[-1]) for r in inert)
    top = sorted(inert, key=lambda r: int(r[-1]))[:5]
    lr = sorted((math.log(int(r[-1])) / math.log(int(r[0])) if int(r[-1]) > 1 else 0.0) for r in inert)
    print(fn, "distinct n", len(seen), "inert", len(inert), "sprp2", s2, "slprp", sl, "korselt", kor)
    print("  sum over distinct inert n of 1/dLt = %.4e" % exph, " median log dL/log n %.3f" % lr[len(lr)//2])
    for r in top: print("   dLt", r[-1], "n", r[0], "rho", r[4], "d", r[5], "X", r[6], "D", r[7])
