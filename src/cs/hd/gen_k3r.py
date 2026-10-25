# Generate Singular files: k = 3 Korselt families in root form, f_j = 1 + c_j A(x) B_j(x), A = prod (x - al_m) (a roots, common
# to all f_j - 1), B_j = prod (x - th_{j,m}) (b_j roots), c_1 = 1, al_1 = 0 (translation) or th_{1,1} = 0 if a = 0.
# Korselt: f_i(th) f_l(th) = 1 at every root th of B_j.  Lucas incidences: singles (f_j(s) = -1, f_i f_l (s) = 1), triples.
# Saturate by c2 c3, all root differences, sigma differences, sigma - roots.  Output dim (-1 = no solution over C).
import itertools, sys, os
Dmax = int(sys.argv[1]); outdir = sys.argv[2]; mode = sys.argv[3] if len(sys.argv) > 3 else 'both'
EXPL = sys.argv[4] if len(sys.argv) > 4 else ''
os.makedirs(outdir, exist_ok=True)
jobs = []
for a in range(0, 4):
    for b in itertools.combinations_with_replacement(range(0, 6), 3):
        D = 3*a + sum(b)
        if D > Dmax or D < 3: continue
        if a == 0 and b[0] == 0: continue
        d = [a + x for x in b]; dmin = min(d)
        targets = set()
        if mode in ('both', 'hd'): targets.add(D//2 + 1)
        if mode in ('both', 'a2'): targets.add(D - 2*dmin + 1)
        pats = []
        for t in range(0, dmin + 1):
            for s in itertools.product(*[range(0, d[j] - t + 1) for j in range(3)]):
                if sum(s) + 3*t in targets: pats.append((s, t))
        # members with equal shape are symmetric: canonicalise s within equal-b blocks
        seen = set(); P2 = []
        for s, t in pats:
            key = []
            for j in range(3): key.append((b[j], s[j]))
            ck = (tuple(sorted(key)), t)
            if ck in seen: continue
            seen.add(ck); P2.append((s, t))
        for s, t in P2:
            jobs.append((a, b, s, t, D))
if EXPL:
    jobs = []
    for spec in EXPL.split(';'):
        a_, b_, s_, t_ = spec.split('/')
        b_ = tuple(int(c) for c in b_); s_ = tuple(int(c) for c in s_)
        jobs.append((int(a_), b_, s_, int(t_), 3*int(a_) + sum(b_)))
for (a, b, s, t, D) in jobs:
    name = f"k3_a{a}_b{''.join(map(str,b))}_s{''.join(map(str,s))}_t{t}"
    vars_ = []
    al = [f"al{m}" for m in range(1, a+1)]
    th = [[f"th{j+1}_{m}" for m in range(1, b[j]+1)] for j in range(3)]
    nsig = sum(s) + t
    sg = [f"sg{q}" for q in range(1, nsig+1)]
    fixed = {}
    if a >= 1: fixed[al[0]] = "0"
    else: fixed[th[0][0]] = "0"
    allv = [v for v in al + sum(th, []) + ["c2", "c3"] + sg if v not in fixed]
    L = []
    L.append(f"// {name}: shape a={a} b={b}, singles {s}, triples {t}, D={D}")
    L.append('LIB "elim.lib";')
    L.append('proc SAT(ideal I, poly g) { def S = sat(I, g); if (typeof(S) == "list") { ideal J = S[1]; return(J); } return(S); }')
    L.append(f"ring R = 0, ({', '.join(allv)}, TT), dp;")
    def r(v): return fixed.get(v, v)
    Aexpr = "*".join([f"(x-{r(v)})" for v in al]) or "1"
    L.append("proc Aof(poly x) { return(" + Aexpr + "); }")
    for j in range(3):
        Bexpr = "*".join([f"(x-{r(v)})" for v in th[j]]) or "1"
        L.append(f"proc B{j+1}(poly x) {{ return({Bexpr}); }}")
    L.append("proc f1(poly x) { return(1 + Aof(x)*B1(x)); }")
    L.append("proc f2(poly x) { return(1 + c2*Aof(x)*B2(x)); }")
    L.append("proc f3(poly x) { return(1 + c3*Aof(x)*B3(x)); }")
    eqs = []
    fn = ["f1", "f2", "f3"]
    for j in range(3):
        o = [l for l in range(3) if l != j]
        for v in th[j]:
            eqs.append(f"{fn[o[0]]}({r(v)})*{fn[o[1]]}({r(v)}) - 1")
    q = 0
    for j in range(3):
        o = [l for l in range(3) if l != j]
        for _ in range(s[j]):
            v = sg[q]; q += 1
            eqs.append(f"{fn[j]}({v}) + 1"); eqs.append(f"{fn[o[0]]}({v})*{fn[o[1]]}({v}) - 1")
    for _ in range(t):
        v = sg[q]; q += 1
        for j in range(3): eqs.append(f"{fn[j]}({v}) + 1")
    L.append("ideal I = " + ", ".join(eqs) + ";")
    roots = [r(v) for v in al] + [r(v) for v in sum(th, [])]
    nzf = ["c2", "c3"]
    for i1 in range(len(roots)):
        for i2 in range(i1):
            if roots[i1] != roots[i2]: nzf.append(f"({roots[i1]}-{roots[i2]})")
    for i1 in range(len(sg)):
        for i2 in range(i1): nzf.append(f"({sg[i1]}-{sg[i2]})")
    # members pairwise distinct: c2 != 1 or roots differ ... enforce via c2 - 1, c3 - 1, c2 - c3 only when shapes equal
    if b[0] == b[1]: nzf.append("(c2-1)")
    if b[0] == b[2]: nzf.append("(c3-1)")
    if b[1] == b[2]: nzf.append("(c2-c3)")
    L.append("poly nz = " + "*".join(nzf) + ";")
    L.append("int t0 = timer; ideal S = std(I + ideal(1 - TT*nz));")
    L.append(f'print("RESULT {name} D={D} G={sum(s)+3*t} dim=" + string(dim(S)) + " time=" + string(timer - t0));')
    L.append("quit;")
    open(os.path.join(outdir, name + ".sing"), "w").write("\n".join(L) + "\n")
print(len(jobs), "jobs")
