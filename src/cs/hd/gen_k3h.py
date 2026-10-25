# k = 3 Korselt families f_j = 1 + c_j A B_j (A monic deg a, common to all f_j - 1; B_j monic deg b_j, pairwise coprime and
# coprime to A: "non-degenerate"), c_1 = 1, translation fixed (A has no x^(a-1) term; if a = 0, B_1 has none).
# Lucas side EXACT with multiplicities: h_j = gcd(F+1, f_j+1) = gcd(f_j + 1, L_j), L_j = (c_j B_j)^2 - prod_{i != j}(c_j B_j - 2 c_i B_i)
# (A-reduction lemma).  Target: deg h_j >= g_j via a monic H_j of degree g_j dividing f_j + 1 and L_j.
# Output: dim of the saturated ideal (-1: no such family over C).
import itertools, sys, os
outdir = sys.argv[1]; specs = sys.argv[2]
os.makedirs(outdir, exist_ok=True)
def gen(a, b, g):
    D = 3*a + sum(b)
    name = f"h3_a{a}_b{''.join(map(str,b))}_g{''.join(map(str,g))}"
    V = []
    Acoef = [f"p{m}" for m in range(a)]            # A = x^a + p_{a-1} x^{a-1} + ... + p_0, with p_{a-1} = 0
    if a >= 1: Afree = Acoef[:a-1]
    Bc = [[f"b{j+1}_{m}" for m in range(b[j])] for j in range(3)]
    Hc = [[f"h{j+1}_{m}" for m in range(g[j])] for j in range(3)]
    fixed0 = set()
    if a >= 1: fixed0.add(f"p{a-1}")
    else: fixed0.add(f"b1_{b[0]-1}")
    allv = [v for v in Acoef + sum(Bc, []) + ["c2", "c3"] + sum(Hc, []) if v not in fixed0]
    def poly_of(coefs, deg):
        terms = [f"x^{deg}"] + [f"{c}*x^{m}" for m, c in enumerate(coefs) if c not in fixed0]
        return "+".join(terms)
    L = [f"// {name}: a={a} b={b} target gcd degrees {g} (sum {sum(g)}), D={D}",
         'LIB "elim.lib"; option(redSB);',
         'proc SAT(ideal I, poly g) { def S = sat(I, g); if (typeof(S) == "list") { ideal J = S[1]; return(J); } return(S); }',
         f"ring R = 0, (x, {', '.join(allv)}), dp;",
         f"poly A = {poly_of(Acoef, a) if a >= 1 else '1'};"]
    for j in range(3):
        L.append(f"poly B{j+1} = {poly_of(Bc[j], b[j]) if b[j] >= 1 else '1'};")
        L.append(f"poly H{j+1} = {poly_of(Hc[j], g[j]) if g[j] >= 1 else '1'};")
    L += ["poly c1 = 1;",
          "poly f1 = 1 + c1*A*B1; poly f2 = 1 + c2*A*B2; poly f3 = 1 + c3*A*B3;",
          "poly L1 = (c1*B1)^2 - (c1*B1 - 2*c2*B2)*(c1*B1 - 2*c3*B3);",
          "poly L2 = (c2*B2)^2 - (c2*B2 - 2*c1*B1)*(c2*B2 - 2*c3*B3);",
          "poly L3 = (c3*B3)^2 - (c3*B3 - 2*c1*B1)*(c3*B3 - 2*c2*B2);",
          "proc rem(poly P, poly M) { matrix C = coeffs(reduce(P, std(M)), x); ideal I = C; return(I); }",
          "// remainders w.r.t. the monic (in x) polynomial M: use division in x with parameters as coefficients",
          "proc remx(poly P, poly M) { int dM = deg(M, intvec(1, 0:(nvars(basering)-1))); poly Q = P; int dQ = deg(Q, intvec(1, 0:(nvars(basering)-1)));",
          "  while (dQ >= dM && Q != 0) { matrix CQ = coeffs(Q, x); poly t = CQ[dQ+1, 1]; Q = Q - t*x^(dQ-dM)*M; dQ = deg(Q, intvec(1, 0:(nvars(basering)-1))); }",
          "  matrix C = coeffs(Q, x); ideal I = C; return(I); }",
          "ideal I = 0;"]
    for j in range(3):
        o = [l for l in range(3) if l != j]
        if b[j] >= 1: L.append(f"I = I + remx(f{o[0]+1}*f{o[1]+1} - 1, B{j+1});")
        if g[j] >= 1:
            L.append(f"I = I + remx(f{j+1} + 1, H{j+1});")
            L.append(f"I = I + remx(L{j+1}, H{j+1});")
    nz = ["c2", "c3"]
    if a >= 1:
        for j in range(3):
            if b[j] >= 1: nz.append(f"resultant(A, B{j+1}, x)")
    for i in range(3):
        for j in range(i):
            if b[i] >= 1 and b[j] >= 1: nz.append(f"resultant(B{i+1}, B{j+1}, x)")
    if b[0] == 0 and b[1] == 0: nz.append("(c2-1)")
    if b[0] == 0 and b[2] == 0: nz.append("(c3-1)")
    if b[1] == 0 and b[2] == 0: nz.append("(c2-c3)")
    L += ["poly nz = " + "*".join(f"({t})" for t in nz) + ";",
          "ring S0 = 0, (" + ", ".join(allv) + ", TT), dp;",
          "ideal I = imap(R, I); poly nz = imap(R, nz);",
          "int t0 = timer; ideal S = std(I + ideal(1 - TT*nz));   // Rabinowitsch: V(S) = V(I) minus V(nz)",
          f'print("RESULT {name} D={D} G={sum(g)} dim=" + string(dim(S)) + " time=" + string(timer - t0));',
          "quit;"]
    open(os.path.join(outdir, name + ".sing"), "w").write("\n".join(L) + "\n")
    return name
for spec in specs.split(';'):
    a_, b_, g_ = spec.split('/')
    print(gen(int(a_), tuple(int(c) for c in b_), tuple(int(c) for c in g_)))
