# AUDIT generator (new formulation, independent of the line's H-variable scripts).
# k = 3, f_1 = 1 + A*B1, f_2 = 1 + A*B2, f_3 = 1 + A*B3, A monic deg a with no x^(a-1) term (translation), c_1 = 1 (scaling),
# B1 monic (or 1), B2, B3 with free leading coefficients (absorbing c_2, c_3).  Korselt: B_j | (F-1)/A via resultant.
# Lucas: deg h_j >= 1 via Res_x(f_j + 1, F + 1) (direct, mode 'R') or Res_x(f_j+1, L_j) (A-reduction, mode 'L');
# deg h_j = full deg L_j via pseudo-remainder of f_j+1 by L_j ('Ldiv'), or f_j+1 | F+1 via pseudo-remainder ('Fdiv').
# Non-degeneracy saturated by Rabinowitsch: lc's, Res(A, B_j), Res(B_i, B_j).
import sys
def script(name, a, b, conds):
    v = []
    if a >= 1:
        Ac = [f"p{m}" for m in range(a-1)]   # x^(a-1) coefficient fixed to 0
        v += Ac
        A = "+".join([f"x^{a}"] + [f"{c}*x^{m}" for m, c in enumerate(Ac)])
    else:
        A = "1"
    Bs = []
    for j in range(3):
        if b[j] == 0:
            if j == 0: Bs.append("1")
            else: v.append(f"k{j+1}"); Bs.append(f"k{j+1}")
        else:
            cs = [f"b{j+1}_{m}" for m in range(b[j]+1)]
            if j == 0: cs[-1] = "1"
            v += [c for c in cs if c != "1"]
            Bs.append("+".join(f"{c}*x^{m}" for m, c in enumerate(cs)))
    L = [f"// {name}", 'option(noredefine);',
         f"ring R = 0, (x, {', '.join(v)}, TT), dp;",
         f"poly A = {A};", f"poly B1 = {Bs[0]}; poly B2 = {Bs[1]}; poly B3 = {Bs[2]};",
         "poly f1 = 1 + A*B1; poly f2 = 1 + A*B2; poly f3 = 1 + A*B3; poly F = f1*f2*f3;",
         "poly Q = division(F - 1, A)[1][1,1];",
         "if (Q*A != F - 1) { print(\"division error\"); quit; }",
         "list Bl = B1, B2, B3; list fl = f1, f2, f3;",
         "poly L1 = 2*B1*(B2+B3) - 4*B2*B3; poly L2 = 2*B2*(B1+B3) - 4*B1*B3; poly L3 = 2*B3*(B1+B2) - 4*B1*B2;",
         "list Ll = L1, L2, L3;",
         "proc dx(poly P) { return(deg(P, intvec(1, 0:(nvars(basering)-1)))); }",
         "proc lcx(poly P) { matrix C = coeffs(P, x); return(C[nrows(C), 1]); }",
         "proc pdivrem(poly P, poly M) { int dM = dx(M); poly l = lcx(M); poly Qr = P; while (Qr != 0 && dx(Qr) >= dM) { Qr = l*Qr - lcx(Qr)*x^(dx(Qr)-dM)*M; } matrix C = coeffs(Qr, x); ideal I = C; return(I); }",
         "ideal I = 0;"]
    for j in range(3):
        if b[j] >= 1: L.append(f"I = I + resultant(Bl[{j+1}], Q, x);")
    for (j, kind) in conds:
        if kind == 'R': L.append(f"I = I + resultant(fl[{j}] + 1, F + 1, x);")
        elif kind == 'L': L.append(f"I = I + resultant(fl[{j}] + 1, Ll[{j}], x);")
        elif kind == 'Ldiv': L.append(f"I = I + pdivrem(fl[{j}] + 1, Ll[{j}]);")
        elif kind == 'Fdiv': L.append(f"I = I + pdivrem(F + 1, fl[{j}] + 1);")
    nz = ["lcx(f1)", "lcx(f2)", "lcx(f3)"]
    if a >= 1:
        for j in range(3):
            if b[j] >= 1: nz.append(f"resultant(A, Bl[{j+1}], x)")
    for i in range(3):
        for j in range(i):
            if b[i] >= 1 and b[j] >= 1: nz.append(f"resultant(Bl[{i+1}], Bl[{j+1}], x)")
    for (j, kind) in conds:
        if kind in ('L', 'Ldiv'): nz.append(f"lcx(Ll[{j}])")
    L += ["poly nz = " + "*".join(f"({t})" for t in nz) + ";",
          f"ring S0 = 0, ({', '.join(v)}, TT), dp;", "ideal I = imap(R, I); poly nz = imap(R, nz);",
          "int t0 = timer; ideal S = std(I + ideal(1 - TT*nz));",
          f'print("RESULT {name} dim=" + string(dim(S)) + " time=" + string(timer - t0));', "quit;"]
    open(name + ".sing", "w").write("\n".join(L) + "\n")
jobs = {
 # (1;0,1,1): D = 5, G = 3 forces deg h_j = 1 for all j (bounds min(d_j, D - 2d_j) = 1,1,1)
 "x_a1_011_R":   (1, (0,1,1), [(1,'R'), (2,'R'), (3,'R')]),
 # control: drop member 3's condition -> must be nonempty (dim >= 0) if the stratum with G = 2 exists
 "x_a1_011_ctl": (1, (0,1,1), [(1,'R'), (2,'R')]),
 "x_a1_011_K":   (1, (0,1,1), []),
 "x_a3_011_K":   (3, (0,1,1), []),
 "x_a1_111_K":   (1, (1,1,1), []),
 # (2;0,1,1): D = 8, G >= 5: patterns (2,2,1), (2,1,2), (1,2,2) with deg L_j = 2 exactly
 "x_a2_011_221": (2, (0,1,1), [(1,'Ldiv'), (2,'Ldiv'), (3,'L')]),
 "x_a2_011_212": (2, (0,1,1), [(1,'Ldiv'), (2,'L'), (3,'Ldiv')]),
 "x_a2_011_122": (2, (0,1,1), [(1,'L'), (2,'Ldiv'), (3,'Ldiv')]),
 "x_a2_011_ctl": (2, (0,1,1), [(1,'Ldiv'), (2,'L'), (3,'L')]),
 # (3;0,1,1): D = 11, G = 6 forces deg h_j = 2 = deg L_j for all j
 "x_a3_011_222": (3, (0,1,1), [(1,'Ldiv'), (2,'Ldiv'), (3,'Ldiv')]),
 # (1;1,1,1): D = 6, G >= 4: (2,2,*) with full members 1,2 or (2,1,1) with full member 1 (scaling WLOG), direct form
 "x_a1_111_22":  (1, (1,1,1), [(1,'Fdiv'), (2,'Fdiv')]),
 "x_a1_111_211": (1, (1,1,1), [(1,'Fdiv'), (2,'R'), (3,'R')]),
}
for n, (a, b, c) in jobs.items(): script(n, a, b, c)
print(" ".join(jobs))
