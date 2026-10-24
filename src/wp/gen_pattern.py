"""Generate Singular scripts for W' (k = 5, coprime normal form) with a given degree pattern.

Members f_j = sum_i c[j][i] t^i, deg f_j = d_j. Normalisation: the member of largest degree D is
f_1 = t^D + c_{1,D-2} t^(D-2) + ... (affine change of t kills the leading and next coefficient).
Conditions: e1 + e3 = 0 and e2 + e4 = 0 identically in t. Nonconstancy and exact degrees: the leading
coefficients are inverted with one Rabinowitsch variable y.
mode 'p': structure over GF(32003): dim/vdim of I', and whether H (the coefficients of e1 and e2) lies in
          radical(I') (H in the radical  <=>  every solution is monomial).
mode 'q': exact over Q with plain Buchberger (std): is I' the unit ideal; and for each h in H, is
          1 in I' + (1 - z h) (Rabinowitsch: h vanishes on V(I')).
Usage: python3 gen_pattern.py d1,d2,d3,d4,d5 mode [crit] > out.sing\ncrit e: H = coefficients of e1, e2 (equal exponents, f_j = zeta^j u); crit m (default for unequal degrees):\nH = non-leading coefficients of all members, i.e. every member is c_j t^(d_j) (generalised monomial with\nu = t, which is forced after normalisation whenever a linear member exists).
"""
import sys

def main():
    degs = sorted((int(x) for x in sys.argv[1].split(',')), reverse=True)
    mode = sys.argv[2]
    D = degs[0]
    names, lead, members, nonlead = [], [], [], []
    for j, d in enumerate(degs, 1):
        terms = []
        for i in range(d + 1):
            if j == 1 and i == D: terms.append(f"t^{i}"); continue
            if j == 1 and i == D - 1 and D >= 1: continue
            v = f"c{j}_{i}"; names.append(v); terms.append(f"{v}*t^{i}")
            if i == d: lead.append(v)
            else: nonlead.append(v)
        members.append(" + ".join(terms))
    crit = sys.argv[3] if len(sys.argv) > 3 else ('e' if len(set(degs)) == 1 else 'm')
    char = "32003" if mode == 'p' else "0"
    out = []
    out.append('LIB "primdec.lib"; LIB "modstd.lib";')
    out.append(f"ring R = {char}, ({','.join(names)},t,X), dp;")
    for j, m in enumerate(members, 1): out.append(f"poly f{j} = {m};")
    out.append("poly P = (X - f1)*(X - f2)*(X - f3)*(X - f4)*(X - f5);")
    out.append("matrix C = coeffs(P, X);")
    out.append("poly E1 = -C[5,1]; poly E2 = C[4,1]; poly E3 = -C[3,1]; poly E4 = C[2,1];")
    out.append("matrix G1 = coeffs(E1 + E3, t); matrix G2 = coeffs(E2 + E4, t);")
    out.append("matrix H1 = coeffs(E1, t); matrix H2 = coeffs(E2, t);")
    out.append("ideal I; ideal H; int i;")
    out.append("for (i = 1; i <= nrows(G1); i++) { I = I, G1[i,1]; }")
    out.append("for (i = 1; i <= nrows(G2); i++) { I = I, G2[i,1]; }")
    out.append("for (i = 1; i <= nrows(H1); i++) { H = H, H1[i,1]; }")
    out.append("for (i = 1; i <= nrows(H2); i++) { H = H, H2[i,1]; }")
    lp = "*".join(lead) if lead else "1"
    out.append(f"ring S = {char}, (z,y,{','.join(names)}), dp;")
    out.append("ideal I = simplify(imap(R, I), 2); ideal H = simplify(imap(R, H), 2);")
    if crit == 'm':
        out.append("H = " + (",".join(nonlead) if nonlead else "0") + ";   // criterion: every member is c_j t^(d_j)")
    out.append(f"ideal Ip = I, 1 - y*{lp};")
    out.append('"pattern ' + ",".join(map(str, degs)) + ': unknowns ' + str(len(names)) + ', equations", size(I);')
    out.append("int tt = timer;")
    if mode == 'p':
        out.append("ideal G = std(Ip);")
        out.append('"std(Ip): dim", dim(G) - 1, "(z is free, subtract 1), time", timer - tt;')
        out.append("ideal rG = radical(Ip);")
        out.append('"H in radical(Ip):", size(reduce(H, std(rG))) == 0, "time", timer - tt;')
        out.append("list mp = minAssGTZ(Ip);")
        out.append('"minimal primes:", size(mp);')
        out.append('for (i = 1; i <= size(mp); i++) { "prime", i, "dim", dim(std(mp[i])) - 1, "contains H:", size(reduce(H, std(mp[i]))) == 0; }')
    else:
        out.append("ideal G0 = std(Ip);")
        out.append('"exact std(Ip) over Q: unit ideal (no solutions):", (size(G0) == 1 && deg(G0[1]) == 0), "dim", dim(G0) - 1, "vdim", vdim(G0), "time", timer - tt;')
        out.append("int allok = 1; ideal T; ideal GT;")
        out.append("for (i = 1; i <= ncols(H); i++) {")
        out.append("  T = Ip, 1 - z*H[i]; GT = std(T);")
        out.append('  "h", i, ": 1 in I + (1 - z h):", (size(GT) == 1 && deg(GT[1]) == 0), "time", timer - tt;')
        out.append("  if (!(size(GT) == 1 && deg(GT[1]) == 0)) { allok = 0; }")
        out.append("}")
        out.append('"EVERY TEST POLYNOMIAL IN H VANISHES ON V(I) (criterion e: e1 = e2 = 0; criterion m: every member c_j t^(d_j)):", allok;')
    out.append("quit;")
    print("\n".join(out))

if __name__ == "__main__":
    main()
