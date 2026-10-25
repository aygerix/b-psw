"""Singular script for the k = 5 real-cubic "-" half in Laurent form with the cofactor V FIXED
(V = y^4 + 15/4 y + 45/16, forced for every family by the expansion at the shared "+" point; see
laurent_fixed.gp). Members a_i = 2 + z d_i + z^2 b_i, z = y + 1; psi_i = 1 + y a_i.
Identity: sum_{r<5} y^(r-1) e_r(a) = e_5(a) (15/4 y + 45/16).  Unknowns b_i, d_i (10); square system.
Option --full adds the "+" half: by t -> -1 - t, f -> -f it is the same identity for the members
(b_i, -b_i - d_i). Writes hVfull5.sing (with --sat: hVfull5_sat.sing).
Options: --sat adds 1 - w * prod(b_i) * prod(alpha_i) with alpha_i = a_i(0) = 2 + d_i + b_i (members of
degree 3 and t = 0 a simple special point). Writes hV5.sing (or hV5_sat.sing).
Usage: python3 gen_half_fixedV.py [--sat]
"""
import sys
k = 5
sat = "--sat" in sys.argv
full = "--full" in sys.argv
B = ",".join(f"b{i}" for i in range(1, k + 1)); Dd = ",".join(f"d{i}" for i in range(1, k + 1))
L = []
L.append(f"ring R = 0, ({B},{Dd},{'w,' if sat else ''}y,X), dp;")
L.append("poly G = 1;")
for i in range(1, k + 1):
    L.append(f"poly a{i} = 2 + (y+1)*d{i} + (y+1)^2*b{i}; G = G*(1 + X*a{i});")
L.append("matrix E = coeffs(G, X);")            # E[r+1,1] = e_r
L.append("poly S = 0; int r;")
L.append(f"for (r = 1; r <= {k-1}; r++) {{ S = S + y^(r-1)*E[r+1,1]; }}")
L.append(f"S = S - E[{k+1},1]*(15/4*y + 45/16);")
L.append("matrix C = coeffs(S, y);")
L.append("ideal I; int j; for (j = 1; j <= nrows(C); j++) { I = I, C[j,1]; }")
if full:
    L.append("poly Gp = 1;")
    for i in range(1, k + 1):
        L.append(f"poly ap{i} = 2 + (y+1)*(-b{i}-d{i}) + (y+1)^2*b{i}; Gp = Gp*(1 + X*ap{i});")
    L.append("matrix Ep = coeffs(Gp, X);")
    L.append("poly Sp = 0;")
    L.append(f"for (r = 1; r <= {k-1}; r++) {{ Sp = Sp + y^(r-1)*Ep[r+1,1]; }}")
    L.append(f"Sp = Sp - Ep[{k+1},1]*(15/4*y + 45/16);")
    L.append("matrix Cp = coeffs(Sp, y);")
    L.append("for (j = 1; j <= nrows(Cp); j++) { I = I, Cp[j,1]; }")
if sat:
    pb = "*".join(f"b{i}" for i in range(1, k + 1)); pa = "*".join(f"(2 + d{i} + b{i})" for i in range(1, k + 1))
    L.append(f"I = I, 1 - w*{pb}*{pa};")
L.append(f"ring S2 = 0, ({B},{Dd}{',w' if sat else ''}), dp;")
L.append("ideal J = imap(R, I); J = simplify(J, 2);")
L.append("option(redSB); option(prot);")
L.append("ideal GB = std(J);")
L.append('if (GB[1] == 1) { "UNIT IDEAL: no solutions"; } else { "NOT UNIT: dim", dim(GB), "vdim", vdim(GB); }')
L.append("quit;")
name = ("hVfull5" if full else "hV5") + ("_sat" if sat else "") + ".sing"
open(name, "w").write("\n".join(L) + "\n")
