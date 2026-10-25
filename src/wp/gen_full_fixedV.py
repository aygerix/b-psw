"""Singular script: real-cubic Williams families (normal form, coprime members) for odd k >= 5, both halves, with
the universal part of the Laurent cofactor built in.
Normal form psi_i(y) = 1 + y a_i(y), a_i = 2 + z d_i + z^2 b_i (z = y + 1). "-" half: F - 1 = y prod(a_i) V,
V = y^(k-1) + R, deg R <= k - 4, with R(-1) = 2^(1-k) - 1 and R'(-1) = (k-1)(1 - 2^(1-k)) for every family
(laurent_fixed.gp; NOTEBOOK session 8). So R = A + B z + z^2 Q(y), deg Q <= k - 6 (free coefficients q_j).
Identity:  sum_{r<k} y^(r-1) e_r(a) = e_k(a) R(y).
"+" half: the same identity for the members (b_i, -b_i - d_i) with its own free Q+ (coefficients p_j).
Option --sat adds 1 - w * prod(b_i) * prod(alpha_i), alpha_i = 2 + d_i + b_i.
Option --lin: member k is the LINEAR member 2t + 1 (a_k = 2, i.e. b_k = d_k = 0; Lemma 6(b) allows at most one).
Then deg R <= k - 2 (not k - 4), so Q has k - 3 free coefficients per half (the two universal constraints at
y = -1 still hold: they use only a_i(-1) = 2).
Usage: python3 gen_full_fixedV.py k [--sat] [--lin]    writes hVfull<k>[_lin][_sat].sing
"""
import sys
from fractions import Fraction as Fr
k = int(sys.argv[1]); sat = "--sat" in sys.argv; lin = "--lin" in sys.argv
assert k % 2 == 1 and k >= (3 if lin else 5)
A = Fr(1, 2 ** (k - 1)) - 1
Bc = (k - 1) * (1 - Fr(1, 2 ** (k - 1)))
nq = (k - 3) if lin else (k - 5)             # free coefficients of Q (deg Q <= deg R - 2)
kc = k - 1 if lin else k                     # number of cubic members (variables b_i, d_i)
B = ",".join(f"b{i}" for i in range(1, kc + 1)); Dd = ",".join(f"d{i}" for i in range(1, kc + 1))
Q = ",".join(f"q{j}" for j in range(nq)); P = ",".join(f"p{j}" for j in range(nq))
extra = (("," + Q + "," + P) if nq else "") + (",w" if sat else "")
def rpoly(c):
    s = f"({A.numerator}/{A.denominator}) + ({Bc.numerator}/{Bc.denominator})*(y+1)"
    if nq: s += " + (y+1)^2*(" + " + ".join(f"{c}{j}*y^{j}" for j in range(nq)) + ")"
    return s
L = [f"ring R = 0, ({B},{Dd}{extra},y,X), dp;"]
for tag, dexpr, rp in (("m", lambda i: f"d{i}", rpoly("q")), ("p", lambda i: f"(-b{i}-d{i})", rpoly("p"))):
    L.append(f"poly G{tag} = 1;")
    for i in range(1, kc + 1):
        L.append(f"G{tag} = G{tag}*(1 + X*(2 + (y+1)*{dexpr(i)} + (y+1)^2*b{i}));")
    if lin:
        L.append(f"G{tag} = G{tag}*(1 + 2*X);")
    L.append(f"matrix E{tag} = coeffs(G{tag}, X); poly S{tag} = 0; int r{tag};")
    L.append(f"for (r{tag} = 1; r{tag} <= {k-1}; r{tag}++) {{ S{tag} = S{tag} + y^(r{tag}-1)*E{tag}[r{tag}+1,1]; }}")
    L.append(f"S{tag} = S{tag} - E{tag}[{k+1},1]*({rp});")
    L.append(f"matrix C{tag} = coeffs(S{tag}, y);")
L.append("ideal I; int j;")
L.append("for (j = 1; j <= nrows(Cm); j++) { I = I, Cm[j,1]; }")
L.append("for (j = 1; j <= nrows(Cp); j++) { I = I, Cp[j,1]; }")
if sat:
    pb = "*".join(f"b{i}" for i in range(1, kc + 1)); pa = "*".join(f"(2 + d{i} + b{i})" for i in range(1, kc + 1))
    L.append(f"I = I, 1 - w*{pb}*{pa};")
L.append(f"ring S2 = 0, ({B},{Dd}{extra}), dp;")
L.append("ideal J = imap(R, I); J = simplify(J, 2);")
L.append("option(redSB);")
L.append("ideal GB = std(J);")
L.append('if (GB[1] == 1) { "UNIT IDEAL: no solutions"; } else { "NOT UNIT: dim", dim(GB), "vdim", vdim(GB); }')
L.append("quit;")
name = f"hVfull{k}" + ("_lin" if lin else "") + ("_sat" if sat else "") + ".sing"
open(name, "w").write("\n".join(L) + "\n"); print("wrote", name)
