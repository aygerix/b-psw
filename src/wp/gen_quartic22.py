"""Singular script: REAL quartic Williams families with sharing pattern (2,2): every member has exactly two real "-"
points (normalised to t = -1, 1) and two real "+" points u1, u2 (roots of H = t^2 - sg t + rh), shared by all members
(Lemma 6); all other special points non-real. Then f_j = f_0 + c_j G H (a PENCIL), G = t^2 - 1, where
f_0 = 1 + G (al t + be) and f_0 + 1 = H (ga t + de) (polynomial identity: 4 equations).
Reduction (PROVED, elementary): modulo K-_j = al t + be + c_j H (member j's non-real "-" points),
  c_j P == -G(al t + be) = -(f_0 - 1), so f_i == 1 + (1 - c_i/c_j) X with X = f_0 - 1; the condition F == 1 becomes
  prod_{i != j} (c_j + (c_j - c_i) X) - c_j^(k-1) == 0 mod K-_j.
Modulo K+_j = ga t + de + c_j G ("+" points): f_i == -1 + (1 - c_i/c_j) Y, Y = f_0 + 1, and F == -1 becomes (k odd)
  prod_{i != j} (c_j - (c_j - c_i) Y) - c_j^(k-1) == 0 mod K+_j.
Since X(w) = -c_j P(w) != 0 (P vanishes only at the real points), the factor X (resp. Y) of R_j is divided out.
Each "== 0 mod K" is a pseudo-remainder (in t) with 2 coefficients. Saturation: c_j, c_i - c_j, disc H, H(1) H(-1).
Usage: python3 gen_quartic22.py k     writes q22_k<k>.sing
"""
import sys
k = int(sys.argv[1])
C = ",".join(f"c{i}" for i in range(1, k + 1))
L = [f"ring R = 0, (t,Z,al,be,ga,de,sg,rh,{C},w), (dp(2),dp);",
     "proc degt(poly f) { if (f == 0) { return(-1); } return(deg(f, intvec(1,0:(nvars(basering)-1)))); }",
     "proc lct(poly f) { int d = degt(f); matrix M = coeffs(f, t); return(M[d+1,1]); }",
     "proc prem(poly f, poly g) { int dg = degt(g); poly lc = lct(g); int df = degt(f);",
     "  while (df >= dg) { f = lc*f - lct(f)*t^(df-dg)*g; df = degt(f); } return(f); }",
     "poly G = t^2 - 1; poly H = t^2 - sg*t + rh;",
     "poly f0 = 1 + G*(al*t + be); poly X = f0 - 1; poly Y = f0 + 1;",
     "ideal I = coeffs(Y - H*(ga*t + de), t);",
     "poly Km; poly Kp; poly Pm; poly Pp; matrix Mc;"]
for j in range(1, k + 1):
    others = [i for i in range(1, k + 1) if i != j]
    # R_j(Z)/Z with Z = X resp. Y: X(w) = -c_j P(w) != 0 at non-real points, so the factor Z is dropped
    pm = "subst((" + "*".join(f"(c{j} + (c{j} - c{i})*Z)" for i in others) + f" - c{j}^{k-1})/Z, Z, X)"
    pp = "subst((" + "*".join(f"(c{j} - (c{j} - c{i})*Z)" for i in others) + f" - c{j}^{k-1})/Z, Z, Y)"
    L.append(f"Km = al*t + be + c{j}*H; Kp = ga*t + de + c{j}*G;")
    L.append(f"Pm = prem({pm}, Km); Pp = prem({pp}, Kp);")
    L.append("I = I, coeffs(Pm, t), coeffs(Pp, t);")
dist = "*".join(f"(c{i} - c{j})" for i in range(1, k + 1) for j in range(i + 1, k + 1))
pc = "*".join(f"c{i}" for i in range(1, k + 1))
L.append(f"I = I, 1 - w*{pc}*{dist}*(sg^2 - 4*rh)*(1 - sg + rh)*(1 + sg + rh);")
L.append(f"ring S = 0, (al,be,ga,de,sg,rh,{C},w), dp; ideal J = imap(R, I); J = simplify(J, 2);")
L.append('"generators:", size(J);')
L.append("option(redSB); ideal GB = std(J);")
L.append('if (GB[1] == 1) { "UNIT IDEAL: no (2,2) quartic family"; } else { "NOT UNIT: dim", dim(GB), "vdim", vdim(GB); }')
L.append("quit;")
open(f"q22_k{k}.sing", "w").write("\n".join(L) + "\n"); print("wrote", f"q22_k{k}.sing")
