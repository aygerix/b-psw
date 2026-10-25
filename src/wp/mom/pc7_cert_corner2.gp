\\ pc7_cert_corner2.gp -- leading ep-coefficient of EVERY coefficient of D1 (and D5) at the corner, with w unscaled:
\\ r1 = -1 + ep a, r2 = 1 - ep b.  D_i = kap + sum_m A_m p_m + sum_n c_n h_n + lam3 e3.
read("pc7_cert_shape.gp");
ep; a; b;
PREC = 16;
{
  my(K = SYM0, A = vector(8, i, SYM1[i] - SYM0[i]));
  my(C = matrix(8, 5, i, n, polcoef(polcoef(SYM2[i] - SYM0[i] - A[i] - subst(A[i], w, v), 1, w), n, v)));
  my(e = ep + O(ep^PREC), R1 = -1 + e*a, R2 = 1 - e*b, Dl = (1 - R1^2)*(1 - R2^2));
  my(SUB = [R1 + R2, R1*R2, 2*(R1 + R2)/Dl, 2*(1 - (R1 + R2)^2 + R1*R2)/Dl]);
  my(sb(P) = substvec(P, [e1, e2, lam, mu], SUB));
  my(ld(P) = if(P == 0, "0", my(Q = sb(P), nu = valuation(Q, ep)); my(L0 = polcoef(Q, nu, ep), nn = numerator(L0), dd = denominator(L0)); Str("ep^", nu, " * (", content(nn)/content(dd), ") * ", factor(nn / content(nn)), " / ", factor(dd / content(dd)))));
  foreach([1, 5], i,
    print("D", i, ":  kap ~ ", ld(K[i]));
    for(m = 1, 9, print("   A", m, " ~ ", ld(polcoef(A[i], m, w))));
    for(n = 2, 6, print("   c", n, " ~ ", ld(C[i, n - 1]))));
}
