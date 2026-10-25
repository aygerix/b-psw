\\ vrf line (adversarial verifier): INDEPENDENT exact check of the claimed k = 9 symmetric point of m9_symexact.
\\ No use of cubic_moments.gp / m9_lib.gp: V is computed by exact polynomial division
\\   V = quotient of (F - 1) by t * prod g_i,   g_i = (f_i - 1)/t,   f = 2t + 1 + t(t+1)(t + u)/v
\\ (normal form f = 2t+1 + t(t+1)(a t + c), a = 1/v, c = u/v).  The quotient of the division IS the polynomial part
\\ of the Laurent expansion at infinity.  Boundary equations: V(-1) = 2^(1-k), V'(-1) = (1-k) 2^(1-k), same for the
\\ reflected family u -> 1 - u.
t; r1; r2; r3; r4; Y; Y1; Y2; Y3; Y4;
fm(u, v) = 2*t + 1 + t*(t+1)*(t + u)/v;
gm(u, v) = 2 + (t+1)*(t + u)/v;
\\ L = vector of [u, v]; returns [V, R] with F - 1 = (t prod g) V + R
Vdiv(L, ex = 1) = {my(F = 1, G = t);
  for(i = 1, #L, my(f = fm(L[i][1], L[i][2]), g = gm(L[i][1], L[i][2]));
    if(ex && f - 1 - t*g != 0, error("g mismatch")); F *= f; G *= g);
  my(qr = divrem(F - 1, G, t)); if(poldegree(qr[1], t) != #L - 1, error("deg V")); if(ex && F - 1 - G*qr[1] - qr[2] != 0, error("division")); qr};
bnd(L, ex = 1) = {my(k = #L, qr = Vdiv(L, ex), V = qr[1], c0 = 2^(1-k), c1 = (1-k)*2^(1-k));
  [subst(V, t, -1) - c0, subst(deriv(V, t), t, -1) - c1, qr[2]]};
refl(L) = vector(#L, i, [1 - L[i][1], L[i][2]]);
toY(e) = {my(z = e); z = substpol(z, r1^2, Y1); z = substpol(z, r2^2, Y2); z = substpol(z, r3^2, Y3); z = substpol(z, r4^2, Y4); z};
tdeg(P) = poldegree(substvec(P, [Y1, Y2, Y3, Y4], [Y*Y1, Y*Y2, Y*Y3, Y*Y4]), Y);
\\ ---------------- claimed data (from the claim text and m9_symexact.out, pair choice [1,3]) ----------------
PV = [15129/25000, 671/12500, 6201/100000, 1193/20000]; vc = 80537/50000;
w2 = 2863/10^5; w4 = 1143/10^4;
D = 40873015522726052437652941104413136372114406461/426655372299127257812500000000000000000000000;
N = -3908730585275433747927744939801/12548687420562566406250000000;
Dn = -286107999639694964577/2362105867400012500;
\\ y1 = (N - sqrt D)/Dn  =>  minimal polynomial (Dn Y - N)^2 - D
qY = (Dn*Y - N)^2 - D;
y3of = 3640457795177541991529/138017997877500000000 - (3585828512/360831367)*Y;
{
  print("D square? ", issquare(D), "   numerator/denominator: ", issquare(numerator(D)), " ", issquare(denominator(D)), "  D*den^2 square? ", issquare(numerator(D)*denominator(D)));
  print("y2 = ", w2^2, "  y4 = ", w4^2);
  my(Lfam = [[1/2 + r1, PV[1]], [1/2 - r1, PV[1]], [1/2 + w2, PV[2]], [1/2 - w2, PV[2]],
             [1/2 + r3, PV[3]], [1/2 - r3, PV[3]], [1/2 + w4, PV[4]], [1/2 - w4, PV[4]], [1/2, vc]]);
  my(B = bnd(Lfam), Br = bnd(refl(Lfam)), E = [B[1], B[2], Br[1], Br[2]], res = vector(4));
  for(j = 1, 4, my(e = E[j]);
    if(subst(e, r1, -r1) != e || subst(e, r3, -r3) != e, error("odd powers"));
    e = substpol(e, r1^2, Y); e = substpol(e, r3^2, y3of);
    if(#setminus(Set(variables(e)), Set([Y])) > 0, error("leftover vars ", variables(e)));
    res[j] = divrem(e, qY, Y)[2]);
  print("INDEPENDENT exact residuals [E1, E2, E1r, E2r] mod minimal polynomial of y1: ", res);
  \\ numeric sanity at the real point
  default(realprecision, 200);
  my(y1 = (N - sqrt(D))/Dn, y3 = subst(y3of, Y, y1));
  print("y1 = ", precision(y1, 20), "  y3 = ", precision(y3, 20));
  my(Ln = [[1/2 + sqrt(y1), PV[1]], [1/2 - sqrt(y1), PV[1]], [1/2 + w2, PV[2]], [1/2 - w2, PV[2]],
           [1/2 + sqrt(y3), PV[3]], [1/2 - sqrt(y3), PV[3]], [1/2 + w4, PV[4]], [1/2 - w4, PV[4]], [1/2, vc]]);
  my(Bn = bnd(Ln, 0), Bnr = bnd(refl(Ln), 0));
  print("numeric residuals (200 digits): ", precision([Bn[1], Bn[2], Bnr[1], Bnr[2]], 5));
  print("full-family remainder R = (F-1) mod t prod g: max |coef| = ", precision(vecmax(abs(Vec(Bn[3]))), 6), " (nonzero expected: only the 4 boundary eqs are claimed)");
  \\ domain numerics: (|u - 1/2| + 1/2)^2 < 8 v  margins
  print("margins 8v - (|w|+1/2)^2: ", precision(vector(9, i, 8*Ln[i][2] - (abs(Ln[i][1] - 1/2) + 1/2)^2), 8));
  \\ ---------------- mechanism: symbolic y's in the symmetric slice ----------------
  my(Ls = [[1/2 + r1, PV[1]], [1/2 - r1, PV[1]], [1/2 + r2, PV[2]], [1/2 - r2, PV[2]],
           [1/2 + r3, PV[3]], [1/2 - r3, PV[3]], [1/2 + r4, PV[4]], [1/2 - r4, PV[4]], [1/2, vc]]);
  default(realprecision, 38);
  my(Bs = bnd(Ls), E1 = Bs[1], E2 = Bs[2]);
  E1 = toY(E1); E2 = toY(E2);
  \\ Z0 = E1, Z1 = T'(-1) = -8 E1 - E2  => 70 Z0 + 9 Z1 = -2 E1 - 9 E2
  my(Lc = -2*E1 - 9*E2);
  print("total degree in Y's: E1 ", tdeg(E1), ", E2 ", tdeg(E2), ", 70 Z0 + 9 Z1 = -2E1 - 9E2: ", tdeg(Lc));
  \\ re-derive y3(y1) and the quadratic from the claimed fixed data
  my(Lf = substvec(Lc, [Y2, Y4], [w2^2, w4^2]), Qf = substvec(E1, [Y2, Y4], [w2^2, w4^2]));
  my(a3 = polcoef(Lf, 1, Y3), y3der = -(Lf - a3*Y3)/a3);
  print("derived y3(y1) - claimed y3(y1): ", substpol(y3der, Y1, Y) - y3of);
  my(q1 = substpol(subst(Qf, Y3, y3der), Y1, Y));
  print("derived quadratic degree ", poldegree(q1, Y), "; proportional to claimed minimal poly? ", polcoef(q1, 2, Y)*qY - polcoef(qY, 2, Y)*q1 == 0);
}
quit;
