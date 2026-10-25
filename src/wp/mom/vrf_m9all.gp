\\ vrf line: independent exact check of EVERY point printed in m9_symexact.out (pair choices with formulas).
read("vrf_m9all_data.gp");
PV = [15129/25000, 671/12500, 6201/100000, 1193/20000]; vc = 80537/50000;
RV = [r1, r2, r3, r4];
{
  for(ci = 1, #CASES, my(C = CASES[ci], i1 = C[1][1], i2 = C[1][2], fx = C[2], D = C[3], N = C[4], e = C[5], Dn = C[6], c0 = C[7], c1 = C[8]);
    my(qY = (Dn*Y - N)^2 - D, yof = vector(4), W = vector(4), L, res = vector(4), ok = 1);
    \\ y_i1 = (N + e sqrt D)/Dn  <=>  root of qY ; y_i2 = c0 + c1 y_i1 ; others fixed rational
    yof[i1] = Y; yof[i2] = c0 + c1*Y; for(j = 1, 2, yof[fx[j][1]] = fx[j][2]);
    for(j = 1, 4, W[j] = if(j == i1 || j == i2, RV[j], my(s = sqrtint(numerator(yof[j])), d = sqrtint(denominator(yof[j]))); if(s^2 != numerator(yof[j]) || d^2 != denominator(yof[j]), error("fixed y not square")); s/d));
    L = concat(vector(4, j, [[1/2 + W[j], PV[j]], [1/2 - W[j], PV[j]]]));
    L = concat(L, [[1/2, vc]]);
    my(B = concat(bnd(L), bnd(refl(L))));
    for(j = 1, 4, my(z = B[j]);
      foreach([i1, i2], ii, if(subst(z, RV[ii], -RV[ii]) != z, error("odd")); z = substpol(z, RV[ii]^2, yof[ii]));
      if(#setminus(Set(variables(z)), Set([Y])) > 0, error("vars"));
      res[j] = divrem(z, qY, Y)[2]);
    \\ exact domain check via rational enclosures
    my(sD = sqenc(D), yI = vector(4));
    for(j = 1, 4, if(j == i1, my(a = (N + e*sD[1])/Dn, b = (N + e*sD[2])/Dn); yI[j] = [min(a, b), max(a, b)],
                  if(j == i2, 0, yI[j] = [yof[j], yof[j]])));
    my(a = c0 + c1*yI[i1][1], b = c0 + c1*yI[i1][2]); yI[i2] = [min(a, b), max(a, b)];
    my(marg = vector(4));
    for(j = 1, 4, if(yI[j][1] <= 0, ok = 0; marg[j] = "y<=0", my(sh = sqenc(yI[j][2])[2]); marg[j] = 8*PV[j] - (sh + 1/2)^2; if(marg[j] <= 0, ok = 0); marg[j] = precision(1.0*marg[j], 6)));
    if(PV[1] <= 1/32 || PV[2] <= 1/32 || PV[3] <= 1/32 || PV[4] <= 1/32 || #Set(concat(PV, vc)) < 5, ok = 0);
    print("case ", C[1], ": exact residuals ", res, "; D square? ", issquare(D), "; y enclosures ", vector(4, j, precision(1.0*yI[j][1], 12)), "; certified margins ", marg, "; DOMAIN ", if(ok, "OK", "FAIL")));
}
quit;
