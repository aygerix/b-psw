\\ m9 line, k = 9 (nine cubics): an EXPLICIT ALGEBRAIC point (in a real quadratic field) of the four boundary
\\ equations inside the admissible domain, in the reflection-symmetric slice:
\\   4 reflected pairs (1/2 +- w_j, v_j) + one self-dual member (1/2, v_c).  The multiset is reflection invariant,
\\   so E1r = E1, E2r = E2 and only Z0+ = Z1+ = 0 remain (odd moments vanish).  Z0+, Z1+ involve the pair data only
\\   through y_j = w_j^2 via M1_2, M2_2 (linear) and M1_4 (quadratic); 70 Z0+ + 9 Z1+ is free of M1_4.
\\ All v's and w_3, w_4 rational; (y_1, y_2) solve a linear + a quadratic equation -> y_1 in Q(sqrt(D)).
\\ Verification (EXACT): the brute-force eqs(9, 0) of cubic_moments.gp vanish after u = 1/2 +- r_j, r_j^2 -> y_j and
\\ reduction modulo the minimal polynomial of y_1; the domain inequalities are decided exactly in Q(sqrt D).
read("m9_lib.gp");
Zs(m, lin) = {my(r = mform(m, lin), Z = r[1], Zr = [reflect(Z[1]), reflect(Z[2])]);
  [(Z[1] + Zr[1])/2, (Z[2] + Zr[2])/2, (Z[1] - Zr[1])/2, (Z[2] - Zr[2])/2]};
\\ substitute a symmetric configuration: pairs PV (v_j), PY (y_j), centre vc
symsub(P, PV, PY, vc) = {my(R = P); foreach(syms(P), z, my(st = strsplit(Str(z), "_"), a = eval(strjoin(Vec(st[1])[2..#Vec(st[1])])), b = eval(st[2]));
  R = subst(R, z, if(b % 2, 0, sum(j = 1, #PV, 2 * PV[j]^a * PY[j]^(b/2)) + if(b == 0, vc^a, 0)))); R};
\\ exact sign of p + q sqrt(D), D > 0 non-square
sgnq(p, q, D) = {if(q == 0, return(sign(p))); if(p == 0, return(sign(q)));
  if(sign(p) == sign(q), return(sign(p))); sign(p) * sign(p^2 - q^2 * D)};
{
  my(ZZ = Zs(9, 0), Z0 = ZZ[1], Z1 = ZZ[2], den = 10^5);
  \\ numerical symmetric solution (certified by m9_kcert.gp): pairs (w, v) and centre
  my(Wn = [1.6286121386777186, 0.0286312273614805, 0.1347922858160165, 0.1143013292143569],
     Vn = [0.605156204392249, 0.053680237637796, 0.0620060685662747, 0.059648427249592], vcn = 1.6107368111675027);
  my(PV = vector(4, j, round(Vn[j] * den) / den), vc = round(vcn * den) / den, w3 = round(Wn[3] * den) / den, w4 = round(Wn[4] * den) / den);
  \\ try every choice of the two unknown pairs
  forsubset([4, 2], ss, my(i1 = ss[1], i2 = ss[2], PY = vector(4), L, Qd, a1, a2, b0, Y2of, q1, D, cand);
    my(fixw = setminus([1, 2, 3, 4], Set([i1, i2])));
    for(j = 1, 4, PY[j] = round(Wn[j] * den)^2 / den^2);
    PY[i1] = Y1; PY[i2] = Y2;
    L = simplify(symsub(70 * Z0 + 9 * Z1, PV, PY, vc));
    Qd = simplify(symsub(Z0, PV, PY, vc));
    if(poldegree(L, Y1) > 1 || poldegree(L, Y2) > 1, error("L not linear"));
    a2 = polcoef(L, 1, Y2); if(a2 == 0, next);
    Y2of = simplify(-(L - a2 * Y2) / a2);                 \\ Y2 as affine function of Y1
    q1 = simplify(subst(Qd, Y2, Y2of));                   \\ quadratic in Y1
    if(poldegree(q1, Y1) != 2, next);
    D = poldisc(q1);
    if(D <= 0 || issquare(D), next);
    my(A2 = polcoef(q1, 2, Y1), A1 = polcoef(q1, 1, Y1), y1n = Wn[i1]^2, y2n = Wn[i2]^2, best = 0, bestd = oo);
    \\ roots Y1 = (-A1 +- sqrt(D)) / (2 A2): as p + q sqrt(D)
    foreach([1, -1], e, my(p = -A1 / (2 * A2), q = e / (2 * A2), y1 = p + q * sqrt(1.0 * D), y2 = subst(Y2of, Y1, y1));
      my(dd = abs(y1 - y1n) + abs(y2 - y2n)); if(dd < bestd, bestd = dd; best = [p, q]));
    my(p = best[1], q = best[2]);
    \\ Y2 = c0 + c1 Y1 = (c0 + c1 p) + c1 q sqrt(D)
    my(c1 = polcoef(Y2of, 1, Y1), c0 = polcoef(Y2of, 0, Y1), p2 = c0 + c1 * p, q2 = c1 * q);
    print("==== unknown pairs ", [i1, i2], ": D = ", D, " (", #Str(D), " digits); y1 = ", precision(p + q * sqrt(1.0*D), 25), ", y2 = ", precision(p2 + q2 * sqrt(1.0*D), 25), "  (distance to numerics ", bestd, ")");
    \\ exact domain checks: y > 0 and sqrt(y) + 1/2 < sqrt(8 v)  <=>  8v - 1/4 - y > 0  and  y < (8v - 1/4 - y)^2
    my(okd = 1, YP = vector(4, j, if(j == i1, [p, q], if(j == i2, [p2, q2], [PY[j], 0]))));
    for(j = 1, 4, my(yp = YP[j][1], yq = YP[j][2], cc = 8 * PV[j] - 1/4);
      if(sgnq(yp, yq, D) <= 0, okd = 0);
      if(sgnq(cc - yp, -yq, D) <= 0, okd = 0);
      \\ (cc - y)^2 - y = (cc - yp)^2 + yq^2 D - yp + (-2 (cc - yp) yq - yq) sqrt(D)
      if(sgnq((cc - yp)^2 + yq^2 * D - yp, -2 * (cc - yp) * yq - yq, D) <= 0, okd = 0));
    if(!(8 * vc > 1/4), okd = 0);
    \\ distinct members: pairs have distinct v's here (checked) and y > 0; centre distinct from pairs (w = 0 vs y > 0)
    if(#Set(PV) < 4, okd = 0);
    print("  exact domain inequalities and distinctness: ", if(okd, "HOLD", "FAIL"));
    \\ exact verification with the brute-force boundary equations
    my(E = eqs(9, 0), res = vector(4));
    for(r = 1, 4, my(e = E[r]);
      for(j = 1, 4, my(wj = if(j == i1, r1, if(j == i2, r2, round(Wn[j] * den) / den)));
        e = subst(subst(e, eval(Str("u", j)), 1/2 + wj), eval(Str("v", j)), PV[j]);
        e = subst(subst(e, eval(Str("u", j + 4)), 1/2 - wj), eval(Str("v", j + 4)), PV[j]));
      e = subst(subst(e, u9, 1/2), v9, vc);
      \\ only even powers of r1, r2 may survive
      if(substvec(e, [r1, r2], [-r1, r2]) != e || substvec(e, [r1, r2], [r1, -r2]) != e, error("odd powers survive"));
      e = substpol(substpol(e, r1^2, Y1), r2^2, Y2);
      e = subst(e, Y2, Y2of);
      e = divrem(e, q1, Y1)[2];
      res[r] = e);
    print("  EXACT residuals of the four boundary equations (brute force, reduced mod the minimal polynomial of y1): ", res);
    if(okd && res == [0, 0, 0, 0],
      print("  => EXPLICIT ALGEBRAIC FEASIBLE POINT.  v (pairs) = ", PV, ", v_c = ", vc, ", fixed y = ", vector(#fixw, j, [fixw[j], PY[fixw[j]]]));
      print("     y_", i1, " = (", -A1, " + (", 2 * q * A2, ") sqrt(D)) / (", 2 * A2, "),  y_", i2, " = ", c0, " + (", c1, ") y_", i1);
      print("     D = ", D)));
}
