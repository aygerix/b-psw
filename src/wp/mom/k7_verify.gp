\\ k = 7 real cubic Williams families: exact verification of the advisor's moment contradiction.
\\ (1) the four boundary equations equal explicit moment polynomials (identity check, exact);
\\ (2) forced M = sum v^2 and K = sum v (u - 1/2)^2 as functions of S = sum v;
\\ (3) the separator: with v_i > 1/32 (reality: c^2 < 8a, (a - c)^2 < 8a => a < 32) and m cubic members,
\\     M < (S - (m-1)/32)^2 + (m-1)/1024; K >= 0 (v_i > 0). Both functions are monotone on the relevant range.
read("cubic_moments.gp");
mexpr(m, cv, cpair, cst) = {my(e = cst);
  for(i = 1, m, my(u = eval(Str("u", i)), v = eval(Str("v", i)));
    e += cv[1]*v + cv[2]*v*u + cv[3]*v*u^2 + cv[4]*v*u^3 + cv[5]*v^2 + cv[6]*v^2*u);
  for(i = 1, m, for(j = i + 1, m, e += cpair * eval(Str("v", i)) * eval(Str("v", j)))); e};
{
  my(T7 = [[-4,-3,-2,-1,6,4], 1, 63/64; [6,3,1,0,-2,0], 0, -189/32; [-10,10,-5,1,10,-4], 1, 63/64; [10,-5,1,0,-2,0], 0, -189/32]);
  my(T6 = [[-5/2,-2,-3/2,-1,5,4], 1, 31/64; [9/2,5/2,1,0,-2,0], 0, -109/32; [-7,8,-9/2,1,9,-4], 1, 31/64; [8,-9/2,1,0,-2,0], 0, -109/32]);
  foreach([[7, 0, T7], [6, 1, T6]], c, my(m = c[1], lin = c[2], T = c[3], E = eqs(m, lin));
    for(j = 1, 4, my(d = E[j] - mexpr(m, T[j,1], T[j,2], T[j,3])); print("m=", m, " lin=", lin, " E", j, " identity residual: ", d)));
  \\ symbolic moments: S, P1 = sum vu, P2 = sum vu^2, P3, M, N = sum v^2 u; pair sum = (S^2 - M)/2
  my(red(T, j) = my(r = T[j,1]); r[1]*S + r[2]*P1 + r[3]*P2 + r[4]*P3 + r[5]*M + r[6]*N + T[j,2]*(S^2 - M)/2 + T[j,3]);
  foreach([[7, T7], [6, T6]], c, my(m = c[1], T = c[2], e = vector(4, j, red(T, j)));
    \\ e2 - e4 is linear in P1, S: solve P1; then e2 for P2; then e1 + e3 (P3, N cancel) for M
    my(p1 = polcoef(e[2] - e[4], 0, P1) / -polcoef(e[2] - e[4], 1, P1));
    my(s2 = subst(e[2], P1, p1), p2 = polcoef(s2, 0, P2) / -polcoef(s2, 1, P2));
    my(s13 = subst(subst(e[1] + e[3], P1, p1), P2, p2));
    if(polcoef(s13, 1, P3) != 0 || polcoef(s13, 1, N) != 0, error("P3/N do not cancel"));
    my(Mf = polcoef(s13, 0, M) / -polcoef(s13, 1, M), Kf = subst(p2, M, Mf) - p1 + S/4);
    print("m=", m, ": P1 = ", p1, ";  M = ", Mf, ";  K = ", Kf);
    \\ separator
    my(r = if(m == 7, 181/200, 669/1000), A = subst(Kf, S, r), B = (r - (m-1)/32)^2 + (m-1)/1024 - subst(Mf, S, r));
    print("   r = ", r, ":  K(r) = ", A, "  [Mmax - M](r) = ", B);
    \\ monotonicity on (m/32, oo): K' < 0 and (Mmax - M)' > 0 there
    my(dK = deriv(Kf, S), dB = deriv((S - (m-1)/32)^2 + (m-1)/1024 - Mf, S));
    print("   K'(S) = ", dK, "  at S = m/32: ", subst(dK, S, m/32), ";  (Mmax-M)'(S) = ", dB, "  at S = m/32: ", subst(dB, S, m/32));
    \\ real roots of the two boundary polynomials (for the record)
    print("   root of K: ", polrootsreal(simplify(subst(Kf, S, xx))), "  root of Mmax - M: ", polrootsreal(simplify(subst((S - (m-1)/32)^2 + (m-1)/1024 - Mf, S, xx)))));
}
