\\ k = 7 families: m cubics, lin in {0,1} linear member 2t+1, and 7 - m - lin "invisible" members of odd degree >= 7
\\ (under (H1), (H2) of Theorem 11'), which contribute nothing to the polynomial part of V (degree 6).
\\ For each (m, lin): the four boundary equations, forced moments, and the k = 7 separator test.
read("cubic_moments.gp");
VpolyK(K, m, lin, U, W) = {
  my(G = 1);
  for(i = 1, m, G = trunc(G * (1 + xser(U[i], W[i], K)), K));
  if(lin, G = trunc(G * (1 + s/2), K));
  sum(j = 0, K - 1, polcoef(G, j, s) * t^(K-1-j))};
eqsK(K, m, lin) = {
  my(U = vector(m, i, eval(Str("u", i))), W = vector(m, i, eval(Str("v", i))), Ur = vector(m, i, 1 - U[i]),
     V = VpolyK(K, m, lin, U, W), Vr = VpolyK(K, m, lin, Ur, W), c0 = 2^(1-K), c1 = (1-K)*2^(1-K));
  [subst(V, t, -1) - c0, subst(deriv(V, t), t, -1) - c1, subst(Vr, t, -1) - c0, subst(deriv(Vr, t), t, -1) - c1]};
{
  my(K = 7);
  for(lin = 0, 1, for(m = 0, K - lin,
    my(E = eqsK(K, m, lin));
    \\ moment reduction: coefficients by type for member 1 (others zeroed), pair coefficient, constant (as in k7_verify.gp)
    my(Z = vector(4, j, my(e = E[j]); for(i = 2, m, e = subst(subst(e, eval(Str("u", i)), 0), eval(Str("v", i)), 0)); e));
    my(ty(e) = if(m == 0, vector(6), [polcoef(polcoef(e,1,v1),0,u1), polcoef(polcoef(e,1,v1),1,u1), polcoef(polcoef(e,1,v1),2,u1), polcoef(polcoef(e,1,v1),3,u1), polcoef(polcoef(e,2,v1),0,u1), polcoef(polcoef(e,2,v1),1,u1)]));
    my(pr = vector(4, j, if(m < 2, 0, my(e = E[j]); for(i = 3, m, e = subst(subst(e, eval(Str("u", i)), 0), eval(Str("v", i)), 0)); polcoef(polcoef(polcoef(polcoef(e,1,v1),1,v2),0,u1),0,u2))));
    my(cs = vector(4, j, my(e = E[j]); for(i = 1, m, e = subst(subst(e, eval(Str("u", i)), 0), eval(Str("v", i)), 0)); e));
    \\ verify the moment form exactly
    my(ok = 1);
    for(j = 1, 4, my(r = ty(Z[j]), ex = cs[j]);
      for(i = 1, m, my(u = eval(Str("u", i)), v = eval(Str("v", i))); ex += r[1]*v + r[2]*v*u + r[3]*v*u^2 + r[4]*v*u^3 + r[5]*v^2 + r[6]*v^2*u);
      for(i = 1, m, for(l = i + 1, m, ex += pr[j]*eval(Str("v", i))*eval(Str("v", l))));
      if(E[j] - ex != 0, ok = 0));
    \\ symbolic elimination as in k7_verify.gp
    my(red(j) = my(r = ty(Z[j])); r[1]*S + r[2]*P1 + r[3]*P2 + r[4]*P3 + r[5]*M + r[6]*N + pr[j]*(S^2 - M)/2 + cs[j]);
    my(e = vector(4, j, red(j)));
    print("m=", m, " lin=", lin, " moment-form exact: ", ok, "  e = ", e);
  ));
}
