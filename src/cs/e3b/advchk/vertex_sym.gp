\\ src/cs/e3b/advchk/vertex_sym.gp -- exact check of the symmetric "vertex FSP" Sum m = 6 points (the 3 "local dim 2" hits
\\ of out/sym_21111.log).  A = x^2 - 1, reflection pairs (2,3), (4,5): f_3(x) = f_2(-x), f_5(x) = f_4(-x).
\\ Vertex FSP: f_j(0) = -1 for all j  <=>  c = 2, l_2 r_2 = -2, l_4 r_4 = -2.  Then K_2 = K_3 and K_4 = K_5 (mirror), so the
\\ family is cut out by two equations in (r, s) = (r_2, r_4).  Check: 0-dimensional?  rational solutions?  Sum m?
A = x^2 - 1;
fam(r, s) = my(l2 = -2/r, l4 = -2/s); [1 + 2*A, 1 + l2*A*(x - r), 1 - l2*A*(x + r), 1 + l4*A*(x - s), 1 - l4*A*(x + s)];
\\ K-form numerators
Kform(f, j, rj) = my(a = subst(A, x, rj), P = prod(i = 1, 5, if (i == j, 1, subst(f[i], x, rj)))); (P - 1) / a;
{
  my(a2 = r^2 - 1, a4 = s^2 - 1, n2, n4, res, fa);
  \\ s^2 (prod_{i != 2} f_i(r) - 1) / a2  and  r^2 (prod_{i != 4} f_i(s) - 1) / a4, explicit polynomials
  n2 = ((2*r^2 - 1)*(4*r^2 - 3)*(s - 2*a2*(r - s))*(s + 2*a2*(r + s)) - s^2) / a2;
  n4 = ((2*s^2 - 1)*(4*s^2 - 3)*(r - 2*a4*(s - r))*(r + 2*a4*(s + r)) - r^2) / a4;
  if (type(n2) != "t_POL" || type(n4) != "t_POL", error("not polynomial"));
  \\ cross-check against the family definition at a random rational point
  my(f = fam(3/7, -5/11), K2 = Kform(f, 2, 3/7), K4 = Kform(f, 4, -5/11));
  print("explicit forms match family: ", K2 * (-5/11)^2 == subst(subst(n2, r, 3/7), s, -5/11), " ", K4 * (3/7)^2 == subst(subst(n4, r, 3/7), s, -5/11));
  print("deg_r, deg_s of n2: ", poldegree(n2, r), " ", poldegree(n2, s), ";  n4: ", poldegree(n4, r), " ", poldegree(n4, s));
  res = polresultant(n2, n4, s);
  print("resultant in s is nonzero: ", res != 0, ", degree in r: ", poldegree(res, r));
  fa = factor(res);
  print("factor degrees (r): ", vector(#fa~, i, [poldegree(fa[i, 1], r), fa[i, 2]]));
  print("linear factors (rational r): ", select(p -> poldegree(p, r) == 1, fa[, 1]~));
  \\ numeric solutions with r, s nonzero and r^2 != 1, s^2 != 1, r^2 != s^2; recompute Sum m from the polynomials
  my(sols = List(), nreal = 0);
  for (i = 1, #fa~, my(p = fa[i, 1]); if (poldegree(p, r) < 1, next);
    foreach(polroots(p), rr,
      if (abs(rr) < 1e-8 || abs(rr^2 - 1) < 1e-8, next);
      my(q = gcd(subst(n2, r, rr), subst(n4, r, rr)), q2 = substpol(subst(n2, r, rr), s, s), ss);
      foreach(polroots(subst(n2, r, rr)), ss,
        if (abs(ss) < 1e-8 || abs(ss^2 - 1) < 1e-8 || abs(ss^2 - rr^2) < 1e-8, next);
        if (abs(subst(subst(n4, r, rr), s, ss)) > 1e-12 * (1 + abs(subst(subst(n4, r, rr), s, 0))), next);
        listput(sols, [rr, ss])));
  );
  print("candidate (r, s) solutions: ", #sols);
  my(best = 0, sums = List());
  foreach(sols, v,
    my(f = fam(v[1], v[2]), F = prod(i = 1, 5, f[i]), m = vector(5), real);
    for (j = 1, 5,
      my(N = f[j] + 1, rts = polroots(N), cnt = 0);
      foreach(rts, t, my(val = abs(subst(F + 1, x, t)), sc = 1 + vecprod(vector(5, i, abs(subst(f[i], x, t)))));
        if (val / sc < 1e-15, cnt++));
      m[j] = cnt);
    real = (abs(imag(v[1])) < 1e-15 && abs(imag(v[2])) < 1e-15);
    listput(sums, [vecsum(m), real]);
    if (vecsum(m) > best, best = vecsum(m)));
  print("max Sum m over vertex-sym points: ", best, ";  Sum m histogram: ", vecsort(apply(v -> v[1], Vec(sums))));
  print("real points: ", #select(v -> v[2], Vec(sums)));
}
