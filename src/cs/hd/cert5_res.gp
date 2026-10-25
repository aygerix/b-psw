\\ Independent check of the k = 5 positive certificate (cert5.py, a1 = 1, a5 = 135): eliminate a4, a3 by resultants and
\\ locate the real root of the eliminant near a2 = 1.99924256189045; then back-solve a3, a4 and evaluate all five members.
default(realprecision, 80);
E(v, j) = prod(i = 1, #v, if (i == j, 1, v[j] - 2*v[i])) - v[j]^(#v - 1);
backsolve(a2, r23, e2) =
{
  my(cand = polroots(subst(r23, 'a2, a2)));
  for (u = 1, #cand,
    my(a3 = cand[u]);
    if (abs(imag(a3)) > 1e-30 || abs(real(a3) - 4.0436871780509) > 1e-6, next);
    a3 = real(a3);
    my(c4 = polroots(substvec(e2, ['a2, 'a3], [a2, a3])));
    for (w = 1, #c4,
      my(a4 = c4[w]);
      if (abs(imag(a4)) > 1e-30 || abs(real(a4) - 7.4630764501653) > 1e-6, next);
      a4 = real(a4);
      my(vv = [1, a2, a3, a4, 135]);
      print("  a = ", vv);
      print("  prod_(i!=j)(1 - 2a_i/a_j), j = 1..5: ", vector(5, j, prod(i = 1, 5, if (i == j, 1, 1 - 2*vv[i]/vv[j]))))));
}
{
  my(v = [1, 'a2, 'a3, 'a4, 135]);
  my(e2 = E(v, 2), e3 = E(v, 3), e4 = E(v, 4));
  my(r23 = polresultant(e2, e3, 'a4), r24 = polresultant(e2, e4, 'a4));
  my(R = polresultant(r23, r24, 'a3));
  print("eliminant degree in a2: ", poldegree(R, 'a2));
  my(fa = factor(R), hit = 0);
  for (i = 1, #fa~,
    my(P = fa[i, 1]);
    if (poldegree(P, 'a2) < 1, next);
    my(rt = polrootsreal(P));
    for (t = 1, #rt,
      if (abs(rt[t] - 1.99924256189045) < 1e-9,
        hit = 1;
        print("irreducible factor of degree ", poldegree(P), " has the real root a2 = ", rt[t]);
        backsolve(rt[t], r23, e2))));
  if (!hit, print("NO real root of the eliminant near the certified a2"));
}
