\\ For rational points (M, a2) on C5, recover a3, a4 exactly (common rational roots of E_2, E_3, E_4) and print full 5-tuples.
E(v, j) = prod(i = 1, #v, if (i == j, 1, v[j] - 2*v[i])) - v[j]^(#v - 1);
P(a, M) = (7*M^3 + 57*M^2 + 57*M + 7)*a^3 + (-98*M^3 - 308*M^2 - 98*M)*a^2 + (392*M^3 + 392*M^2)*a - 448*M^3;
H = 400;
{
  my(cnt = 0, pos = 0);
  for (u = -H, H, for (w = 1, H, if (gcd(u, w) != 1 || u == 0, next);
    my(M = u/w, Q = P('a, M)); if (Q == 0, next);
    my(fa = factor(Q));
    for (i = 1, #fa~, if (poldegree(fa[i,1]) != 1, next);
      my(a2 = -polcoeff(fa[i,1], 0)/polcoeff(fa[i,1], 1));
      my(v = [1, a2, 'a3, 'a4, M], e2 = E(v, 2), e3 = E(v, 3), e4 = E(v, 4));
      my(r = polresultant(e2, e3, 'a4), g = gcd(r, polresultant(e2, e4, 'a4)));
      if (poldegree(g, 'a3) < 1, next);
      my(fg = factor(g));
      for (s = 1, #fg~, if (poldegree(fg[s,1], 'a3) != 1, next);
        my(a3 = -polcoeff(fg[s,1], 0)/polcoeff(fg[s,1], 1));
        my(h = gcd(gcd(subst(e2, 'a3, a3), subst(e3, 'a3, a3)), subst(e4, 'a3, a3)));
        if (poldegree(h, 'a4) < 1, next);
        my(fh = factor(h));
        for (t = 1, #fh~, if (poldegree(fh[t,1], 'a4) != 1, next);
          my(a4 = -polcoeff(fh[t,1], 0)/polcoeff(fh[t,1], 1), vv = [1, a2, a3, a4, M]);
          if (#Set(vv) < 5 || vecmin(apply(abs, vv)) == 0, next);
          my(good = vector(5, j, E(vv, j) == 0));
          cnt++;
          my(allpos = vecmin(vv) > 0 || vecmax(vv) < 0);
          if (allpos, pos++);
          if (cnt <= 25 || allpos, print(if (allpos, "SAME-SIGN ", "mixed     "), "a = ", vv, "  good = ", good)))))));
  print("full rational 5-tuples: ", cnt, "   same-sign: ", pos);
}
