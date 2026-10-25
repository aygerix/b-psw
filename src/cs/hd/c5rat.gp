\\ rational points on the k = 5 two-dead component C5: P(a2, M) = 0 (a1 = 1, a5 = M); search M = u/v, |u|,|v| <= 300.
P(a, M) = (7*M^3 + 57*M^2 + 57*M + 7)*a^3 + (-98*M^3 - 308*M^2 - 98*M)*a^2 + (392*M^3 + 392*M^2)*a - 448*M^3;
{
  my(cnt = 0);
  for (u = -300, 300, for (v = 1, 300, if (gcd(u, v) != 1 || u == 0, next);
    my(M = u/v, Q = P('a, M));
    if (Q == 0, next);
    my(fa = factor(Q));
    for (i = 1, #fa~, if (poldegree(fa[i,1]) == 1, my(a2 = -polcoeff(fa[i,1], 0)/polcoeff(fa[i,1], 1));
       cnt++; if (cnt <= 40, print("M = ", M, "  a2 = ", a2))))));
  print("rational points found: ", cnt);
}
