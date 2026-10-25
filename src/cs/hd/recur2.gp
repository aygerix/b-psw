\\ Symbolic two-dead families for k = 5, 7, 9, 11: dead slopes 1 and M (d = 1, 1/M), constraint in q1 factored over Q(M);
\\ for each factor linear in q1, the good-member polynomial G_k(a; M) := reciprocal of P (roots = good slopes a = -1/root).
recP(k, b0, b1, q1) = {
  my(w = vector(k, i, 2^(k-2-(i-1))), q0 = 2^(3-k)*b0, p = vector(k+1));
  p[k-1] = 1; p[k] = 0; p[k+1] = 0;
  forstep (n = k-1, 2, -1,
    my(t = p[n] * (w[n] * q1 - b1) + if (n <= k-2, p[n+1] * (w[n+1] * q0 - b0), 0));
    p[n-1] = -t / (w[n-1] - 1));
  my(con = p[1] * (w[1] * q1 - b1) + p[2] * (w[2] * q0 - b0));
  [Pol(vector(k-1, i, p[k-i]), 'x), con] };
E(v, j) = prod(i = 1, #v, if (i == j, 1, v[j] - 2*v[i])) - v[j]^(#v - 1);
{
  for (k = 5, 11, if (k % 2 == 0, next);
    my(r = recP(k, 1/'M, 1 + 1/'M, 'q), con = numerator(r[2]), fa = factor(con));
    print("k = ", k, ": constraint factors (deg in q): ", vector(#fa~, i, poldegree(fa[i,1], 'q)));
    for (i = 1, #fa~, if (poldegree(fa[i,1], 'q) != 1, next);
      my(q1 = -polcoeff(fa[i,1], 0, 'q)/polcoeff(fa[i,1], 1, 'q));
      my(P = subst(r[1], 'q, q1));
      \\ good slopes a = -1/root: G(a) = a^(k-2) P(-1/a)
      my(G = numerator(subst(P, 'x, -1/'a) * 'a^(k-2)));
      G = G / content(G);
      \\ numeric test at M = 2^(k+6): are the roots real, positive, distinct, and all good?
      my(MM = 2^(k+6), Gn = subst(G, 'M, MM));
      if (poldegree(Gn, 'a) != k-2, print("   q1 = ", q1, ": degenerate"); next);
      my(rt = polroots(Gn), isreal = vecmax(apply(z -> abs(imag(z)), rt)) < 1e-20);
      if (!isreal, print("   q1 = ", q1, ": roots not all real at M = ", MM); next);
      my(v = concat(concat([1], vector(#rt, j, real(rt[j]))), [MM]));
      my(res = vector(k, j, E(v, j) / v[j]^(k-1)));
      print("   q1 = ", q1, ":  at M = ", MM, " good slopes ", vector(#rt, j, real(rt[j])), "  residuals E_j/a_j^(k-1) = ", apply(z -> if (abs(z) < 1e-25, 0, z), res));
      if (k <= 7, print("   G_", k, "(a; M) = ", G))));
}
