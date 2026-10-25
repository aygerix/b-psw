\\ Exact two-dead families G_k(a; M) for k = 5..13 (odd): derive G_k from the recurrence with q1 = b1/2^(k-2), then PROVE
\\ symbolically (in Q(M)[y]) that every root y of G_k is Lucas-good:  2^(k-2) (y-2)(y-2M) G(y/2) + c y^k == 0 mod G(y),
\\ and give an exact positive instance (Sturm count of roots in (1, M) at M = 2^(k+3)).
recP(k, b0, b1, q1) = {
  my(w = vector(k, i, 2^(k-2-(i-1))), q0 = 2^(3-k)*b0, p = vector(k+1));
  p[k-1] = 1; p[k] = 0; p[k+1] = 0;
  forstep (n = k-1, 2, -1,
    my(t = p[n] * (w[n] * q1 - b1) + if (n <= k-2, p[n+1] * (w[n+1] * q0 - b0), 0));
    p[n-1] = -t / (w[n-1] - 1));
  my(con = p[1] * (w[1] * q1 - b1) + p[2] * (w[2] * q0 - b0));
  [Pol(vector(k-1, i, p[k-i]), 'x), con] };
{
  for (k = 5, 13, if (k % 2 == 0, next);
    my(b1 = 1 + 1/'M, b0 = 1/'M, q1 = b1/2^(k-2), r = recP(k, b0, b1, q1));
    my(G = numerator(subst(r[1], 'x, -1/'y) * 'y^(k-2))); G = G / content(G);
    my(c = pollead(G, 'y), R = 2^(k-2)*('y - 2)*('y - 2*'M)*subst(G, 'y, 'y/2) + c*'y^k);
    my(MM = 2^(k+3), Gn = subst(G, 'M, MM));
    print("k = ", k, ": constraint at q1 = b1/2^(k-2): ", simplify(r[2]), ";  identity remainder: ", R % G,
          ";  deg_y G = ", poldegree(G, 'y), ";  M = ", MM, ": real roots in (1, M) = ", polsturm(Gn, [1, MM]),
          ", G(1) != 0: ", subst(Gn, 'y, 1) != 0, ", G(M) != 0: ", subst(Gn, 'y, MM) != 0, ", squarefree: ", poldisc(Gn) != 0));
}
