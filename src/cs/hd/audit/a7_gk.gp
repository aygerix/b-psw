\\ AUDIT: take G_k from the line's recurrence (recP copied verbatim from gk_exact.gp), then check NUMERICALLY with own code
\\ that slopes {1, roots of G_k(.;M), M} have members 2..k-1 good (R_j = 1) at M = 2^(k+3) and M = 10^5, and that the
\\ roots coincide with the Newton/continuation doubling chain of a2_chain.gp.
default(realprecision, 200);
recP(k, b0, b1, q1) = {
  my(w = vector(k, i, 2^(k-2-(i-1))), q0 = 2^(3-k)*b0, p = vector(k+1));
  p[k-1] = 1; p[k] = 0; p[k+1] = 0;
  forstep (n = k-1, 2, -1,
    my(t = p[n] * (w[n] * q1 - b1) + if (n <= k-2, p[n+1] * (w[n+1] * q0 - b0), 0));
    p[n-1] = -t / (w[n-1] - 1));
  my(con = p[1] * (w[1] * q1 - b1) + p[2] * (w[2] * q0 - b0));
  [Pol(vector(k-1, i, p[k-i]), 'x), con] };
Rj(s, j) = prod(i=1, #s, if(i==j, 1, 1 - 2*s[i]/s[j]));
{
  for (k = 5, 13, if (k % 2 == 0, next);
    my(b1 = 1 + 1/'M, b0 = 1/'M, r = recP(k, b0, b1, b1/2^(k-2)));
    my(G = numerator(subst(r[1], 'x, -1/'y) * 'y^(k-2))); G = G / content(G);
    foreach([2^(k+3), 10^5, 3/7, -5], MM,
      my(Gn = subst(G, 'M, MM)); if (poldegree(Gn) < k-2, printf("k=%2d M=%s: degree drops to %d (leading coefficient vanishes)\n", k, MM, poldegree(Gn)); next); my(rr = polroots(Gn), re = [real(z) | z <- rr, abs(imag(z)) < 1e-60]);
      my(s = concat([1], concat(vecsort(re), [MM])));
      if (#re < k-2, printf("k=%2d M=%s: only %d real roots (complex family)\n", k, MM, #re);
         my(sc = concat([1], concat(Vec(rr~), [MM])), Rc = vector(k, j, Rj(sc, j)));
         printf("      complex check: max |R_j - 1| over middle = %.3g ; |R_1-1| = %.3g ; |R_k-1| = %.3g\n",
                vecmax(vector(k-2, j, abs(Rc[j+1]-1))), abs(Rc[1]-1), abs(Rc[k]-1)); next);
      my(R = vector(k, j, Rj(s, j)));
      printf("k=%2d M=%-8s real roots=%d  max|R_j-1| middle=%.3g  R_1-1=%.3g R_k-1=%.3g  positive=%d\n",
             k, MM, #re, vecmax(vector(k-2, j, abs(R[j+1]-1))), R[1]-1, R[k]-1, vecmin(re) > 0)));
}
