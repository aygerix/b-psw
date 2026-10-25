\\ divisibility of the top coefficients by powers of lam (polynomials in e1, e2, lam, mu, before lam, mu are specialised):
\\ for every D_i, A_m (m >= 5) and c_n (n >= 5): largest k with lam^k | coefficient.
read("pc7_cert_shape.gp");
{
  my(K = SYM0, A = vector(8, i, SYM1[i] - SYM0[i]));
  my(C = matrix(8, 5, i, n, polcoef(polcoef(SYM2[i] - SYM0[i] - A[i] - subst(A[i], w, v), 1, w), n, v)));
  my(vl(P) = if(P == 0, oo, valuation(P, lam)));
  for(i = 1, 8, print("D", i, ": lam-valuations A1..A9 = ", vector(9, m, vl(polcoef(A[i], m, w))), "   c2..c6 = ", vector(5, n, vl(C[i, n]))));
}
