\\ pc7_cert_dsym.gp -- SYMBOLIC D-basis data (polynomials in e1, e2, lam, mu): kap (8), A (8 x 9, coeff of w^1..w^9),
\\ C (8 x 5, Hankel c_2..c_6) with D_i = kap_i + sum_j a_i(w_j) + sum_{i<j} b_i(w_i, w_j) + lam3_i e3, lam3 = [1,0,0,0,-1,0,0,0].
read("pc7_cert_shape.gp");
{
  my(f = "pc7_cert_dsym.txt", E0 = SYM0, A = vector(8, i, SYM1[i] - SYM0[i]));
  my(C = matrix(8, 5, i, n, polcoef(polcoef(SYM2[i] - SYM0[i] - A[i] - subst(A[i], w, v), 1, w), n, v)));
  for(i = 1, 8, my(Bp = SYM2[i] - SYM0[i] - A[i] - subst(A[i], w, v));
    if(Bp != sum(n = 2, 6, C[i, n - 1] * sum(a = 1, n - 1, w^a * v^(n - a))), error("Hankel structure fails, i = ", i)));
  write(f, "K ", E0);
  for(i = 1, 8, write(f, "A ", vector(9, m, polcoef(A[i], m, w))));
  for(i = 1, 8, write(f, "C ", C[i, ]));
  print("Hankel structure verified symbolically; w-degrees ", vector(8, i, poldegree(A[i], w)));
}
