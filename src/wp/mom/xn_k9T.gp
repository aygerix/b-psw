\\ xn line: k = 9 even parts T_n^e = (T_n + T~_n)/2 (n = 3..8) of the Vpoly coefficients in moment form.
\\ Vpoly(t) = sum_n T_n t^(8-n); T~ = reflected family (w -> -w).  Moments Mp_j = sum v^p w^j (w = u - 1/2).
read("xn_k9mom.gp");
{
  my(k = 9, ell = ellv(w + 1/2, v, k), ellr = ellv(1/2 - w, v, k), T = Tpoly(k));
  foreach([0, 1], lin,
    print("==== lin = ", lin);
    my(Te = vector(k - 1, n, 0));
    foreach([0, 1], refl,
      my(Lv = vector(k - 1, n, tomom(if(refl, ellr[n], ell[n])) + lin * ellin(k)[n]));
      for(n = 1, k - 1, my(e = polcoef(T, n, s)); for(j = 1, k - 1, e = subst(e, eval(Str("L", j)), Lv[j])); Te[n] += e/2));
    for(n = 1, k - 1, print("T", n, "e = ", Te[n])));
}
