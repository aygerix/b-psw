read("q22_moments.gp");
{
  my(k = 7, B = bnd(k, e1, e2, lam, mu));
  \\ symmetric slice e1 = 0 (=> lam = 0, mu = 2/(1+e2))
  my(Bs = vector(8, i, subst(subst(B[i], e1, 0), lam, 0)));
  for(i = 1, 8, print("Bsym", i, " = ", Bs[i]));
}
