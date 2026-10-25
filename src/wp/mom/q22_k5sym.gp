read("q22_moments.gp");
{
  my(k = 5, B = bnd(k, e1, e2, lam, mu));
  for(i = 1, #B, print("B", i, " vars ", variables(B[i]), "  deg in p's: ", vector(5, a, poldegree(B[i], eval(Str("p", a))))));
  for(i = 1, #B, print("B", i, " = ", B[i]));
}
