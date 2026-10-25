\\ Jacobian ranks of the M = 1 (2,2)-pencil system in the free power sums p1..p9 at a fixed shape and random p's
read("pc7_mono_hs.gp");
{
  my(k = 7, r1 = 1/3, r2 = -1/5, E1 = r1 + r2, E2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2), lam0 = 2*E1/Dl, mu0 = 2*(1 - E1^2 + E2)/Dl);
  my(B = bnd(k, E1, E2, lam0, mu0), Hs = hsym(k, E1, E2, lam0, mu0, 1, [P1]));
  my(pv = vector(9, a, eval(Str("p", a))), pt = vector(9, a, random(1000)/37 + 1));
  my(jac(L) = matrix(#L, 9, i, a, my(d = deriv(L[i], pv[a])); for(b = 1, 9, d = subst(d, pv[b], pt[b])); d));
  print("rank B1..B8: ", matrank(jac(B)), "  V-side B1..B4: ", matrank(jac(B[1..4])), "  W-side B5..B8: ", matrank(jac(B[5..8])));
  print("rank hsym M=1: ", matrank(jac(Hs)), "  V: ", matrank(jac(Hs[1..4])), "  W: ", matrank(jac(Hs[5..8])));
  print("P1 present in hsym M=1: ", #select(x -> deriv(x, P1) != 0, Hs));
  print("joint rank B with hsym: ", matrank(jac(concat(B, Hs))));
  \\ which p's appear in each B
  for(i = 1, 8, print("B", i, " vars: ", select(a -> deriv(B[i], pv[a]) != 0, vector(9, a, a))));
}
