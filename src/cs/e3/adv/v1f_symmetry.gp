\\ src/cs/e3/adv/v1f_symmetry.gp -- the rho-family has the involution rho -> 1/2 - rho (x -> 1 - x swaps members 2 and 3),
\\ so the scan's (rho, d) grid revisits the same n (also d = 2 with X even repeats d = 1).
fam(r) = my(B = 2*x^2 - (4*r+1)*x + 2*r - 2); [1 + B, 1 + B*x/(2*r), 1 + B*(x-1)/(2*r-1)];
{
  my(f = fam(r), g = fam(1/2 - r), fs = subst(f, x, 1 - x));
  print("f(rho)(1 - x) == [g1, g3, g2] with g = fam(1/2 - rho): ", fs == [g[1], g[3], g[2]]);
  print("rho = 1, x = 6: ", subst(fam(1), x, 6), "   rho = -1/2, x = -5: ", subst(fam(-1/2), x, -5));
}
