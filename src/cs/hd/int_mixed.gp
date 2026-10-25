\\ Exact integral (mixed-sign) HD counterexample from the rational point M = -119/54 of the k = 5 two-dead family:
\\ slopes {1, 238/41, -476/95, 952/481, -119/54} scaled to integers; f_j = 1 + A_j x in Z[x]. Direct polynomial checks.
{
  my(s = [1, 238/41, -476/95, 952/481, -119/54], L = lcm(apply(denominator, s)), A = s * L / content(s * L));
  my(f = vector(5, j, 1 + A[j]*'x), F = prod(j = 1, 5, f[j]));
  print("integer slopes A_j = ", A);
  print("Korselt  (F-1) mod (f_j-1): ", vector(5, j, (F - 1) % (f[j] - 1)));
  my(gd = vector(5, j, poldegree(gcd(F + 1, f[j] + 1))));
  print("deg gcd(F+1, f_j+1):        ", gd);
  print("total Lucas defect = ", 5 - vecsum(gd), "  vs  deg F / 2 = ", poldegree(F)/2);
}
\\ Exact positive instance, k = 5: slopes 1, the three roots of G_5(a; 256), 256 (Sturm: 3 roots in (1, 256)).
{
  my(M = 256, G = (7*M^3 + 57*M^2 + 57*M + 7)*'a^3 + (-98*M^3 - 308*M^2 - 98*M)*'a^2 + (392*M^3 + 392*M^2)*'a - 448*M^3);
  print("G_5(a; 256) = ", G, "   irreducible over Q: ", polisirreducible(G), "   roots in (1,256): ", polsturm(G, [1, 256]));
  \\ exact check in the splitting field is the symbolic identity of c5exact.gp / gk_exact.gp; numeric echo:
  my(r = polrootsreal(G), v = [1, r[1], r[2], r[3], M]);
  print("  slopes: ", v);
  print("  prod_(i!=j)(1 - 2a_i/a_j): ", vector(5, j, prod(i = 1, 5, if (i == j, 1, 1 - 2*v[i]/v[j]))));
}
