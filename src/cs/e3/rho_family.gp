\\ src/cs/e3/rho_family.gp -- E3 confirm branch: the rational shared-point family ("rho-family") of k = 3 shared-A
\\ Korselt families with all three members Lucas-good.  Polynomial variable x (= t), parameter r (= rho).
\\   B = 2 x^2 - (4 r + 1) x + 2 r - 2  (= 2A),  f1 = 1 + B,  f2 = 1 + B x/(2 r),  f3 = 1 + B (x - 1)/(2 r - 1).
\\ Shared Lucas point x = 2 r (f1 = f2 = f3 = -1).
fam(r) = my(B = 2*x^2 - (4*r+1)*x + 2*r - 2); [1 + B, 1 + B*x/(2*r), 1 + B*(x-1)/(2*r-1)];
{
  my(f = fam('r), F = f[1]*f[2]*f[3]);
  print("Symbolic checks over Q(r)[x]:");
  for (j = 1, 3, print("  Korselt: (F-1) mod (f", j, "-1) = ", (F - 1) % (f[j] - 1)));
  for (j = 1, 3, print("  f", j, "(2r) = ", subst(f[j], x, 2*'r)));
  for (j = 1, 3, my(h = gcd(F + 1, f[j] + 1)); print("  h", j, " = gcd(F+1, f", j, "+1) = ", h / pollead(h), "   [deg ", poldegree(h), " of ", poldegree(f[j] + 1), "]"));
  print("  f1+1 = ", factor(f[1] + 1));
  print("  f2+1 = ", factor(f[2] + 1));
  print("  f3+1 = ", factor(f[3] + 1));
  my(r2 = 2*x^2 - x - 2, r3 = 2*x^2 - 3*x - 1);
  my(e2 = numerator(polresultant(r2, numerator(1 - f[1]*f[3]), x)), e3 = numerator(polresultant(r3, numerator(1 - f[1]*f[2]), x)));
  print("  extra Lucas point for f2 (root of 2x^2-x-2 with F = -1) <=> ", factor(e2));
  print("  extra Lucas point for f3 (root of 2x^2-3x-1 with F = -1) <=> ", factor(e3));
  print("  disc(B) = ", factor(poldisc(2*x^2 - (4*'r+1)*x + 2*'r - 2, x)));
  print("  (F+1) / ((f1+1) * (x - 2r)) has content-free form: ", factor(numerator((F+1)/((f[1]+1)))));
}
