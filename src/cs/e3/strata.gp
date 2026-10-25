\\ src/cs/e3/strata.gp -- E3: rational (twisted) families on each Lucas gcd-degree stratum found by e3_full.sing.
\\ Invariants under x -> beta x + gamma:  I1 = c r3^2,  I2 = q / r3^2,  I3 = mu r3^3  (g2 = mu (x - r2), r2 = 0).
\\ A rational family with invariants (I1, I2, I3): r3 = 1, c = I1, q = I2, mu = I3, lambda = c/(1 + c q),
\\ p = -(c + mu)/(c mu) - 1 - q  (Korselt K2, K3 with mu).  f1 = 1 + c A, f2 = 1 + mu A x, f3 = 1 + lambda A (x - 1).
famI(I1, I2, I3) = {
  my(c = I1, q = I2, mu = I3, la = c/(1 + c*q), p = -(c + mu)/(c*mu) - 1 - q, A = x^2 + p*x + q);
  [1 + c*A, 1 + mu*A*x, 1 + la*A*(x - 1)] };
info(f, tag) = {
  my(F = f[1]*f[2]*f[3], kor = vector(3, j, (F - 1) % (f[j] - 1) == 0),
     gd = vector(3, j, poldegree(gcd(F + 1, f[j] + 1))), df = vector(3, j, poldegree(f[j] + 1) - gd[j]));
  print(tag, ":  Korselt ", kor, "  gcd degrees ", gd, "  Lucas defect degrees ", df, "  total ", vecsum(df),
        "  disc(A)!=0: ", poldisc((f[1] - 1)/pollead(f[1] - 1)) != 0);
  df };
{
  print("-- rho-family (I1 = 2, I2 = rho - 1, I3 = 1/rho), symbolic rho:");
  info(famI(2, 'r - 1, 1/'r), "  rho generic");
  print("-- f1-full component 2 (I1 = s, I3 = s(s+2)/2, I2 = -(s+4)/(s(s+2))), symbolic s:");
  info(famI('s, -('s + 4)/('s*('s + 2)), 's*('s + 2)/2), "  s generic");
  for (s = 1, 4, info(famI(s, -(s + 4)/(s*(s + 2)), s*(s + 2)/2), Str("  s = ", s)));
  print("-- f2 gcd>=2 curve (I1 = 1, I2 = 1/(2 I3) - 3/2), symbolic u = I3; Lemma 5 is u = -1:");
  info(famI(1, 1/(2*'u) - 3/2, 'u), "  u generic");
  my(L5 = famI(1, 1/(2*(-1)) - 3/2, -1)); info(L5, "  u = -1");
  print("     Lemma 5 check: f_j(x) at x -> substitution? f = ", L5);
  print("-- f3 gcd>=2 curve (I1 = 1, I2 = 2/I3 - 2), symbolic u:");
  info(famI(1, 2/'u - 2, 'u), "  u generic");
  print("-- f2>=2 and f3>=2 point (I1 = 1, I3 = 3, I2 = -4/3):");
  my(P23 = famI(1, -4/3, 3)); info(P23, "  P23"); print("     P23 members: ", P23);
  print("     P23 factorisations of f_j + 1: ", vector(3, j, factor(P23[j] + 1)));
  my(F = P23[1]*P23[2]*P23[3]); print("     P23 gcds: ", vector(3, j, gcd(F + 1, P23[j] + 1)));
}
