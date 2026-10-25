\\ src/cs/e3/adv/v1_basic.gp -- adversarial re-check of the E3 symbolic claims (independent code, general mu).
\\ Shape: f1 = 1 + c A, f2 = 1 + mu A (x - r2), f3 = 1 + la A (x - r3), A = x^2 + p x + q.
default(parisize, 200000000);
print("== V1a: Korselt equations with general mu (r2 = 0)");
{
  my(A = x^2 + p*x + q, f1 = 1 + c*A, f2 = 1 + mu*A*x, f3 = 1 + la*A*(x - r3), F = f1*f2*f3,
     H = (F - 1) / A);
  print("  (F-1) divisible by A: ", (F - 1) % A == 0);
  print("  H(0)  = ", factor(subst(H, x, 0)));
  print("  H(r3) = ", factor(subst(H, x, r3)));
  print("  claimed K2 - H(0) relation: H(0) == c - la r3 (1 + c q): ", subst(H, x, 0) == c - la*r3*(1 + c*q));
  print("  claimed K3: H(r3) == c (1 + mu r3 A(r3)) + mu r3: ", subst(H, x, r3) == c*(1 + mu*r3*subst(A, x, r3)) + mu*r3);
  \\ r2 = r3 case: H(r) = c
  my(f2b = 1 + mu*A*(x - r), f3b = 1 + la*A*(x - r), Hb = ((f1*f2b*f3b) - 1) / A);
  print("  r2 = r3 = r: H(r) = ", subst(Hb, x, r));
  \\ invariants under x -> beta x
  my(Ab = subst(A, x, beta*x) / beta^2, f2s = subst(f2, x, beta*x));
  print("  scaling: f2(beta x) - 1 == mu beta^3 * Abar * x : ", f2s - 1 == mu*beta^3*Ab*x);
}
print("== V1b: rho-family from invariants (I1, I2, I3) = (2, rho - 1, 1/rho), r3 = 1; compare with stated B-form");
famI(I1, I2, I3) = {
  my(c = I1, q = I2, mu = I3, la = c/(1 + c*q), p = (-1/c - 1/mu - 1 - q), A = x^2 + p*x + q);
  [1 + c*A, 1 + mu*A*x, 1 + la*A*(x - 1)] };
kor(f) = my(F = f[1]*f[2]*f[3]); vector(#f, j, (F - 1) % (f[j] - 1) == 0);
gdeg(f) = my(F = f[1]*f[2]*f[3]); vector(#f, j, poldegree(gcd(F + 1, f[j] + 1)));
{
  my(fr = famI(2, r - 1, 1/r), B = 2*x^2 - (4*r + 1)*x + 2*r - 2, fs = [1 + B, 1 + B*x/(2*r), 1 + B*(x - 1)/(2*r - 1)]);
  print("  famI(2, r-1, 1/r) == stated family: ", fr == fs);
  print("  Korselt over Q(r)[x]: ", kor(fs), "   gcd degrees: ", gdeg(fs));
  print("  f_j(2r): ", vector(3, j, subst(fs[j], x, 2*r)));
  my(F = fs[1]*fs[2]*fs[3]);
  print("  (F+1) % (f1+1) = ", (F + 1) % (fs[1] + 1));
  print("  f2+1 = ", factor(fs[2] + 1), "   f3+1 = ", factor(fs[3] + 1));
  \\ extra Lucas points: common roots of residual and F+1, via resultant in x
  my(e2 = polresultant(2*x^2 - x - 2, numerator(F + 1), x), e3 = polresultant(2*x^2 - 3*x - 1, numerator(F + 1), x));
  print("  Res(2x^2-x-2, F+1) = ", factor(numerator(e2)), "   Res(2x^2-3x-1, F+1) = ", factor(numerator(e3)));
  \\ double root at x = 2r of f_j + 1 also raises the gcd degree
  print("  f2+1 double root at 2r <=> ", factor(numerator(subst(deriv(fs[2] + 1, x), x, 2*r))),
        "   f3+1 double root at 2r <=> ", factor(numerator(subst(deriv(fs[3] + 1, x), x, 2*r))));
  print("  disc(B) = ", factor(poldisc(B, x)));
  \\ mu = 1 normalisation as stated
  my(w = 'w, pw = -(4*w^3 + 1)/(2*w), qw = (w^3 - 1)/w^2, cw = 2*w^2, lw = 2*w^3/(2*w^3 - 1), r3w = 1/w,
     Aw = x^2 + pw*x + qw, fw = [1 + cw*Aw, 1 + Aw*x, 1 + lw*Aw*(x - r3w)]);
  print("  mu=1 normalisation: Korselt ", kor(fw), "  gcd degrees ", gdeg(fw),
        "  invariants (c r3^2, q/r3^2, mu r3^3) = ", [cw*r3w^2, qw/r3w^2, r3w^3]);
  \\ example rho = 1, x = 2T
  my(f1e = subst(subst(fs, r, 1), x, 2*T));
  print("  rho=1, x=2T: ", f1e, "  at T=1: ", subst(f1e, T, 1));
}
print("== V1c: Lemma 5 curve (I1 = 1, I2 = 1/(2u) - 3/2, I3 = u) and P23");
{
  my(fu = famI(1, 1/(2*u) - 3/2, u));
  print("  generic u: Korselt ", kor(fu), " gcd degrees ", gdeg(fu));
  my(L = famI(1, 1/(2*(-1)) - 3/2, -1), Ly = subst(L, x, -Y));
  print("  u = -1, x -> -Y: ", Ly);

  print("  minpolys of 2cos(2pi/10), 2cos(2pi/14), 2cos(2pi/18): ",
        [algdep(2*cos(2*Pi/10), 2), algdep(2*cos(2*Pi/14), 3), algdep(2*cos(2*Pi/18), 3)]);
  \\ special rational u where the gcd pattern differs from the generic [0,2,0]
  my(bad = List());
  for (nu = -60, 60, for (de = 1, 30, if (gcd(nu, de) != 1 || nu == 0, next); my(uu = nu/de, ok = 1);
     my(q0 = 1/(2*uu) - 3/2); if (1 + q0 == 0, next);
     my(ff = famI(1, q0, uu), A0 = (ff[1] - 1));
     if (poldisc(A0) == 0, next);
     my(g = gdeg(ff)); if (g != [0, 2, 0], listput(bad, [uu, g]))));
  print("  rational u = nu/de (|nu| <= 60, de <= 30) with pattern != [0,2,0]: ", bad);
  my(P23 = famI(1, -4/3, 3)); print("  P23 members ", P23, "  Korselt ", kor(P23), "  gcd degrees ", gdeg(P23));
  \\ f1 = x^2 - x - 1/3 never integral on Q: brute check X/d, d <= 200, |X| <= 2000
  my(cnt = 0); for (d = 1, 200, for (X = -2000, 2000, if (gcd(X, d) == 1 && denominator(subst(P23[1], x, X/d)) == 1, cnt++)));
  print("  integral values of P23 f1 at X/d (d <= 200, |X| <= 2000): ", cnt);
}
print("== V1d: isolated-point polynomials");
{
  my(s6 = 47*y^6 + 102*y^5 - 36*y^4 + 189*y^3 - 162*y^2 - 2187*y - 729, s2 = y^2 - 42*y + 9);
  print("  sextic irreducible over Q: ", polisirreducible(s6), "  rational roots: ", nfroots(, s6));
  print("  quadratic disc = ", poldisc(s2), " = ", factor(poldisc(s2)), "  issquare: ", issquare(poldisc(s2)));
  print("  degree-18 in r3 = sextic(r3^3): ", subst(s6, y, r3^3) == 47*r3^18+102*r3^15-36*r3^12+189*r3^9-162*r3^6-2187*r3^3-729);
}
print("== V1e: dominant scan row n = 1152271");
{
  my(v = [43, 127, 211], n = 1152271);
  print("  product: ", vecprod(v) == n, "  primes: ", vector(3, j, isprime(v[j])), "  Korselt p-1 | n-1: ", vector(3, j, (n - 1) % (v[j] - 1) == 0));
  print("  rho = 1, x = 6 values: ", subst(subst([1 + (2*x^2 - 5*x), 1 + (2*x^2 - 5*x)*x/2, 1 + (2*x^2 - 5*x)*(x - 1)], r, 1), x, 6));
}
