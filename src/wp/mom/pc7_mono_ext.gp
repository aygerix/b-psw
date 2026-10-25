\\ explicit extreme-member sign functions for real (2,2) pencils.
\\ K-(c) = c H + L, K+(c) = c G + L' (multiplied by Dl > 0 to keep polynomial coefficients; same roots).
\\ For a real quadratic Q = a t^2 + b t + q0 with non-real roots z, zbar and a real polynomial A:
\\ Re A(z) = alpha * (-b/(2a)) + beta, where alpha t + beta = A mod Q.
default(parisizemax, 1000000000);
t; c; r1; r2;
reA(A, Q) = {my(R = A % Q, a = polcoef(Q, 2, t), b = polcoef(Q, 1, t)); polcoef(R, 1, t) * (-b/(2*a)) + polcoef(R, 0, t)};
{
  my(e1 = r1 + r2, e2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2), lamN = 2*e1, muN = 2*(1 - e1^2 + e2));
  my(G = t^2 - 1, H = t^2 - e1*t + e2, LN = lamN*t + muN, LpN = lamN*t + muN + lamN*e1);
  my(KmN = c*Dl*H + LN, KpN = c*Dl*G + LpN);
  my(Rm = reA(G*H, KmN), Rp = reA(G*H, KpN));
  print("rho-(c) = ", factor(numerator(Rm)), "  /  ", factor(denominator(Rm)));
  print("rho+(c) = ", factor(numerator(Rp)), "  /  ", factor(denominator(Rp)));
  print("disc K- (times Dl^2) = ", factor(poldisc(KmN, t))); print("disc K+ (times Dl^2) = ", factor(poldisc(KpN, t)));
  my(sub(X) = subst(subst(subst(X, r1, -3/10), r2, 6/10), c, 13/10));
  print("check rho- = ", sub(Rm) * 1., "  rho+ = ", sub(Rp) * 1., "   (python: 3.894655997, 0.817732920)");
}
