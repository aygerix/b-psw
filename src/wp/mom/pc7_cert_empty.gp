\\ exact check that the reality interval I(shape) is empty: the sets {q1 < 0} and {q2 < 0} (w > 0) are disjoint,
\\ q1 = (w lam - e1)^2 - 4 (e2 + w mu)  [K- = cH + L has no real root],  q2 = w^2 lam^2 - 4 w (mu + lam e1) + 4  [K+ = cG + L'].
chk(R1, R2) = {my(E1 = R1 + R2, E2 = R1*R2, Dl = (1 - R1^2)*(1 - R2^2), LA = 2*E1/Dl, MU = 2*(1 - E1^2 + E2)/Dl);
  my(q1 = (x*LA - E1)^2 - 4*(E2 + x*MU), q2 = x^2*LA^2 - 4*x*(MU + LA*E1) + 4);
  my(i1 = polrootsreal(q1), i2 = polrootsreal(q2));
  \\ both open upward (LA != 0); negativity sets are (roots) intervals
  my(lo = max(i1[1], i2[1]), hi = min(i1[#i1], i2[#i2]));
  print([R1, R2], ": {q1<0} = ", i1, "  {q2<0} = ", i2, "  intersection empty: ", lo >= hi,
        "  (exact: resultant sign check q1 at the q2-interval ends) ", [sign(subst(q1, x, i2[1])), sign(subst(q1, x, i2[2]))])};
chk(13/20, 17/20); chk(17/20, 17/20);
