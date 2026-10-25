\\ independent check of the B&B result at (-3/10, 4/5): evaluate sum alpha_i D_i DIRECTLY from Dof (polynomial division,
\\ no 3-body data) at the corner w = (LO,...,LO) and at random rational points of [LO, HI]^7, exactly.
read("pc7_cert_shape.gp");
read("pc7_cert_alpha1.gp");
SHAPE = [-3/10, 4/5]; LO = 6/25; HI = 8/5;
{
  my(R1 = SHAPE[1], R2 = SHAPE[2], Dl = (1 - R1^2)*(1 - R2^2), sub = [R1 + R2, R1*R2, 2*(R1 + R2)/Dl, 2*(1 - (R1 + R2)^2 + R1*R2)/Dl]);
  my(Ev(ws) = my(D = substvec(Dof(ws), [e1, e2, lam, mu], sub)); sum(i = 1, 8, ALPHA[i] * D[i]));
  print("E_alpha(LO,...,LO) = ", Ev(vector(7, j, LO)) * 1.);
  my(mn = oo);
  for(trial = 1, 300, my(ws = vector(7, j, LO + (HI - LO) * random(10^6) / 10^6), e = Ev(ws)); if(e < mn, mn = e));
  print("min over 300 random points of [LO,HI]^7: ", mn * 1.);
  \\ and outside the reality box: the all-tiny configuration (c -> infinity) has E -> kap_alpha < 0
  print("E_alpha(1/1000,...) = ", Ev(vector(7, j, 1/1000)) * 1., "   (negative: alpha is NOT a certificate on (0, oo)^7)");
}
