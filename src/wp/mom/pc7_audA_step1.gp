\\ pc7_audA_step1: Step 1 product formulas, relative errors per trial (debug of absolute-error spike)
default(realprecision, 120);
{
  setrand(12345);
  for(trial = 1, 20, my(k = if(trial % 2, 7, 9));
    my(R1 = (random(1999) - 999)/1000, R2 = (random(1999) - 999)/1000);
    if(R1 == R2 || abs(R1) == 1 || abs(R2) == 1, next);
    my(E1 = R1 + R2, E2 = R1*R2, D = (1 - R1^2)*(1 - R2^2), la = 2*E1/D, m = 2*(1 - E1^2 + E2)/D);
    my(cs = vector(k, j, (random(4001) - 2000)/317), G = x^2 - 1, H = (x - R1)*(x - R2), f0 = 1 + G*(la*x + m));
    my(F = prod(j = 1, k, f0 + cs[j]*G*H), worst = 0., big = 0.);
    for(j = 1, k, if(cs[j] == 0, next);
      foreach(polroots(cs[j]*H + la*x + m), z, my(u = subst(G*H, x, z), fz = subst(F, x, z), pz = prod(l = 1, k, if(l == j, 1, 1 + (cs[l] - cs[j])*u)));
        worst = max(worst, abs(fz - pz)/(1 + abs(pz))); big = max(big, abs(pz)));
      foreach(polroots(cs[j]*G + la*x + m + la*E1), z, my(u = subst(G*H, x, z), fz = subst(F, x, z), pz = (-1)^k*prod(l = 1, k, if(l == j, 1, 1 - (cs[l] - cs[j])*u)));
        worst = max(worst, abs(fz - pz)/(1 + abs(pz))); big = max(big, abs(pz))));
    \\ exact check: F - prod formula reduced mod K-_j and K+_j as polynomials (exact rationals)
    my(ex = 0);
    for(j = 1, k, if(cs[j] == 0, next); my(Km = cs[j]*H + la*x + m, Kp = cs[j]*G + la*x + m + la*E1, u = G*H);
      ex += ((F - prod(l = 1, k, if(l == j, 1, 1 + (cs[l] - cs[j])*u))) % Km != 0);
      ex += ((F - (-1)^k*prod(l = 1, k, if(l == j, 1, 1 - (cs[l] - cs[j])*u))) % Kp != 0));
    print("trial ", trial, " k=", k, " r=", [R1, R2], " rel err ", worst, "  max|prod| ", big, "  exact mod-K failures: ", ex));
}
