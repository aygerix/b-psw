\\ test phiVW against VW (power-sum moment formulas) at random numeric data: shape and w's -> V, W; phiVW(V) must equal W.
read("pc7_mono_phi.gp"); read("pc7_mono_hs.gp");
{
  my(k = 7, r1 = 1/3, r2 = -1/5, E1 = r1 + r2, E2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2), lam0 = 2*E1/Dl, mu0 = 2*(1 - E1^2 + E2)/Dl);
  my(ws = [1/2, 1/3, 2/7, 1/5, 3/11, 7/13, 4/9]);
  my(vw = VW(k, E1, E2, lam0, mu0));
  vw = vector(2, i, my(b = vw[i]); for(a = 1, 2*k, b = subst(b, eval(Str("p", a)), sum(j = 1, k, ws[j]^a))); b);
  my(Wp = phiVW(k, vw[1], E1, E2, lam0, mu0, 1/lam0));
  print("phi(V) - W = ", Wp - vw[2]);
  \\ also with free (formal) p's: random rational p's, not from members
  my(pr = vector(2*k, a, random(100)/7 - 5));
  my(vw2 = VW(k, E1, E2, lam0, mu0)); vw2 = vector(2, i, my(b = vw2[i]); for(a = 1, 2*k, b = subst(b, eval(Str("p", a)), pr[a])); b);
  print("formal p's: phi(V) - W = ", phiVW(k, vw2[1], E1, E2, lam0, mu0, 1/lam0) - vw2[2]);
  \\ k = 5 too
  my(vw5 = VW(5, E1, E2, lam0, mu0)); vw5 = vector(2, i, my(b = vw5[i]); for(a = 1, 10, b = subst(b, eval(Str("p", a)), pr[a])); b);
  print("k = 5 formal: ", phiVW(5, vw5[1], E1, E2, lam0, mu0, 1/lam0) - vw5[2]);
}
