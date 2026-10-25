\\ test: at a random rational shape and random P's, a's, the reduced system equals the W-part of hsym evaluated at the recovered p's,
\\ and the V-part of hsym vanishes there.
read("pc7_mono_red.gp");
{
  my(k = 7, M = 3, r1 = 1/3, r2 = -1/5, E1 = r1 + r2, E2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2), lam0 = 2*E1/Dl, mu0 = 2*(1 - E1^2 + E2)/Dl);
  my(Pn = [3/2, 5/7, -2/3], an = [7/11]);
  my(R = redsys(k, M, E1, E2, lam0, mu0, Pn));
  print("free a's: ", R[3]);
  my(sol = R[2]); for(i = 1, #R[3], sol = subst(sol, R[3][i], an[i]));
  my(eqr = R[1]); for(i = 1, #R[3], eqr = subst(eqr, R[3][i], an[i]));
  my(vw = VW(k, E1, E2, lam0, mu0)); vw = vector(2, i, my(b = vw[i]); for(a = 1, 2*k - 1, b = subst(b, eval(Str("p", a)), sol[a])); b);
  my(hs = hsym(k, E1, E2, lam0, mu0, M, Pn, vw));
  print("V-part of hsym at recovered p: ", hs[1..2*M+2]);
  print("W-part minus reduced: ", hs[2*M+3..4*M+4] - eqr);
}
