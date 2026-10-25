\\ test of pc7_mono_hs.gp: random rational pencil (k = 7, not Williams): formal equations (p's, P's substituted) must equal
\\ the remainders (V_dir - R) mod H^(M+1), (W_dir - R') mod G^(M+1) of the exact rational functions R = (F-1)/(G Km), R' = (F+1)/(H Kp).
read("pc7_mono_hs.gp");
{
  my(k = 7, M = 4, r1 = 1/3, r2 = -1/5, cs = [2, 3, 7/2, 5, 11/3, 13/7, 9/4]);
  my(E1 = r1 + r2, E2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2), lam0 = 2*E1/Dl, mu0 = 2*(1 - E1^2 + E2)/Dl);
  my(Pv = vector(M, m, sum(j = 1, k, cs[j]^m)));
  my(vw = VW(k, E1, E2, lam0, mu0));
  vw = vector(2, i, my(b = vw[i]); for(a = 1, 2*k, b = subst(b, eval(Str("p", a)), sum(j = 1, k, (1/cs[j])^a))); b);
  my(fo = hsym(k, E1, E2, lam0, mu0, M, Pv, vw));
  my(G = t^2 - 1, H = (t - r1)*(t - r2), f0 = 1 + G*(lam0*t + mu0), f = vector(k, j, f0 + cs[j]*G*H), F = prod(j = 1, k, f[j]));
  my(Km = prod(j = 1, k, (f[j] - 1)/G), Kp = prod(j = 1, k, (f[j] + 1)/H));
  my(Vd = (F - 1) \ (G*Km), Wd = (F + 1) \ (H*Kp), HM = H^(M+1), GM = G^(M+1));
  my(RVt = lift(Mod(F - 1, HM) / Mod(G*Km, HM)), RWt = lift(Mod(F + 1, GM) / Mod(H*Kp, GM)));
  my(dV = (Vd - RVt) % HM, dW = (Wd - RWt) % GM);
  my(dir = concat(vector(2*M + 2, i, polcoef(dV, i - 1, t)), vector(2*M + 2, i, polcoef(dW, i - 1, t))));
  print("formal - direct: ", fo - dir);
  print("direct nonzero entries: ", #select(x -> x != 0, dir), " of ", #dir);
  \\ M = 1 symbolic vs the original B1..B8: same ideal? check that each hsym M=1 equation is a combination of B's (rank test at random p's)
}
