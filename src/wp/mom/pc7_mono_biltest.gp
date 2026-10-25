\\ test: random rational pencil (k = 7, not Williams); bilinear residual must equal ((V - R) * D) mod H^(M+1) etc., R the exact rational fn.
read("pc7_mono_bil.gp");
{
  my(k = 7, M = 4, r1 = 1/3, r2 = -1/5, cs = [2, 3, 7/2, 5, 11/3, 13/7, 9/4]);
  my(E1 = r1 + r2, E2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2), lam0 = 2*E1/Dl, mu0 = 2*(1 - E1^2 + E2)/Dl);
  my(Ev = vector(M, m, polcoef(prod(j = 1, k, 1 + cs[j]*x), m, x)));
  my(vw = VW(k, E1, E2, lam0, mu0));
  vw = vector(2, i, my(b = vw[i]); for(a = 1, 2*k, b = subst(b, eval(Str("p", a)), sum(j = 1, k, (1/cs[j])^a))); b);
  my(fo = bil(k, E1, E2, lam0, mu0, M, Ev, vw));
  my(G = t^2 - 1, H = (t - r1)*(t - r2), f0 = 1 + G*(lam0*t + mu0), f = vector(k, j, f0 + cs[j]*G*H), F = prod(j = 1, k, f[j]));
  my(Km = prod(j = 1, k, (f[j] - 1)/G), Kp = prod(j = 1, k, (f[j] + 1)/H), HM = H^(M+1), GM = G^(M+1));
  my(Vd = (F - 1) \ (G*Km), Wd = (F + 1) \ (H*Kp));
  my(RVt = lift(Mod(F - 1, HM) / Mod(G*Km, HM)), RWt = lift(Mod(F + 1, GM) / Mod(H*Kp, GM)));
  my(DV = prod(j = 1, k, 1 - (f[j] + 1)/2), DW = prod(j = 1, k, 1 + (f[j] - 1)/2));
  my(dV = ((Vd - RVt) * DV) % HM, dW = ((Wd - RWt) * DW) % GM);
  my(dir = concat(vector(2*M + 2, i, polcoef(dV, i - 1, t)), vector(2*M + 2, i, polcoef(dW, i - 1, t))));
  print("bilinear formal - direct: ", fo - dir, "   (direct nonzero: ", #select(x -> x != 0, dir), ")");
}
