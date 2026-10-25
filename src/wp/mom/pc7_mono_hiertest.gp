\\ test of pc7_mono_hier.gp against a direct construction (random rational pencil, k = 7, not a Williams family):
\\ formal residual (moments + local P-series) must equal V_dir(r + h) - [(F - 1)/(G Km)](r + h), i.e. minus the proper part.
read("pc7_mono_hier.gp");
{
  my(k = 7, M = 4, r1 = 1/3, r2 = -1/5, cs = [2, 3, 7/2, 5, 11/3, 13/7, 9/4]);
  my(E1 = r1 + r2, E2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2), lam0 = 2*E1/Dl, mu0 = 2*(1 - E1^2 + E2)/Dl);
  my(Pv = vector(M, m, sum(j = 1, k, cs[j]^m)));
  my(vw = VW(k, E1, E2, lam0, mu0));
  vw = vector(2, i, my(b = vw[i]); for(a = 1, 2*k, b = subst(b, eval(Str("p", a)), sum(j = 1, k, (1/cs[j])^a))); b);
  my(fo = heqs(k, r1, r2, lam0, mu0, M, Pv, vw));
  my(G = t^2 - 1, H = (t - r1)*(t - r2), f0 = 1 + G*(lam0*t + mu0), f = vector(k, j, f0 + cs[j]*G*H), F = prod(j = 1, k, f[j]));
  my(Km = prod(j = 1, k, (f[j] - 1)/G), Kp = prod(j = 1, k, (f[j] + 1)/H));
  my(Vd = (F - 1) \ (G*Km), Wd = (F + 1) \ (H*Kp), RV = (F - 1)/(G*Km), RW = (F + 1)/(H*Kp), dir = List());
  print("V formal = direct: ", vw[1] == Vd, "  W: ", vw[2] == Wd);
  foreach([r1, r2], r, my(D = subst(Vd, t, r + h + O(h^(M+1))) - subst(RV, t, r + h + O(h^(M+1)))); for(m = 0, M, listput(dir, polcoef(D, m, h))));
  foreach([1, -1], r, my(D = subst(Wd, t, r + h + O(h^(M+1))) - subst(RW, t, r + h + O(h^(M+1)))); for(m = 0, M, listput(dir, polcoef(D, m, h))));
  print("formal - direct (all orders 0..", M, "): ", fo - Vec(dir));
  print("residual sizes: ", apply(x -> if(x == 0, 0, round(log(abs(x))/log(10))), Vec(dir)));
}
