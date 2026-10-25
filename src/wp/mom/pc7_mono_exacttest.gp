\\ test: random rational pencil (k = 7): with V = (F-1) \ (G Km), W = (F+1) \ (H Kp) and the true E's, the residuals RV, RW must be
\\ minus the remainders of F -+ 1 modulo G Km, H Kp; and the top coefficients must vanish identically.
read("pc7_mono_exact.gp");
{
  my(k = 7, r1 = 1/3, r2 = -1/5, cs = [2, 3, 7/2, 5, 11/3, 13/7, 9/4]);
  my(E1 = r1 + r2, E2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2), lam0 = 2*E1/Dl, mu0 = 2*(1 - E1^2 + E2)/Dl);
  my(G = t^2 - 1, H = (t - r1)*(t - r2), f0 = 1 + G*(lam0*t + mu0), f = vector(k, j, f0 + cs[j]*G*H), F = prod(j = 1, k, f[j]));
  my(Km = prod(j = 1, k, (f[j] - 1)/G), Kp = prod(j = 1, k, (f[j] + 1)/H));
  my(Vd = (F - 1) \ (G*Km), Wd = (F + 1) \ (H*Kp), Ev = vector(k, m, polcoef(prod(j = 1, k, 1 + cs[j]*x), m, x)));
  my(R = exsys(k, E1, E2, lam0, mu0, Vd, Wd, Ev));
  my(remV = (F - 1) % (G*Km), remW = (F + 1) % (H*Kp));
  my(dirV = vector(4*k, i, -polcoef(remV, i - 1, t)), dirW = vector(4*k, i, -polcoef(remW, i - 1, t)));
  print("exact form - direct: ", R - concat(dirV, dirW));
  print("V monic deg 12: ", poldegree(Vd), " ", pollead(Vd), "   W: ", poldegree(Wd), " ", pollead(Wd));
}
