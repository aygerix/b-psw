\\ direct check of the quad = 1 variant: symmetric H (e1 = 0), members f0 (quadratic) and f0 + c_j G H, k = 5 total
read("q22_moments.gp");
{
  my(k = 5, r = 2/5, E2 = -r^2, cs = [2, 3, 7/2, 5]);
  my(Dl0 = (1 - r^2)^2, MU = 2*(1 + E2)/Dl0, AA = 2/(1 + E2));
  my(B = subst(bnd(k, 0, E2, 0, MU, 1), aa, AA));
  my(vals = vector(#B, i, my(b = B[i]); for(a = 1, 2*k, b = subst(b, eval(Str("p", a)), sum(j = 1, k - 1, (1/cs[j])^a))); b));
  my(G = t^2 - 1, H = t^2 + E2, f0 = 1 + MU*G);
  if(f0 - 1 != AA*G || f0 + 1 != AA*H, error("f0 setup"));
  my(f = concat([f0], vector(k - 1, j, f0 + cs[j]*G*H)), F = prod(j = 1, k, f[j]));
  my(Km = prod(j = 2, k, (f[j] - 1)/G), Kp = prod(j = 2, k, (f[j] + 1)/H));
  \\ lcm(f_j - 1) = AA G Km, so V := (F - 1)/(G Km) / AA has the formal normalisation (product over all k members)
  my(V = ((F - 1) \ (G*Km)) / AA, W = ((F + 1) \ (H*Kp)) / AA, c0 = 2^(1-k));
  my(R1 = (V - c0*G^(k-1)) % H, R2 = (deriv(V, t) - c0*deriv(G^(k-1), t)) % H, R3 = (W - c0*H^(k-1)) % G, R4 = (deriv(W, t) - c0*deriv(H^(k-1), t)) % G);
  my(dir = [polcoef(R1,0,t), polcoef(R1,1,t), polcoef(R2,0,t), polcoef(R2,1,t), polcoef(R3,0,t), polcoef(R3,1,t), polcoef(R4,0,t), polcoef(R4,1,t)]);
  print("quad: formal - direct: ", vals - dir, "   direct (nonzero check): ", dir);
}
