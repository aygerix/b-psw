\\ end-to-end numeric control of the member-level conditions used in the min-member argument (k = 7, random real pencil,
\\ NOT a Williams family): at every root z of K-_j,  F(z) = prod_{l != j} (1 + (c_l - c_j) GH(z));
\\ at every root z' of K+_j,  F(z') = (-1)^k prod_{l != j} (1 - (c_l - c_j) GH(z')).
\\ So Williams (F(z) = 1, F(z') = -1, k odd) <=> both products equal 1.
default(realprecision, 60);
{
  my(k = 7, r1 = -3/10, r2 = 11/20, E1 = r1 + r2, E2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2), lam = 2*E1/Dl, mu = 2*(1 - E1^2 + E2)/Dl);
  my(cs = [3/7, 1/2, 2/3, 5/6, 1, 4/3, 3/2]);
  my(G = x^2 - 1, H = x^2 - E1*x + E2, L = lam*x + mu, Lp = lam*x + mu + lam*E1, f0 = 1 + G*L);
  my(f = vector(k, j, f0 + cs[j]*G*H), F = prod(j = 1, k, f[j]), err = 0.);
  for(j = 1, k,
    foreach(polroots(cs[j]*H + L), z, my(u = subst(G*H, x, z));
      err = max(err, abs(subst(F, x, z) - prod(l = 1, k, if(l == j, 1, 1 + (cs[l] - cs[j])*u)))));
    foreach(polroots(cs[j]*G + Lp), z, my(u = subst(G*H, x, z));
      err = max(err, abs(subst(F, x, z) - (-1)^k * prod(l = 1, k, if(l == j, 1, 1 - (cs[l] - cs[j])*u))))));
  print("max |F - product formula| over all 28 roots: ", err);
  \\ the obstruction quantity for the min member j0 = 1: Re mean_{K-} GH - Re mean_{K+} GH = 4/c_1 > 0
  my(um = apply(z -> subst(G*H, x, z), polroots(cs[1]*H + L)), up = apply(z -> subst(G*H, x, z), polroots(cs[1]*G + Lp)));
  print("min member: mean GH over K- roots = ", real(vecsum(um))/2, ", over K+ roots = ", real(vecsum(up))/2, ", difference ", real(vecsum(um) - vecsum(up))/2, " vs 4/c = ", 4/cs[1]*1.);
}
