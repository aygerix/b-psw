\\ m9v: evaluate the INDEPENDENT (Euclidean-quotient) boundary equations at the Krawczyk centre of m9_kcert_9_1 (sanity check only)
\p 60
Vq(k, lin, L, refl) = {my(N = t^(k-1), D = 1);
  for(i = 1, #L, my(z = L[i], u = if(refl, 1 - z[1], z[1]), v = z[2]);
    N *= v*(2*t + 1) + t*(t + 1)*(t + u); D *= v*(2*t) + t*(t + 1)*(t + u));
  if(lin, N *= 2*t + 1; D *= 2*t);
  divrem(N, D, t)[1]};
bndL(k, lin, L) = {my(c0 = 2^(1-k), c1 = (1-k)*2^(1-k), V = Vq(k, lin, L, 0), Vr = Vq(k, lin, L, 1));
  [subst(V, t, -1) - c0, subst(deriv(V, t), t, -1) - c1, subst(Vr, t, -1) - c0, subst(deriv(Vr, t), t, -1) - c1]};
{
  my(fixv = [29087857/50000000, 128975363/100000000, 16755043/50000000, -83766877/100000000, 399093/2500000, 79580323/100000000, 6350827/50000000, 10615079/25000000, 4121457/50000000, 89485417/25000000, 1843907/12500000, 966509/10000000]);
  my(cen = [1.82953149302044206865672848354433378183, 4.2344194351053724574030595731241311242, 1.77272152572783112375930024201608557615, 2.26648178560688801734576671779062568235]);
  my(sel = [4, 6, 13, 14], fix = [1,2,3,5,7,8,9,10,11,12,15,16], X = vector(16));
  for(j = 1, 12, X[fix[j]] = fixv[j]); for(j = 1, 4, X[sel[j]] = cen[j]);
  my(L = vector(8, i, [X[i], X[8 + i]]));
  print("members [u, v]: ", L);
  print("independent boundary residuals at the certified point (8 cubics + 2t+1): ", bndL(9, 1, L));
  print("perturbed control (v1 += 1e-6): ", bndL(9, 1, concat([[L[1][1], L[1][2] + 1e-6]], L[2..8])));
  my(Mv(a, b) = sum(i = 1, 8, L[i][2]^a * (L[i][1] - 1/2)^b));
  print("M1_1 = ", Mv(1,1), "   S = ", Mv(1,0), "   (S - 111/8) M11 - (9/2 M13 - 19 M21) = ", (Mv(1,0) - 111/8)*Mv(1,1) - (9/2*Mv(1,3) - 19*Mv(2,1)));
}
