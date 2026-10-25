\\ pc7_mono_thm: exact verification of the identities behind the MIN-MEMBER argument for real (2,2) quartic pencils.
\\ Setting: G = t^2 - 1, H = t^2 - e1 t + e2, Dl = 1 - e1^2 + 2 e2 + e2^2, lam = 2 e1/Dl, mu = 2(1 - e1^2 + e2)/Dl,
\\ L = lam t + mu, L' = lam t + mu + lam e1, f0 = 1 + G L; member f = f0 + c G H; K- = c H + L, K+ = c G + L'.
\\ (I1) f0 + 1 = H L', so f - 1 = G K-, f + 1 = H K+.
\\ (I2) at roots z of K-: G H(z) = (1 - f0(z))/c;  at roots z' of K+: G H(z') = -(1 + f0(z'))/c   (as remainders mod K-+).
\\ (I3) [sum over roots of K+ of f0] - [sum over roots of K- of f0] = 4, identically in c.
\\ (I4) hence  mean_{K-}(GH) - mean_{K+}(GH) = 4/c.
\\ (I5) reality: K+ has non-real roots iff lam^2 w^2 - 4(mu + lam e1) w + 4 < 0 (w = 1/c), and mu + lam e1 = 2(1 + e2)/Dl.
default(parisizemax, 500000000);
t; c; e1; e2; w;
rsum(A, Q) = {  \\ sum of A over the two roots of the quadratic Q (in t): trace of multiplication by A mod Q
  my(R = A % Q, a = polcoef(Q, 2, t), b = polcoef(Q, 1, t), q0 = polcoef(Q, 0, t));
  polcoef(R, 1, t) * (-b/a) + 2 * polcoef(R, 0, t)};
{
  my(Dl = 1 - e1^2 + 2*e2 + e2^2, lam = 2*e1/Dl, mu = 2*(1 - e1^2 + e2)/Dl);
  my(G = t^2 - 1, H = t^2 - e1*t + e2, L = lam*t + mu, Lp = lam*t + mu + lam*e1, f0 = 1 + G*L);
  my(Km = c*H + L, Kp = c*G + Lp);
  print("(I1) f0 + 1 - H L' = ", simplify(f0 + 1 - H*Lp));
  print("(I1') f0(+-1) = ", simplify([subst(f0, t, 1), subst(f0, t, -1)]), "   f0 mod H + 1 = ", simplify((f0 + 1) % H));
  print("(I2-) (G H - (1 - f0)/c) mod K- = ", simplify((G*H - (1 - f0)/c) % Km));
  print("(I2+) (G H + (1 + f0)/c) mod K+ = ", simplify((G*H + (1 + f0)/c) % Kp));
  print("(I3) sum_{K+} f0 - sum_{K-} f0 = ", simplify(rsum(f0, Kp) - rsum(f0, Km)));
  print("(I4) mean_{K-}(GH) - mean_{K+}(GH) - 4/c = ", simplify(rsum(G*H, Km)/2 - rsum(G*H, Kp)/2 - 4/c));
}
{
  my(Dl = 1 - e1^2 + 2*e2 + e2^2, lam = 2*e1/Dl, mu = 2*(1 - e1^2 + e2)/Dl, G = t^2 - 1, Lp = lam*t + mu + lam*e1);
  my(Kp = (1/w)*G + Lp, D = simplify(poldisc(Kp, t) * w^2));
  print("(I5) disc_t(K+(c = 1/w)) * w^2 = ", D, "   equals lam^2 w^2 - 4 (mu + lam e1) w + 4: ", simplify(D - (lam^2*w^2 - 4*(mu + lam*e1)*w + 4)) == 0);
  print("(I5') mu + lam e1 - 2(1 + e2)/Dl = ", simplify(mu + lam*e1 - 2*(1 + e2)/Dl));
}
\\ numeric control of the lemma: for random d_l > 0 (6 of them), every nonzero root u of prod(1 + d_l u) = 1 has Re u < 0
{
  my(worst = -oo);
  for(trial = 1, 300, my(d = vector(6, l, random(10^6)/10^5 + 10^-3), P = prod(l = 1, 6, 1 + d[l]*x) - 1, rts = polroots(P / x));
    for(i = 1, #rts, worst = max(worst, real(rts[i]))));
  print("lemma control: max Re(u) over nonzero roots, 300 random d's: ", worst);
}
\\ numeric control of (I4) at a random real shape and random c
{
  my(r1 = -0.31, r2 = 0.57, E1 = r1 + r2, E2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2), lam = 2*E1/Dl, mu = 2*(1 - E1^2 + E2)/Dl, cc = 1.7);
  my(G(z) = z^2 - 1, H(z) = z^2 - E1*z + E2);
  my(zm = polroots(cc*(x^2 - E1*x + E2) + lam*x + mu), zp = polroots(cc*(x^2 - 1) + lam*x + mu + lam*E1));
  my(mm = sum(i = 1, 2, G(zm[i])*H(zm[i]))/2, mp = sum(i = 1, 2, G(zp[i])*H(zp[i]))/2);
  print("numeric (I4): mean- - mean+ = ", real(mm - mp), "   4/c = ", 4/cc, "   roots K- ", zm, "  K+ ", zp);
}
