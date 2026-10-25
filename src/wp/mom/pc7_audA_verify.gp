\\ pc7_audA_verify: independent audit of the min-member obstruction for (2,2) quartic pencils.
\\ Independent of pc7_mono_thm.gp: parametrise by the roots r1, r2 directly, and compute root sums by RESULTANTS
\\ (not by remainder traces).
default(realprecision, 80);
default(parisizemax, 1000000000);
t; c; r1; r2; Y;
\\ sum over the roots of the quadratic Q (in t) of A(t), via Res_t(Q, Y - A) = lc(Q)^deg(A) * prod (Y - A(z_i))
rsumres(A, Q) = {my(R = polresultant(Q, Y - A, t)); -polcoef(R, 1, Y) / polcoef(R, 2, Y)};
{
  my(e1 = r1 + r2, e2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2), lam = 2*e1/Dl, mu = 2*(1 - e1^2 + e2)/Dl);
  my(G = t^2 - 1, H = (t - r1)*(t - r2), L = lam*t + mu, Lp = lam*t + mu + lam*e1, f0 = 1 + G*L);
  print("Dl - (1 - e1^2 + 2 e2 + e2^2) = ", simplify(Dl - (1 - e1^2 + 2*e2 + e2^2)));
  print("(I1) f0 + 1 - H L' = ", simplify(f0 + 1 - H*Lp), ";  f0(1) - 1 = ", simplify(subst(f0, t, 1) - 1), ";  f0(-1) - 1 = ", simplify(subst(f0, t, -1) - 1));
  my(f = f0 + c*G*H, Km = c*H + L, Kp = c*G + Lp);
  print("(I1) f - 1 - G K- = ", simplify(f - 1 - G*Km), ";  f + 1 - H K+ = ", simplify(f + 1 - H*Kp));
  print("deg_t K- = ", poldegree(Km, t), ", lc = ", polcoef(Km, 2, t), ";  deg_t K+ = ", poldegree(Kp, t), ", lc = ", polcoef(Kp, 2, t));
  my(Sm = rsumres(G*H, Km), Sp = rsumres(G*H, Kp));
  print("(I4, resultant) (S- - S+)/2 - 4/c = ", simplify((Sm - Sp)/2 - 4/c));
  my(Fm = rsumres(f0, Km), Fp = rsumres(f0, Kp));
  print("(I3, resultant) sum_{K+} f0 - sum_{K-} f0 = ", simplify(Fp - Fm));
  print("(I5) poldisc(K+) * c^-2 - (lam^2 w^2 - 4(mu + lam e1) w + 4) at w = 1/c: ", simplify(poldisc(Kp, t)/c^2 - (lam^2/c^2 - 4*(mu + lam*e1)/c + 4)));
  print("(I5') mu + lam e1 - 2(1 + e2)/Dl = ", simplify(mu + lam*e1 - 2*(1 + e2)/Dl));
}
\\ numeric: random COMPLEX shapes and complex c (the identity is algebraic), direct roots
{
  my(worst = 0.);
  for(trial = 1, 200,
    my(R1 = (random(2001) - 1000)/337. + I*(random(2001) - 1000)/541., R2 = (random(2001) - 1000)/419. + I*(random(2001) - 1000)/613.);
    my(E1 = R1 + R2, E2 = R1*R2, D = (1 - R1^2)*(1 - R2^2), la = 2*E1/D, m = 2*(1 - E1^2 + E2)/D, cc = (random(2001) - 1000)/97. + I*(random(2001)-1000)/89.);
    if(abs(cc) < 1e-3, next);
    my(gh(z) = (z^2 - 1)*(z - R1)*(z - R2));
    my(zm = polroots(cc*(x - R1)*(x - R2) + la*x + m), zp = polroots(cc*(x^2 - 1) + la*x + m + la*E1));
    my(dd = (gh(zm[1]) + gh(zm[2]))/2 - (gh(zp[1]) + gh(zp[2]))/2 - 4/cc);
    worst = max(worst, abs(dd) / (1 + abs(4/cc))));
  print("numeric (I4), 200 random complex shapes and c: max rel err = ", worst);
}
\\ Step 1 end-to-end on random real pencils (mixed signs allowed, k = 7 and 9): F at the roots of K-+_j
{
  my(worst = 0.);
  for(trial = 1, 20, my(k = if(trial % 2, 7, 9));
    my(R1 = (random(1999) - 999)/1000, R2 = (random(1999) - 999)/1000);
    if(R1 == R2 || abs(R1) == 1 || abs(R2) == 1, next);
    my(E1 = R1 + R2, E2 = R1*R2, D = (1 - R1^2)*(1 - R2^2), la = 2*E1/D, m = 2*(1 - E1^2 + E2)/D);
    my(cs = vector(k, j, (random(4001) - 2000)/317), G = x^2 - 1, H = (x - R1)*(x - R2), f0 = 1 + G*(la*x + m));
    my(F = prod(j = 1, k, f0 + cs[j]*G*H));
    for(j = 1, k, if(cs[j] == 0, next);
      foreach(polroots(cs[j]*H + la*x + m), z, my(u = subst(G*H, x, z));
        worst = max(worst, abs(subst(F, x, z) - prod(l = 1, k, if(l == j, 1, 1 + (cs[l] - cs[j])*u)))));
      foreach(polroots(cs[j]*G + la*x + m + la*E1), z, my(u = subst(G*H, x, z));
        worst = max(worst, abs(subst(F, x, z) - (-1)^k*prod(l = 1, k, if(l == j, 1, 1 - (cs[l] - cs[j])*u)))))));
  print("Step 1 product formulas, 20 random real pencils (k = 7, 9, mixed signs): max abs err = ", worst);
}
\\ lemma control with some d_l = 0 and widely spread d_l
{
  my(worst = -oo);
  for(trial = 1, 400, my(mm = 2 + trial % 8, d = vector(mm, l, if(random(4) == 0, 0, 10.^(random(9) - 4) * (1 + random(1000)/1000.))));
    if(vecmax(d) == 0, next);
    my(P = prod(l = 1, mm, 1 + d[l]*x) - 1); if(P == 0, next);
    my(rts = polroots(P / x));
    for(i = 1, #rts, worst = max(worst, real(rts[i]))));
  print("lemma control (zeros allowed): max Re u over nonzero roots: ", worst);
}
