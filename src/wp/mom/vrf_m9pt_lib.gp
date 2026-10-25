\\ vrf line: independent helpers (no cubic_moments.gp): V by exact polynomial division of F - 1 by t prod g_i.
t; r1; r2; r3; r4; Y; Y1; Y2; Y3; Y4;
fm(u, v) = 2*t + 1 + t*(t+1)*(t + u)/v;
gm(u, v) = 2 + (t+1)*(t + u)/v;
Vdiv(L) = {my(F = 1, G = t);
  for(i = 1, #L, my(f = fm(L[i][1], L[i][2]), g = gm(L[i][1], L[i][2]));
    if(f - 1 - t*g != 0, error("g mismatch")); F *= f; G *= g);
  my(qr = divrem(F - 1, G, t)); if(poldegree(qr[1], t) != #L - 1, error("deg V")); if(F - 1 - G*qr[1] - qr[2] != 0, error("division")); qr};
bnd(L) = {my(k = #L, qr = Vdiv(L), V = qr[1], c0 = 2^(1-k), c1 = (1-k)*2^(1-k));
  [subst(V, t, -1) - c0, subst(deriv(V, t), t, -1) - c1]};
refl(L) = vector(#L, i, [1 - L[i][1], L[i][2]]);
\\ exact rational enclosure of sqrt(x), x >= 0 rational
sqenc(x) = {my(S = 10^120, n = sqrtint(floor(x * S^2)), lo = n / S, hi = (n + 1) / S); if(!(lo^2 <= x && x <= hi^2), error("enc")); [lo, hi]};
