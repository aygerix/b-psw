\\ pc7_mono_hier: HIGHER-ORDER boundary equations for (2,2) quartic pencils (k odd).
\\ For a genuine family V = G^(k-1) Q_V exactly, Q_V = prod(1 + x_j) - prod x_j. Near a root r of H write
\\ f_j + 1 = d_j(h) = A(h) + c_j B(h)  (A = H L', B = H G at r + h), so f_j - 1 = -2 + d_j and
\\   Q_V = [1 + prod(1 - d_j)] / prod(2 - d_j) = 2^-k [exp(-S2) + exp(S1 - S2)],
\\   S1 = sum_j log(1 - d_j), S2 = sum_j log(1 - d_j/2): series in h whose coefficients involve the POSITIVE power sums
\\   P_m = sum_j c_j^m (P_0 = k).  Near s = +-1: f_j - 1 = e_j(h) = A' + c_j B' (A' = G L, B' = G H), f_j + 1 = 2 + e_j,
\\   Q_W = [1 + prod(1 + e_j)] / prod(2 + e_j).  Equations: V(r + h) - G^(k-1) Q_V = O(h^(M+1)) at r = r1, r2 and
\\   W(s + h) - H^(k-1) Q_W = O(h^(M+1)) at s = +-1. Orders 0, 1 are the 8 boundary equations (no P's); order m adds P_m.
read("q22_moments.gp");
VW(k, e1, e2, lam, mu) = {
  my(N = 2*k - 1, pv = vector(2*k, a, eval(Str("p", a))));
  my(X = mkx(N, w, lam, mu, e1, e2), Y = mky(N, w, lam, mu, e1, e2));
  my(LX = lg1(X, N), LY = lg1(-Y, N));
  my(SX = sum(n = 0, N - 1, topow(polcoef(LX, n, s), w, pv) * s^n), SY = sum(n = 0, N - 1, topow(polcoef(LY, n, s), w, pv) * s^n));
  my(EX = ex1(SX, N), EY = ex1(SY, N));
  my(BV = trunc((1 - s^2)^(k - 1) * EX, N), BW = trunc((1 - e1*s + e2*s^2)^(k - 1) * EY, N));
  [sum(n = 0, N - 1, polcoef(BV, n, s) * t^(2*k - 2 - n)), sum(n = 0, N - 1, polcoef(BW, n, s) * t^(2*k - 2 - n))]};
\\ sum_j F(a + c_j b) for F(z) = log(1 + sg z) as a series in h; a, b series in h with a(0) = b(0) = 0
plog(k, a, b, sg, M, Pv) = {
  my(r = O(h^(M + 1)));
  for(n = 1, M, my(term = sum(m = 0, n, binomial(n, m) * a^(n - m) * b^m * if(m == 0, k, Pv[m])));
    r += (-1)^(n + 1) * sg^n * term / n);
  r};
\\ local series; Pv = [P1, ..., PM]
locV(k, r, e1, e2, lam, mu, M, Pv) = {
  my(tt = r + h + O(h^(M + 1)), G = tt^2 - 1, H = tt^2 - e1*tt + e2, Lp = lam*tt + mu + lam*e1);
  my(A = H*Lp, B = H*G, S1 = plog(k, A, B, -1, M, Pv), S2 = plog(k, A/2, B/2, -1, M, Pv));
  G^(k - 1) * 2^(-k) * (exp(-S2) + exp(S1 - S2))};
locW(k, sv, e1, e2, lam, mu, M, Pv) = {
  my(tt = sv + h + O(h^(M + 1)), G = tt^2 - 1, H = tt^2 - e1*tt + e2, L = lam*tt + mu);
  my(A = G*L, B = G*H, T1 = plog(k, A, B, 1, M, Pv), T2 = plog(k, A/2, B/2, 1, M, Pv));
  H^(k - 1) * 2^(-k) * (exp(-T2) + exp(T1 - T2))};
\\ all 4(M+1) equations at the shape r1, r2 (e1 = r1 + r2, e2 = r1 r2); lam, mu passed (symbols or numbers)
heqs(k, r1, r2, lam, mu, M, Pv, VWv = 0) = {
  my(e1 = r1 + r2, e2 = r1*r2, vw = if(VWv, VWv, VW(k, e1, e2, lam, mu)), V = vw[1], W = vw[2], out = List());
  foreach([r1, r2], r, my(D = subst(V, t, r + h + O(h^(M + 1))) - locV(k, r, e1, e2, lam, mu, M, Pv));
    for(m = 0, M, listput(out, polcoef(D, m, h))));
  foreach([1, -1], sv, my(D = subst(W, t, sv + h + O(h^(M + 1))) - locW(k, sv, e1, e2, lam, mu, M, Pv));
    for(m = 0, M, listput(out, polcoef(D, m, h))));
  Vec(out)};
