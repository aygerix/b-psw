\\ pc7_mono_hs: SYMMETRIC form of the order-M (2,2)-pencil equations (k odd), in e1, e2 (no individual roots).
\\ V-side: f_j + 1 = H K+_j, K+_j = L' + c_j G; f_j - 1 = -2 + H K+_j. With eta := H,
\\   Q_V = [1 + prod(1 - eta K_j)] / prod(2 - eta K_j) = 2^-k [exp(-S2) + exp(S1 - S2)], S1 = sum log(1 - eta K_j), S2 = sum log(1 - eta K_j/2),
\\   sum_j K_j^n = sum_m binom(n, m) L'^(n - m) G^m P_m, P_m = sum c_j^m (P_0 = k).
\\   Genuine family: V = G^(k-1) Q_V exactly, hence V == T_V := G^(k-1) Q_V(eta = H) mod H^(M+1) for every M.
\\ W-side: f_j - 1 = G J_j, J_j = L + c_j H; f_j + 1 = 2 + G J_j; Q_W = [1 + prod(1 + eta J_j)] / prod(2 + eta J_j), eta := G;
\\   W == T_W := H^(k-1) Q_W(eta = G) mod G^(M+1).
\\ Equations: the 2(M+1) coefficients (t^0..t^(2M+1)) of (V - T_V) mod H^(M+1) and of (W - T_W) mod G^(M+1).
read("pc7_mono_hier.gp");
\\ truncated series in eta as vectors [c_0, ..., c_M] (entries: polynomials in t and parameters)
smul(a, b) = {my(M = #a - 1); vector(M + 1, n, sum(i = 0, n - 1, a[i + 1] * b[n - i]))};
sexp(a) = {my(M = #a - 1, r = vector(M + 1), P = vector(M + 1)); \\ a[1] = 0
  r[1] = 1; P[1] = 1; my(f = 1);
  for(m = 1, M, P = smul(P, a); f *= m; r += P / f); r};
\\ sum_j log(1 + sg * eta * K_j) with K_j = A + c_j B (A, B polynomials in t)
slog(k, A, B, sg, M, Pv) = {
  my(r = vector(M + 1));
  for(n = 1, M, r[n + 1] = (-1)^(n + 1) * sg^n * sum(m = 0, n, binomial(n, m) * A^(n - m) * B^m * if(m == 0, k, Pv[m])) / n);
  r};
TQ(k, A, B, sg, M, Pv) = {   \\ 2^-k [exp(-S2) + exp(S1 - S2)]: S1 = slog(sg), S2 = slog(sg/2)
  my(S1 = slog(k, A, B, sg, M, Pv), S2 = slog(k, A, B, sg/2, M, Pv));
  2^(-k) * (sexp(-S2) + sexp(S1 - S2))};
hsym(k, e1, e2, lam, mu, M, Pv, VWv = 0) = {
  my(vw = if(VWv, VWv, VW(k, e1, e2, lam, mu)), V = vw[1], W = vw[2]);
  my(G = t^2 - 1, H = t^2 - e1*t + e2, L = lam*t + mu, Lp = lam*t + mu + lam*e1);
  my(QV = TQ(k, Lp, G, -1, M, Pv), QW = TQ(k, L, H, 1, M, Pv));
  my(HM = H^(M + 1), GM = G^(M + 1));
  my(TV = sum(n = 0, M, (G^(k - 1) * QV[n + 1] % HM) * H^n) % HM, TW = sum(n = 0, M, (H^(k - 1) * QW[n + 1] % GM) * G^n) % GM);
  my(RV = (V - TV) % HM, RW = (W - TW) % GM);
  concat(vector(2*M + 2, i, polcoef(RV, i - 1, t)), vector(2*M + 2, i, polcoef(RW, i - 1, t)))};
