\\ pc7_mono_bil: E-LINEAR (bilinear) form of the order-M equations. E_m = e_m(c_1..c_k) (elementary symmetric in the c_j).
\\ V-side (eta = H, K_j = L' + c_j G): prod(1 - eta K_j) = sum_n (-eta)^n e_n(K), e_n(K) = sum_m binom(k - m, n - m) L'^(n-m) G^m E_m,
\\   V * prod(1 - H K_j/2) == 2^-k G^(k-1) (1 + prod(1 - H K_j)) mod H^(M+1)   [k odd]
\\ W-side (eta = G, J_j = L + c_j H):  W * prod(1 + G J_j/2) == 2^-k H^(k-1) (1 + prod(1 + G J_j)) mod G^(M+1).
\\ Both are LINEAR in (E_1..E_M) and linear in the coefficients of V, W.
read("pc7_mono_hs.gp");
eK(k, n, A, B, Ev) = sum(m = 0, n, binomial(k - m, n - m) * A^(n - m) * B^m * if(m == 0, 1, Ev[m]));
bil(k, e1, e2, lam, mu, M, Ev, VWv = 0) = {
  my(vw = if(VWv, VWv, VW(k, e1, e2, lam, mu)), V = vw[1], W = vw[2]);
  my(G = t^2 - 1, H = t^2 - e1*t + e2, L = lam*t + mu, Lp = lam*t + mu + lam*e1, HM = H^(M + 1), GM = G^(M + 1));
  my(eV = vector(M, n, eK(k, n, Lp, G, Ev)), eW = vector(M, n, eK(k, n, L, H, Ev)));
  my(DV = 1 + sum(n = 1, M, (-1/2)^n * H^n * eV[n]), NV = 1 + sum(n = 1, M, (-1)^n * H^n * eV[n]));
  my(DW = 1 + sum(n = 1, M, (1/2)^n * G^n * eW[n]), NW = 1 + sum(n = 1, M, G^n * eW[n]));
  my(RV = (V * DV - 2^(-k) * G^(k - 1) * (1 + NV)) % HM, RW = (W * DW - 2^(-k) * H^(k - 1) * (1 + NW)) % GM);
  concat(vector(2*M + 2, i, polcoef(RV, i - 1, t)), vector(2*M + 2, i, polcoef(RW, i - 1, t)))};
