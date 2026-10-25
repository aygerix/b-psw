\\ pc7_mono_exact: EXACT bilinear form of the Williams condition for (2,2) pencils (no moments, no relaxation).
\\ With E_m = e_m(c_1..c_k):  prod_j K-_j = prod (L + c_j H) = sum_m E_m H^m L^(k-m),  prod_j K+_j = sum_m E_m G^m L'^(k-m),
\\ F = prod (f0 + c_j G H) = sum_m E_m (G H)^m f0^(k-m).  Williams (with the K's pairwise coprime and prime to G, H) <=>
\\   (VE)  G V sum_m E_m H^m L^(k-m) - sum_m E_m (GH)^m f0^(k-m) + 1 = 0,
\\   (WE)  H W sum_m E_m G^m L'^(k-m) - sum_m E_m (GH)^m f0^(k-m) - 1 = 0,
\\ for polynomials V, W of degree 2k - 2 (monic). Unknowns v_0..v_(2k-3), w_0..w_(2k-3), E_1..E_k: BILINEAR equations.
exsys(k, e1, e2, lam, mu, Vv = 0, Wv = 0, Ev = 0) = {
  my(G = t^2 - 1, H = t^2 - e1*t + e2, L = lam*t + mu, Lp = lam*t + mu + lam*e1, f0 = 1 + G*L, dg = 2*k - 2);
  my(V = if(Vv, Vv, t^dg + sum(i = 0, dg - 1, eval(Str("v", i)) * t^i)), W = if(Wv, Wv, t^dg + sum(i = 0, dg - 1, eval(Str("w", i)) * t^i)));
  my(E = if(Ev, Ev, vector(k, m, eval(Str("E", m)))), Em(m) = if(m == 0, 1, E[m]));
  my(F = sum(m = 0, k, Em(m) * (G*H)^m * f0^(k - m)));
  my(RV = G*V*sum(m = 0, k, Em(m) * H^m * L^(k - m)) - F + 1, RW = H*W*sum(m = 0, k, Em(m) * G^m * Lp^(k - m)) - F - 1);
  my(D = poldegree(RV, t));
  concat(vector(4*k, i, polcoef(RV, i - 1, t)), vector(4*k, i, polcoef(RW, i - 1, t)))};
