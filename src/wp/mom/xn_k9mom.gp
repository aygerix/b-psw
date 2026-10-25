\\ xn line: k = 9 (9 cubics; 8 cubics + 2t+1) boundary equations and the second-order quantity in MOMENT form.
\\ w = u - 1/2 (reflection: w -> -w, v -> v).  Moments  M(p, j) = sum_i v_i^p w_i^j  (symbols Mp_j).
\\ T(s) = exp(Lam) mod s^k, Lam_n = sum_i ell_n(u_i, v_i) (+ linear member constant); V(t) = t^(k-1) T(1/t).
\\ E0 = V(-1) - 2^(1-k), E1 = V'(-1) + (k-1) 2^(1-k), D2 = V''(-1) - (k-1)(k-2) 2^(1-k)  [family: D2 = 2^(-k) e2(beta) > 0]
t; s; w; v; L1; L2; L3; L4; L5; L6; L7; L8; L9; L10;
trunc(P, n) = sum(j = 0, n - 1, polcoef(P, j, s) * s^j);
xser(u, vv, n) = {my(c = vector(n + 1)); c[1] = 1; if(n >= 1, c[2] = -(1 + u));
  for(j = 2, n, c[j+1] = -(1 + u)*c[j] - (u + 2*vv)*c[j-1]);
  trunc(vv * s^3 * sum(j = 0, n, c[j+1]*s^j), n + 1)};
ellv(u, vv, k) = {my(X = xser(u, vv, k), L = 0, P = 1);
  for(j = 1, k, P = trunc(P * X, k); if(P == 0, break); L += (-1)^(j+1) * P / j);
  vector(k - 1, n, polcoef(L, n, s))};
ellin(k) = vector(k - 1, n, (-1)^(n+1) / (n * 2^n));
\\ T as polynomial in formal L1..L(k-1)
Tpoly(k) = {my(Ls = vector(k - 1, n, eval(Str("L", n))), Lp = sum(n = 1, k - 1, Ls[n] * s^n), E = 1, P = 1);
  for(j = 1, k, P = trunc(P * Lp, k) / j; if(P == 0, break); E += P); E};
\\ V and derivatives at -1 from T(s) = sum T_n s^n: V(t) = sum T_n t^(k-1-n)
Vders(k, T) = {my(V = sum(n = 0, k - 1, polcoef(T, n, s) * t^(k-1-n)));
  [subst(V, t, -1) - 2^(1-k), subst(deriv(V, t), t, -1) + (k-1)*2^(1-k), subst(deriv(deriv(V, t), t), t, -1) - (k-1)*(k-2)*2^(1-k)]};
\\ replace a single-member polynomial in (w, v) by moment symbols
tomom(P) = {my(r = 0); for(p = 0, poldegree(P, v), my(Q = polcoef(P, p, v));
   for(j = 0, poldegree(Q, w), my(cf = polcoef(Q, j, w)); if(cf != 0, r += cf * eval(Str("M", p, "_", j))))); r};
xn_k9mom_main() = {
  my(k = 9, ell = ellv(w + 1/2, v, k), ellr = ellv(1/2 - w, v, k), T = Tpoly(k), D = Vders(k, T));
  for(n = 3, k - 1, print("ell_", n, "(w, v) = ", ell[n]));
  foreach([0, 1], lin,
    print("\n==== k = 9, lin = ", lin, " (", k - lin, " cubic members)");
    foreach([0, 1], refl,
      my(Lv = vector(k - 1, n, tomom(if(refl, ellr[n], ell[n])) + lin * ellin(k)[n]));
      my(Dm = vector(3, r, my(e = D[r]); for(n = 1, k - 1, e = subst(e, eval(Str("L", n)), Lv[n])); e));
      print(if(refl, "reflected", "family"), ":");
      print("  E0 = ", Dm[1]); print("  E1 = ", Dm[2]); print("  D2 = ", Dm[3])));
}
\\ (library: read("xn_k9mom.gp"); xn_k9mom_main() prints the k = 9 moment forms of E0, E1, D2)
