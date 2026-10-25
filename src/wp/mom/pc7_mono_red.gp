\\ pc7_mono_red: REDUCED order-M system. D_V := V - G^6 (deg 8) is triangular in p (t^(8-n) coefficient = c_n lam^n p_(n+1) + lower),
\\ so off the slice (lam != 0) the V-side equations V == T_V mod H^(M+1) are solved by D_V = R_V(P) + H^(M+1) A(t)
\\ (R_V = (T_V - G^6) mod H^(M+1), A of degree 8 - 2(M+1), free coefficients a0..), and p is recovered by back-substitution.
\\ The W-side equations then read  D_W(p) == R_W(P) mod G^(M+1): 2(M+1) equations in P1..PM, a's, shape.
read("pc7_mono_hs.gp");
redsys(k, M, e1, e2, lam, mu, Pv) = {
  my(vw = VW(k, e1, e2, lam, mu), V = vw[1], W = vw[2]);
  my(G = t^2 - 1, H = t^2 - e1*t + e2, L = lam*t + mu, Lp = lam*t + mu + lam*e1);
  my(QV = TQ(k, Lp, G, -1, M, Pv), QW = TQ(k, L, H, 1, M, Pv), HM = H^(M + 1), GM = G^(M + 1));
  my(TV = sum(n = 0, M, (G^(k - 1) * QV[n + 1] % HM) * H^n) % HM, TW = sum(n = 0, M, (H^(k - 1) * QW[n + 1] % GM) * G^n) % GM);
  my(RV = (TV - G^(k - 1)) % HM, RW = (TW - H^(k - 1)) % GM);
  my(na = 2*k - 2 - 4 - 2*(M + 1) + 1, av = vector(max(na, 0), i, eval(Str("a", i - 1))));
  my(DVs = RV + if(na > 0, HM * sum(i = 1, na, av[i] * t^(i - 1)), 0));
  my(DV = V - G^(k - 1), DW = W - H^(k - 1), dV = 2*k - 6, pv = vector(2*k - 1, a, eval(Str("p", a))), sol = vector(2*k - 1));
  \\ back-substitution: coefficient of t^(dV - n) determines p_(n+1)
  my(cur = DV);
  for(n = 0, dV, my(c = polcoef(cur, dV - n, t), q = pv[n + 1], co = deriv(c, q));
    if(co == 0, error("not triangular at n = ", n));
    if(deriv(co, q) != 0, error("nonlinear in p_", n + 1));
    my(val = (polcoef(DVs, dV - n, t) - subst(c, q, 0)) / co);
    sol[n + 1] = val; cur = subst(cur, q, val));
  if(cur != DVs, error("back-substitution check failed"));
  my(DWs = DW); for(a = 1, 2*k - 1, DWs = subst(DWs, pv[a], sol[a]));
  my(Rr = (DWs - RW) % GM);
  [vector(2*(M + 1), i, polcoef(Rr, i - 1, t)), sol, av]};
