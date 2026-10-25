\\ xn line: PROVED a priori bound on S = sum 1/a_i for k = 9 real cubic Williams families (lin = 0: 9 cubics;
\\ lin = 1: 8 cubics + 2t+1).  Ingredients (all PROVED, see xn notes):
\\  (a) even boundary equations Ve(-1) = 2^-8, Ve'(-1) = -2^-5 (Ve = (V + V~)/2, V~ = reflected family's V);
\\  (b) positivity: V(t) >= 2^-8 t^8 and V~(t) >= 2^-8 t^8 for all real t (even-subset identity, h_i/g_i > 0);
\\      with (a): Re(t) := (Ve(t) - 2^-8 t^8)/(t+1)^2 = (1-2^-8) rho8 + sum_n Tne rho_(8-n) >= 0,
\\      rho_d(t) = (t^d - (-1)^d - d (-1)^(d-1) (t+1))/(t+1)^2;
\\  (c) exact second order: Phi := (D2 + D2r)/2 = 2^-10 (e2(beta) + e2(gamma)) <= PhiBar, beta_i, gamma_i in (0,18)
\\      (linear member: beta = gamma = 2);  Phi = 2 Re(-1);
\\  (d) B0 = sum v_i^2 >= 0.
\\ For t0 < -1: Re(t0) - lam Phi = q2 S^2 + qB B0 + q1 S + q0 with lam chosen to cancel A2; if q2 < 0, qB <= 0,
\\ then 0 <= Re(t0) gives q2 S^2 + q1 S + q0 + lam PhiBar >= 0, i.e. S <= root.
rho(d, t) = (t^d - (-1)^d - d*(-1)^(d-1)*(t+1))/(t+1)^2;
{
  foreach([0, 1], lin,
    my(Te, Phi, PhiBar);
    if(lin == 0,
      Te = [0, 0, S, -3/2*S, 7/4*S + A2 - 2*B0, S^2/2 - 15/8*S - 5/2*A2 + 11/2*B0];
      Phi = S^2 - 209/4*S - 11*A2 + 23*B0 + 1785/32;
      PhiBar = 2*binomial(9, 2)*18^2/1024,
      Te = [1/2, 0, S, -S, S + A2 - 2*B0, S^2/2 - S - 2*A2 + 9/2*B0];
      Phi = S^2 - 40*S - 10*A2 + 21*B0 + 1113/32;
      PhiBar = 2*(binomial(8, 2)*18^2 + 8*18*2)/1024);
    \\ sanity: 2 Re(-1) == Phi
    my(Re(t0) = (1 - 2^-8)*subst(rho(8, 'x), 'x, t0) + sum(n = 1, 6, Te[n]*subst(rho(8 - n, 'x), 'x, t0)));
    print("lin = ", lin, ":  2 Re(-1) - Phi = ", simplify(2*Re(-1) - Phi), "   PhiBar = ", PhiBar, " = ", PhiBar*1.);
    my(best = [0, 10^9]);
    forstep(t0num = -400, -11, 1, my(t0 = t0num/40, R0 = Re(t0));
      my(lam = polcoef(R0, 1, A2)/polcoef(Phi, 1, A2));
      if(lam <= 0, next);
      my(D = R0 - lam*Phi, q2 = polcoef(D, 2, S), qB = polcoef(D, 1, B0), q1 = polcoef(polcoef(D, 1, S), 0, B0), q0 = polcoef(polcoef(D, 0, S), 0, B0));
      if(polcoef(D, 1, A2) != 0, error("A2 not cancelled"));
      if(q2 >= 0 || qB > 0, next);
      my(Q = q2*'y^2 + q1*'y + q0 + lam*PhiBar, rts = polrootsreal(Q), Smax = vecmax(rts));
      if(Smax < best[2], best = [t0, Smax]));
    print("   best t0 = ", best[1], ":  S <= ", best[2]);
    \\ exact certificate at the best t0: recompute and print the quadratic
    my(t0 = best[1], R0 = Re(t0), lam = polcoef(R0, 1, A2)/polcoef(Phi, 1, A2), D = R0 - lam*Phi);
    print("   lambda = ", lam, ";  Re(t0) - lambda Phi = ", D);
    my(Q = polcoef(D, 2, S)*'y^2 + polcoef(polcoef(D, 1, S), 0, B0)*'y + polcoef(polcoef(D, 0, S), 0, B0) + lam*PhiBar);
    print("   bound polynomial (must be >= 0): ", Q, "   -> S <= ", vecmax(polrootsreal(Q))));
}
