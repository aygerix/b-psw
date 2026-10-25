\\ m9 line: moment form of the four boundary equations for real cubic Williams families, k = 9, 11
\\ (cases: k cubics; k-1 cubics + 2t+1).  Builds on cubic_moments.gp (xser, Vpoly, eqs).
\\
\\ Reformulation used here (elementary): with T(s) = prod_i (1 + x_i(s)) mod s^k (V(t) = t^(k-1) T(1/t)) and k odd,
\\   V(-1) = T(-1),  V'(-1) = -(k-1) T(-1) - T'(-1),
\\ so the boundary equations are  T(-1) = 2^(1-k),  T'(-1) = 0  (and the same for the reflected family).
\\ Moment form: log T = sum_n Lam_n s^n mod s^k, Lam_n = sum_i ell_n(u_i, v_i) (+ the linear member's constant),
\\ ell_n = coefficient of s^n in log(1 + x(u, v; s)).  So T(-1), T'(-1) are POLYNOMIALS in Lam_3..Lam_(k-1)
\\ (Lam_1, Lam_2 only from the linear member: constants).  The identity is checked exactly against eqs().
read("cubic_moments.gp");
\\ ell(u, v, k): vector [ell_1, ..., ell_(k-1)] (coefficients of log(1 + x) mod s^k)
ellv(u, v, k) = {my(X = xser(u, v, k), L = 0, P = 1);
  for(j = 1, k, P = trunc(P * X, k); if(P == 0, break); L += (-1)^(j+1) * P / j);
  vector(k - 1, n, polcoef(L, n, s))};
ellin(k) = vector(k - 1, n, (-1)^(n+1) / (n * 2^n));
\\ T(-1), T'(-1) as polynomials in symbols Lam = [L1..L(k-1)]
Tform(k, Lam) = {my(Lp = sum(n = 1, k - 1, Lam[n] * s^n), E = 1, P = 1);
  for(j = 1, k, P = trunc(P * Lp, k) / j; if(P == 0, break); E += P);
  [subst(E, s, -1), subst(deriv(E, s), s, -1)]};
Lsym(k) = vector(k - 1, n, eval(Str("L", n)));
\\ identity check: eqs(m, lin) = [E1, E2, E1r, E2r] with E1 = V(-1) - c0, E2 = V'(-1) - c1
checkid(m, lin) = {
  my(k = m + lin, c0 = 2^(1-k), Ls = Lsym(k), TF = Tform(k, Ls), E = eqs(m, lin), res = vector(4));
  for(refl = 0, 1,
    my(Lv = vector(k - 1, n, lin * ellin(k)[n]));
    for(i = 1, m, my(u = eval(Str("u", i)), v = eval(Str("v", i))); if(refl, u = 1 - u);
      Lv += ellv(u, v, k));
    my(T0 = TF[1], T1 = TF[2]);
    for(n = 1, k - 1, T0 = subst(T0, Ls[n], Lv[n]); T1 = subst(T1, Ls[n], Lv[n]));
    res[2*refl + 1] = E[2*refl + 1] - (T0 - c0);
    res[2*refl + 2] = E[2*refl + 2] - (-(k-1)*(T0 - c0) - T1));
  res};
