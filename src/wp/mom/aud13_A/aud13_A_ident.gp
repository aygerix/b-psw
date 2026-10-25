\\ aud13_A_ident.gp -- audit of candidate Theorem 13, item (ii) and the positivity/Taylor identities. New code.
\\ For arbitrary real members (not a family) H := (F - 1)/(t prod g) is a rational function; for a family H = V.
\\ Checks, exactly over Q at random rational members, for k = 5, 7, 9, 11, with and without 2t+1:
\\  (P)  H = 2^(1-k) sum_j t^(k-1-2j) (t+1)^(2j) e_(2j)(h_i/g_i)          (positivity identity)
\\  (T0..T3) Taylor coefficients of H at t = -1:
\\       c0 = 2^(1-k), c1 = -(k-1) 2^(1-k), c2 = C(k-1,2) 2^(1-k) + 2^(-1-k) e2(beta),
\\       c3 = 2^(1-k)[-C(k-1,3) - (k-1)/4 e2(beta) + sum_{i != j} e1_i e2_j],  e1 = -beta/2, e2 = beta + gamma/2 - 3 - beta^2/4
\\       (xn's form; my own form: c3 = 2^(1-k)[-C(k-1,3) - (k-3)/4 e2(beta) + sum_{i != j} r_i r'_j], r = h/g at -1)
\\  (H0) H(0) = sum_i prod_{j != i} 1/g_j(0)
\\  (RF) reflection: -f(-1-t) = 2t+1 + t(t+1)(a t + (a - c)), and beta <-> gamma
\\  (B)  2^10 * (e2(beta) + e2(gamma))/2^10 == (sum b)^2 + (sum om)^2 - sum b^2 - sum om^2, b = 2 + a/2, om = a w
t; e; a; c; z;
e2v(L) = {my(s1 = vecsum(L), s2 = sum(i = 1, #L, L[i]^2)); (s1^2 - s2)/2};
mem(A, C, lin) = {my(L = vector(#A, i, 2*t + 1 + t*(t + 1)*(A[i]*t + C[i]))); if(lin, L = concat(L, [2*t + 1])); L};
chk(k, lin, A, C) = {
  my(fs = mem(A, C, lin), F = prod(i = 1, #fs, fs[i]), gs = vector(#fs, i, (fs[i] - 1)/t), hs = vector(#fs, i, (fs[i] + 1)/(t + 1)));
  my(H = (F - 1)/(t*prod(i = 1, #gs, gs[i])), rs = vector(#fs, i, hs[i]/gs[i]));
  \\ (P) elementary symmetric functions of the r_i via generating polynomial in z
  my(Gz = prod(i = 1, #rs, 1 + z*rs[i]), Pid = 2^(1 - k)*sum(j = 0, (k - 1)/2, t^(k - 1 - 2*j)*(t + 1)^(2*j)*polcoef(Gz, 2*j, z)));
  my(okP = (simplify(H - Pid) == 0));
  \\ Taylor at -1
  my(Hs = subst(H, t, -1 + e) + O(e^4), cf = vector(4, n, polcoef(Hs, n - 1, e)));
  my(be = vector(#fs, i, subst(hs[i], t, -1)), ga = vector(#fs, i, subst(gs[i], t, 0)));
  my(E1 = vector(#fs, i, -be[i]/2), E2 = vector(#fs, i, be[i] + ga[i]/2 - 3 - be[i]^2/4));
  my(cross = vecsum(E1)*vecsum(E2) - sum(i = 1, #fs, E1[i]*E2[i]));
  my(r0 = vector(#fs, i, subst(rs[i], t, -1)), r1 = vector(#fs, i, subst(deriv(rs[i], t), t, -1)));
  my(cross2 = vecsum(r0)*vecsum(r1) - sum(i = 1, #fs, r0[i]*r1[i]));
  my(cl = [2^(1 - k), -(k - 1)*2^(1 - k), binomial(k - 1, 2)*2^(1 - k) + 2^(-1 - k)*e2v(be),
           2^(1 - k)*(-binomial(k - 1, 3) - (k - 1)/4*e2v(be) + cross)]);
  my(c3mine = 2^(1 - k)*(-binomial(k - 1, 3) - (k - 3)/4*e2v(be) + cross2));
  my(okT = vector(4, n, cf[n] == cl[n]), ok3 = (cf[4] == c3mine));
  \\ H(0)
  my(okH0 = (subst(H, t, 0) == sum(i = 1, #fs, prod(j = 1, #fs, if(j == i, 1, 1/ga[j])))));
  \\ reflection: reflected members, check beta <-> gamma
  my(fr = vector(#fs, i, -subst(fs[i], t, -1 - t)));
  my(Ar = concat(A, []), Cr = vector(#A, i, A[i] - C[i]), frc = mem(Ar, Cr, lin), okR = (fr == frc));
  my(ber = vector(#frc, i, subst((frc[i] + 1)/(t + 1), t, -1)), okRB = (ber == ga));
  \\ (B) b / om identity with w = c/a - 1/2
  my(bb = vector(#fs, i, if(i <= #A, 2 + A[i]/2, 2)), om = vector(#fs, i, if(i <= #A, A[i]*(C[i]/A[i] - 1/2), 0)));
  my(okB = (e2v(be) + e2v(ga) == vecsum(bb)^2 + vecsum(om)^2 - sum(i = 1, #bb, bb[i]^2) - sum(i = 1, #om, om[i]^2)));
  [okP, okT, ok3, okH0, okR, okRB, okB]};
{
  setrand(20261025);
  foreach([5, 7, 9, 11], k, for(lin = 0, 1,
    my(m = k - lin, allok = 1, out = 0);
    for(trial = 1, if(k <= 9, 4, 2),
      my(A = vector(m, i, random(300)/17 + 1/7), C = vector(m, i, (random(301) - 150)/13), r = chk(k, lin, A, C));
      if(r != [1, [1, 1, 1, 1], 1, 1, 1, 1, 1], allok = 0; print("k=", k, " lin=", lin, " FAIL ", r)));
    print("k = ", k, ", lin = ", lin, ": all identities (P, T0-T3 xn form, T3 own form, H(0), reflection, beta<->gamma, b/om) exact: ", allok)));
}
