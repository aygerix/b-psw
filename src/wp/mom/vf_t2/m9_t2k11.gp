\\ m9 line: k = 11 order-2 quantity T2+ reduced modulo the even boundary equations (eliminate the two top moments
\\ M1_6 and M4_0, then report which moments remain).  Uses m9_t2.gp functions (identity of T''(-1) checked there).
read("m9_lib.gp");
T2form(k, Lam) = {my(Lp = sum(n = 1, k - 1, Lam[n] * s^n), E = 1, P = 1);
  for(j = 1, k, P = trunc(P * Lp, k) / j; if(P == 0, break); E += P);
  subst(deriv(deriv(E, s), s), s, -1)};
T2mom(m, lin) = {my(k = m + lin, Ls = Lsym(k), el = ellv(1/2 + w, v, k), Z = T2form(k, Ls));
  for(n = 1, k - 1, Z = subst(Z, Ls[n], lin * ellin(k)[n] + if(n >= 3, lin2mom(el[n]), 0))); Z};
T2brute(m, lin, U, W) = {my(k = m + lin, G = 1);
  for(i = 1, m, G = trunc(G * (1 + xser(U[i], W[i], k)), k));
  if(lin, G = trunc(G * (1 + s/2), k));
  subst(deriv(deriv(G, s), s), s, -1)};
{
  foreach([[11, 0], [10, 1]], c, my(m = c[1], lin = c[2], k = m + lin, T2 = simplify(T2mom(m, lin)), T2r = reflect(T2), r = mform(m, lin), Z = r[1], Zr = [reflect(Z[1]), reflect(Z[2])]);
    my(U = vector(m, i, eval(Str("u", i))), W = vector(m, i, eval(Str("v", i))));
    print("==== k = ", k, ", m = ", m, ", lin = ", lin);
    print("identity check T''(-1): ", simplify(T2brute(m, lin, U, W) - evalmom(T2, m)));
    my(Zp = [(Z[1] + Zr[1])/2, (Z[2] + Zr[2])/2], T2p = (T2 + T2r)/2, T2m = (T2 - T2r)/2);
    print("T2+ = ", T2p);
    foreach([[M1_6, M4_0], [M1_6, M2_4], [M1_6, M3_2]], pr, my(a = pr[1], b = pr[2]);
      my(A = simplify([polcoef(Zp[1], 1, a), polcoef(Zp[1], 1, b); polcoef(Zp[2], 1, a), polcoef(Zp[2], 1, b)]));
      if(matdet(A) == 0, print("  pair ", pr, ": singular"); next);
      my(sol = matsolve(A, -simplify([subst(subst(Zp[1], a, 0), b, 0), subst(subst(Zp[2], a, 0), b, 0)]~)));
      my(T2pr = simplify(subst(subst(T2p, a, sol[1]), b, sol[2])));
      print("  eliminating ", pr, ": T2+ reduced = ", T2pr);
      print("     remaining moments: ", syms(T2pr)));
    print("T2- = ", T2m));
}
