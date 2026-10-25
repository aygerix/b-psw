\\ m9 line: moment form of the order-2 quantity T''(-1) (k = 9, 11) and its reduction modulo the boundary equations.
\\ Order-2 Taylor equation at t = -1 on the boundary set:  e_2(s) = 2^k T''(-1),  s_i = h_i(-1) = 2 + (1-u_i)/v_i > 0.
\\ (V''(-1) = (k-1)(k-2) T(-1) + (2k-2) T'(-1) + T''(-1); e_2(s) = 2^k V''(-1) - 2(k-1)(k-2).)  Identity of the moment
\\ form of T''(-1) is checked exactly against the brute-force truncated product.
read("m9_lib.gp");
T2form(k, Lam) = {my(Lp = sum(n = 1, k - 1, Lam[n] * s^n), E = 1, P = 1);
  for(j = 1, k, P = trunc(P * Lp, k) / j; if(P == 0, break); E += P);
  subst(deriv(deriv(E, s), s), s, -1)};
T2mom(m, lin) = {my(k = m + lin, Ls = Lsym(k), el = ellv(1/2 + w, v, k), Z = T2form(k, Ls));
  for(n = 1, k - 1, Z = subst(Z, Ls[n], lin * ellin(k)[n] + if(n >= 3, lin2mom(el[n]), 0))); Z};
\\ brute force T''(-1) from the truncated product
T2brute(m, lin, U, W) = {my(k = m + lin, G = 1);
  for(i = 1, m, G = trunc(G * (1 + xser(U[i], W[i], k)), k));
  if(lin, G = trunc(G * (1 + s/2), k));
  subst(deriv(deriv(G, s), s), s, -1)};
{
  foreach([[9, 0], [8, 1]], c, my(m = c[1], lin = c[2], k = m + lin, T2 = simplify(T2mom(m, lin)), T2r = reflect(T2), r = mform(m, lin), Z = r[1], Zr = [reflect(Z[1]), reflect(Z[2])]);
    my(U = vector(m, i, eval(Str("u", i))), W = vector(m, i, eval(Str("v", i))));
    print("==== k = ", k, ", m = ", m, ", lin = ", lin);
    print("identity check T''(-1): ", simplify(T2brute(m, lin, U, W) - evalmom(T2, m)));
    my(Zp = [(Z[1] + Zr[1])/2, (Z[2] + Zr[2])/2], Zm = [(Z[1] - Zr[1])/2, (Z[2] - Zr[2])/2]);
    my(T2p = (T2 + T2r)/2, T2m = (T2 - T2r)/2);
    print("T2+ = ", T2p); print("T2- = ", T2m);
    \\ reduce T2+ modulo the even boundary equations: eliminate M1_4 and M3_0 (linear in Zp)
    my(sol = matsolve(simplify([polcoef(Zp[1], 1, M1_4), polcoef(Zp[1], 1, M3_0); polcoef(Zp[2], 1, M1_4), polcoef(Zp[2], 1, M3_0)]),
        -simplify([subst(subst(Zp[1], M1_4, 0), M3_0, 0), subst(subst(Zp[2], M1_4, 0), M3_0, 0)]~)));
    my(T2pr = simplify(subst(subst(T2p, M1_4, sol[1]), M3_0, sol[2])));
    print("T2+ reduced (M1_4, M3_0 eliminated with Z0+ = Z1+ = 0): ", T2pr);
    \\ odd: eliminate M1_5 and M3_1 (if possible)
    my(Am = simplify([polcoef(Zm[1], 1, M1_5), polcoef(Zm[1], 1, M3_1); polcoef(Zm[2], 1, M1_5), polcoef(Zm[2], 1, M3_1)]));
    print("odd elimination matrix det: ", matdet(Am));
    my(T2mr = T2m);
    if(matdet(Am) != 0, my(sm = matsolve(Am, -simplify([subst(subst(Zm[1], M1_5, 0), M3_1, 0), subst(subst(Zm[2], M1_5, 0), M3_1, 0)]~)));
       T2mr = simplify(subst(subst(T2m, M1_5, sm[1]), M3_1, sm[2])),
       \\ singular: use Z0- to eliminate M1_5 only, record the S-free relation
       T2mr = simplify(subst(T2m, M1_5, simplify(-subst(Zm[1], M1_5, 0)/polcoef(Zm[1], 1, M1_5)))));
    print("T2- reduced: ", T2mr);
    print("variables of reduced T2+: ", syms(T2pr), "   of reduced T2-: ", syms(T2mr)));
}
