\\ pc7_cert_corner.gp -- blow-up of the D-system at the corner (r1, r2) -> (-1, 1):
\\ r1 = -1 + ep*a, r2 = 1 - ep*b, w_j = ep^SC * om_j (natural scale SC = 2: L ~ 1/ep, H - G ~ ep).
\\ Each D_i = kap_i + sum_m A_im p_m + sum_n c_in h_n + lam3_i e3 (3-body Hankel form); p_m(w) = ep^(SC m) p_m(om),
\\ h_n(w) = ep^(SC n) h_n(om), e3(w) = ep^(3 SC) e3(om).  Prints the ep-order of every coefficient (oo = identically 0)
\\ and the leading-order system.
SC = 2;
read("pc7_cert_shape.gp");
ep; a; b;
PREC = 14;
{
  my(K = SYM0, A = vector(8, i, SYM1[i] - SYM0[i]));
  my(C = matrix(8, 5, i, n, polcoef(polcoef(SYM2[i] - SYM0[i] - A[i] - subst(A[i], w, v), 1, w), n, v)));
  my(e = ep + O(ep^PREC), R1 = -1 + e*a, R2 = 1 - e*b, Dl = (1 - R1^2)*(1 - R2^2));
  my(SUB = [R1 + R2, R1*R2, 2*(R1 + R2)/Dl, 2*(1 - (R1 + R2)^2 + R1*R2)/Dl]);
  my(sb(P) = substvec(P, [e1, e2, lam, mu], SUB));
  my(ordl(P, sh) = if(P == 0, [oo, 0], my(Q = sb(P) * ep^sh, nu = valuation(Q, ep)); [nu, polcoef(Q, nu, ep)]));
  my(lam3 = [1,0,0,0,-1,0,0,0], LEAD = vector(8));
  for(i = 1, 8,
    my(T = [ordl(K[i], 0)]);
    for(m = 1, 9, T = concat(T, [ordl(polcoef(A[i], m, w), SC*m)]));
    for(n = 2, 6, T = concat(T, [ordl(C[i, n - 1], SC*n)]));
    T = concat(T, [if(lam3[i], [3*SC, lam3[i]], [oo, 0])]);
    my(o = vector(#T, j, T[j][1]), mn = vecmin(o));
    print("D", i, ": orders [kap, A1..A9, c2..c6, e3] = ", o, "   leading order ", mn);
    \\ leading-order expression in symbols P1..P9 (power sums of om), Hn (Hankel brackets), E3
    my(nm = concat(concat([1], vector(9, m, eval(Str("P", m)))), concat(vector(5, n, eval(Str("H", n + 1))), [E3])));
    LEAD[i] = sum(j = 1, #T, if(T[j][1] == mn, T[j][2] * nm[j], 0));
    print("   lead: ", LEAD[i]));
}
