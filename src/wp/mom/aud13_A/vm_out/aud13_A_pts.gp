\\ aud13_A_pts.gp -- (a) negative controls for the identity checks; (b) cross-check of the relaxation rows actually used
\\ by the certificates (Z0e, Z1e, CS, PHI, POS at t0 in {-19/5, 0, 1/5, 3/5, 4/5}) at the m9 boundary points.
\\ Uses only my forms aud13_A_forms_{0,1}.gp.
\p 60
t; e; z;
read("aud13_A_forms_0.gp"); read("aud13_A_forms_1.gp");
MS = [M1_0, M1_1, M1_2, M2_0, M1_3, M2_1, M1_4, M2_2, M3_0, M1_5, M2_3, M3_1];
MP = [[1,0],[1,1],[1,2],[2,0],[1,3],[2,1],[1,4],[2,2],[3,0],[1,5],[2,3],[3,1]];
ev(X, U, V) = substvec(X, MS, vector(#MS, i, sum(q = 1, #U, V[q]^MP[i][1]*(U[q] - 1/2)^MP[i][2])));
rowsat(FM, U, V) = {my(Z0 = ev(FM[1], U, V), Z1 = ev(FM[2], U, V), Phi = ev(FM[3], U, V), Te = FM[4],
    S = vecsum(V), M11 = sum(i = 1, #U, V[i]*(U[i] - 1/2)), A2 = sum(i = 1, #U, V[i]*(U[i] - 1/2)^2),
    pos = vector(5, j, my(t0 = [-19/5, 0, 1/5, 3/5, 4/5][j]); (1 - 2^-8)*t0^8 + sum(n = 1, 8, ev(Te[n], U, V)*t0^(8 - n))));
  [S, Z0, Z1, Phi, S*A2 - M11^2, pos]};
{
  foreach(["m9_sol_9_0.gp", "m9_sol_9_1.gp"], fn, read(fn);
    my(lin = CASE[2], r = rowsat(if(lin, FORMS1, FORMS0), U0, V0));
    print(fn, " lin = ", lin, ": S = ", r[1]*1., "  Z0e = ", r[2]*1., " Z1e = ", r[3]*1., "  Phi = ", r[4]*1., "  S*A2 - M11^2 = ", r[5]*1.);
    print("     POS (Ve(t0) - 2^-8 t0^8) at t0 = -19/5, 0, 1/5, 3/5, 4/5: ", precision(r[6]*1., 6)));
  \\ level-2 large-S points (u, v from the m9_l2cert_9.gp start vectors; only S and Phi are relevant: S is beyond the cover)
  my(L2 = [[0, [-1.5022143532841334, -4.123908001403858, 0.5139699491914088, -6.9488109387036845, 0.4557359591171056, 1.9898717712453684, 8.384115532862122, 0.0300074710797625, 0.2116870195358351],
                      [2.1953624588934777, 4.160359143007728, 0.1638235274005706, 12.93704174320486, 0.1504054295230128, 11.332431645079225, 13.503906596594524, 0.8344275287407846, 16.583303646561404]],
               [1, [7.383602524018461, -0.281220887363323, 3.7362686543302996, -6.113476957120817, -1.9894291170655451, 0.9158133556559549, -1.4477281105089814, -2.7740334229151604],
                      [10.187126685349437, 7.09708766558695, 3.420275946834393, 9.28262974093801, 1.3912245328824524, 13.036011292752544, 1.249721086956731, 4.659440794003997]]]);
  foreach(L2, c, my(lin = c[1], r = rowsat(if(lin, FORMS1, FORMS0), c[2], c[3]));
    my(Re75 = (( (1 - 2^-8)*(-7/5)^8 + sum(n = 1, 8, ev(if(lin, FORMS1, FORMS0)[4][n], c[2], c[3])*(-7/5)^(8 - n))) - r[2] - r[3]*(-2/5))/(4/25));
    print("level-2 start lin = ", lin, ": S = ", r[1]*1., " Z0e = ", r[2]*1., " Z1e = ", r[3]*1., " Phi = ", r[4]*1., "  Re(-7/5) = ", Re75*1.));
}
\\ (a) negative controls: perturbed second/third-order formulas must FAIL at a random member set
{
  my(k = 9, A = vector(9, i, random(300)/17 + 1/7), C = vector(9, i, (random(301) - 150)/13));
  my(fs = vector(9, i, 2*t + 1 + t*(t + 1)*(A[i]*t + C[i])), F = prod(i = 1, 9, fs[i]), gs = vector(9, i, (fs[i] - 1)/t), hs = vector(9, i, (fs[i] + 1)/(t + 1)));
  my(H = (F - 1)/(t*prod(i = 1, 9, gs[i])), Hs = subst(H, t, -1 + e) + O(e^4));
  my(be = vector(9, i, subst(hs[i], t, -1)), ga = vector(9, i, subst(gs[i], t, 0)), e2b = (vecsum(be)^2 - sum(i = 1, 9, be[i]^2))/2);
  my(E1 = -be/2, E2 = vector(9, i, be[i] + ga[i]/2 - 3 - be[i]^2/4), cr = vecsum(E1)*vecsum(E2) - sum(i = 1, 9, E1[i]*E2[i]));
  print("negctl c2 with 2^-k instead of 2^(-1-k): equal? ", polcoef(Hs, 2, e) == binomial(8, 2)*2^-8 + 2^-9*e2b);
  print("negctl c3 with (k-3)/4 in xn's form: equal? ", polcoef(Hs, 3, e) == 2^-8*(-binomial(8, 3) - 6/4*e2b + cr));
  print("positive c3 (xn form): equal? ", polcoef(Hs, 3, e) == 2^-8*(-binomial(8, 3) - 8/4*e2b + cr));
}
