\\ m9v (adversarial verification of the m9 order-2 claim).  INDEPENDENT of m9_lib / m9_moments / cubic_moments:
\\ series via t_SER in a variable of highest priority, members built from (a, c) directly, V from exact polynomial division.
\\ Part A (exact, random rational points): r_0, r_1, r_2 of H - V at t = -1, H = t^(k-1)(F-1)/prod(f_j - 1), V = poly quotient;
\\   compare with  r_0 = 2^(1-k) - V(-1),  r_1 = (1-k)2^(1-k) - V'(-1),  r_2 = [e_2(s) - (2^k V''(-1) - 2(k-1)(k-2))]/2^(k+1);
\\   also V(t) = t^(k-1) T(1/t) and the coefficient of T'(-1) in V''(-1).
\\ Part B (exact, symbolic): T2+ - Red  ==  lam1 Z0+ + lam2 Z1+  identically in (u_i, v_i), lam constants,
\\   Red = the claimed reduced form; Z0+ = (T(-1) + T~(-1))/2 - 2^(1-k), Z1+ = (T'(-1) + T~'(-1))/2.
S = varhigher("S");
xs(u, v, k) = v*S^3 / (1 + (1+u)*S + (u + 2*v)*S^2 + O(S^k));
Tpol(U, W, k, lin) = {my(T = 1 + O(S^k)); for(i = 1, #U, T *= 1 + xs(U[i], W[i], k)); if(lin, T *= 1 + S/2); truncate(T)};
Tder(P) = [subst(P, S, -1), subst(deriv(P, S), S, -1), subst(deriv(deriv(P, S), S), S, -1)];
e2(x) = (vecsum(x)^2 - sum(i = 1, #x, x[i]^2))/2;
partA(m, lin, ntr) = {my(k = m + lin, bad = 0);
  for(tr = 1, ntr,
    my(A = vector(m, i, random(97) + 1)/ (random(13) + 1), C = vector(m, i, (random(201) - 100)/(random(17) + 1)));
    my(f = vector(m, i, 2*t + 1 + t*(t+1)*(A[i]*t + C[i])));
    if(lin, f = concat(f, [2*t + 1]));
    my(F = prod(i = 1, k, f[i]), Pi = prod(i = 1, k, f[i] - 1), N = t^(k-1)*(F - 1), qr = divrem(N, Pi), V = qr[1], R = qr[2]);
    my(z = 'z, ser = subst(R/Pi, t, -1 + z + O(z^4)), r = vector(3, j, polcoef(ser, j - 1, z)));
    my(V0 = subst(V, t, -1), V1 = subst(deriv(V, t), t, -1), V2 = subst(deriv(deriv(V, t), t), t, -1));
    my(sv = vector(k, i, subst(deriv(f[i], t), t, -1)));
    my(c0 = 2^(1-k));
    my(ok0 = r[1] == c0 - V0, ok1 = r[2] == (1-k)*c0 - V1, ok2 = r[3] == (e2(sv) - (2^k*V2 - 2*(k-1)*(k-2)))/2^(k+1));
    \\ V from the Laurent part, with u = c/a, v = 1/a
    my(U = vector(m, i, C[i]/A[i]), W = vector(m, i, 1/A[i]), P = Tpol(U, W, k, lin), VL = sum(j = 0, k - 1, polcoef(P, j, S)*t^(k-1-j)));
    my(okV = (VL == V), Td = Tder(P), cT1 = if(Td[2] == 0, "T'=0", (V2 - (k-1)*(k-2)*Td[1] - Td[3])/Td[2]));
    \\ s_i = h_i(-1) = 2 + (1-u)/v for cubics
    my(okS = vector(m, i, sv[i] == 2 + (1 - U[i])/W[i]) == vector(m, i, 1));
    if(!(ok0 && ok1 && ok2 && okV && okS), bad++; print("  FAIL trial ", tr, ": ", [ok0, ok1, ok2, okV, okS]));
    if(tr <= 2, print("  trial ", tr, ": r0,r1,r2 ok = ", [ok0, ok1, ok2], "  V==Laurent ", okV, "  coefficient of T'(-1) in V''(-1): ", cT1)));
  print("partA m=", m, " lin=", lin, ": failures ", bad, " of ", ntr)};
ev(E, p, U, W) = {my(Q = E, m = #U); for(i = 1, m, Q = subst(subst(Q, U[i], p[i]), W[i], p[m + i])); Q};
partB(m, lin, Red) = {my(k = m + lin, U = vector(m, i, eval(Str("u", i))), W = vector(m, i, eval(Str("v", i))), Ur = vector(m, i, 1 - U[i]));
  my(P = Tpol(U, W, k, lin), Pr = Tpol(Ur, W, k, lin), d = Tder(P), dr = Tder(Pr));
  my(Z0 = (d[1] + dr[1])/2 - 2^(1-k), Z1 = (d[2] + dr[2])/2, T2p = (d[3] + dr[3])/2);
  my(Sm = vecsum(W), M12 = sum(i = 1, m, W[i]*(U[i] - 1/2)^2), M20 = sum(i = 1, m, W[i]^2));
  my(M11 = sum(i = 1, m, W[i]*(U[i] - 1/2)), M14 = sum(i = 1, m, W[i]*(U[i] - 1/2)^4), M22 = sum(i = 1, m, W[i]^2*(U[i] - 1/2)^2), M30 = sum(i = 1, m, W[i]^3));
  my(Rd = Red(Sm, M12, M20, M11, M14, M22, M30), D = T2p - Rd);
  \\ solve lam from two random rational points, then test the identity symbolically
  my(pts = vector(3, j, vector(2*m, i, (random(2001) - 1000)/(random(50) + 1))));
  my(Mt = matrix(3, 2, r, c, ev(if(c == 1, Z0, Z1), pts[r], U, W)), rhs = vector(3, r, ev(D, pts[r], U, W))~);
  my(lam = matsolve(Mt[1..2, ], rhs[1..2]));
  my(res = D - lam[1]*Z0 - lam[2]*Z1);
  print("partB m=", m, " lin=", lin, ": lam = ", lam, "   symbolic residual D - lam.Z == 0 ? ", res == 0);
  if(res != 0, print("   residual variables: ", variables(res), "  sample: ", ev(res, pts[3], U, W)))};
default(parisizemax, 2500000000);
setrand(20261025);
partA(9, 0, 6); partA(8, 1, 6); partA(11, 0, 4); partA(10, 1, 4); partA(3, 0, 4);
partB(9, 0, (S1, M12, M20, M11, M14, M22, M30) -> S1^2 - 209/4*S1 - 11*M12 + 23*M20 + 1785/32);
partB(8, 1, (S1, M12, M20, M11, M14, M22, M30) -> S1^2 - 40*S1 - 10*M12 + 21*M20 + 1113/32);
partB(11, 0, (S1, M12, M20, M11, M14, M22, M30) -> 83/4*S1^2 + 2*M12*S1 - 4*M20*S1 - 2561/16*S1 + M11^2 - 169/2*M12 + 785/4*M20 - 13*M14 + 81*M22 - 56*M30 + 46035/512);
partB(10, 1, (S1, M12, M20, M11, M14, M22, M30) -> 71/4*S1^2 + 2*M12*S1 - 4*M20*S1 - 112*S1 + M11^2 - 67*M12 + 631/4*M20 - 12*M14 + 75*M22 - 52*M30 + 27603/512);
quit;
