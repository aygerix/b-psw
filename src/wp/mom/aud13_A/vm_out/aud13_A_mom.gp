\\ aud13_mom.gp -- independent audit of candidate Theorem 13 (k = 9), part (i)/(ii) + row export.
\\ New code; does not read any xn_* or m9_* script.  It READS m9_mform.out (text) only to compare.
\\ Derivation: x_i = 1/(f_i - 1) = v s^3/(1 + (1+u) s + (u+2v) s^2), s = 1/t, u = w + 1/2;
\\ T(s) = prod(1 + x_i) mod s^9 = 1 + e1(x) + e2(x) (e3 = O(s^9)), times (1 + s/2) for 2t+1;
\\ V(t) = sum_n T_n t^(8-n) = polynomial quotient of (F - 1) by t prod g_i.
s; w; v; t; e; z;
Msym(p, j) = eval(Str("M", p, "_", j));
MLIST = concat(vector(4, p, vector(7, j, Msym(p, j - 1))));
tomom(P) = {my(r = 0, Q, c); P = simplify(P);
  for(p = 0, poldegree(P, v), Q = polcoef(P, p, v);
    for(j = 0, poldegree(Q, w), c = polcoef(Q, j, w); if(c, if(p == 0, error("p = 0 monomial")); r += c*Msym(p, j))));
  r};
MPJ = concat(vector(4, p, vector(7, j, [p, j - 1])));
refl(X) = substvec(X, MLIST, vector(#MLIST, i, (-1)^MPJ[i][2]*MLIST[i]));
substr(st, a, b) = if(b < a, "", concat(Vec(st)[a..b]));
\\ per-member series
XS = v*s^3/(1 + (1 + (w + 1/2))*s + ((w + 1/2) + 2*v)*s^2) + O(s^9);
X1 = sum(n = 3, 8, tomom(polcoef(XS, n))*s^n);
X2 = sum(n = 6, 8, tomom(polcoef(XS^2, n))*s^n);
Tser(lin) = {my(T = (1 + X1 + (X1^2 - X2)/2)*(1 + s/2)^lin); vector(9, n, polcoef(T, n - 1, s))};
Vpol(lin) = {my(T = Tser(lin)); sum(n = 0, 8, T[n + 1]*t^(8 - n))};
\\ evaluate moments at a configuration
momval(U, VV) = vector(#MLIST, i, my(p = MPJ[i][1], j = MPJ[i][2]); sum(q = 1, #U, VV[q]^p*(U[q] - 1/2)^j));
evalm(X, U, VV) = substvec(X, MLIST, momval(U, VV));
\\ exact family data from members
famFG(U, VV, lin) = {my(F = 1, G = 1);
  for(i = 1, #U, my(a = 1/VV[i], c = U[i]/VV[i], f = 2*t + 1 + t*(t + 1)*(a*t + c)); F *= f; G *= (f - 1)/t);
  if(lin, F *= 2*t + 1; G *= 2); [F, G]};
quotV(U, VV, lin) = {my(FG = famFG(U, VV, lin), dr = divrem(FG[1] - 1, t*FG[2])); dr[1]};
randcfg(m) = {my(U = vector(m), VV = vector(m)); for(i = 1, m, VV[i] = random(400)/37 + 1/31; U[i] = (random(201) - 100)/53); [U, VV]};
{
for(lin = 0, 1,
  my(m = 9 - lin, Vs = Vpol(lin), ok = 1);
  print("==== lin = ", lin, " (", m, " cubics)");
  \\ (A) exact check of the moment form of V against polynomial division (random rational members, not in D needed)
  for(trial = 1, 6, my(c = randcfg(m), q = quotV(c[1], c[2], lin), d = q - evalm(Vs, c[1], c[2])); if(d != 0, ok = 0; print("MISMATCH ", d)));
  print("(A) V moment form == quotient of (F-1)/(t prod g) at 6 random rational configs: ", ok);
  \\ (B) compare T_n with exp(sum Lam_n s^n) from m9_mform.out
  my(L = readstr("m9_mform.out"), blk = 0, Lam = vector(8), Zp = vector(2));
  for(i = 1, #L, my(l = L[i]);
    if(#l > 10 && substr(l, 1, 9) == "==== k = ", blk = (l == Str("==== k = 9, cubic members m = ", m, ", linear member: ", lin)); next);
    if(!blk, next);
    my(q = strsplit(l, " = "));
    if(#q == 2, my(h = q[1]);
      if(#h >= 7 && substr(h, 3, 6) == "Lam_", Lam[eval(substr(h, 7, #h))] = eval(q[2]));
      if(h == "Z0+", Zp[1] = eval(q[2])); if(h == "Z1+", Zp[2] = eval(q[2]))));
  if(lin, Lam[1] = 1/2; Lam[2] = -1/8);
  my(Ls = sum(n = 1, 8, Lam[n]*s^n), Ex = sum(j = 0, 8, Ls^j/j!), Tm = vector(9, n, polcoef(Ex, n - 1, s)), Tme = Tser(lin));
  print("(B) T_n(mine) - T_n(exp Lam, m9_mform.out), n = 0..8: ", vector(9, n, simplify(Tme[n] - Tm[n])));
  \\ (C) even part, boundary equations
  my(Ve = (Vs + refl(Vs))/2, Q = Ve - 2^-8*t^8, Z0 = subst(Q, t, -1), Z1 = subst(deriv(Q, t), t, -1));
  print("(C) Z0e - Z0+(m9) = ", simplify(Z0 - Zp[1]), ";  Z1e + 8 Z0+ + Z1+ (m9) = ", simplify(Z1 + 8*Zp[1] + Zp[2]));
  my(Rt = (Q - Z0 - Z1*(t + 1)), dv = divrem(Rt, (t + 1)^2));
  if(dv[2] != 0, error("Re not polynomial"));
  my(Re = simplify(dv[1]));
  print("    variables of Re(t) (after simplify): ", variables(Re));
  print("    Re(t) = ", Re);
  my(Phi = 2*subst(Re, t, -1), Phi2 = subst(deriv(deriv(Q, t), t), t, -1));
  my(PhiClaim = if(lin == 0, M1_0^2 - 209/4*M1_0 - 11*M1_2 + 23*M2_0 + 1785/32, M1_0^2 - 40*M1_0 - 10*M1_2 + 21*M2_0 + 1113/32));
  print("(D) Phi = 2 Re(-1) = Q''(-1)? ", Phi == Phi2, ";  Phi - T2+ claim = ", simplify(Phi - PhiClaim));
  print("    Phi = ", Phi);
  \\ (E) S-identity at t0 = -7/5
  my(t0 = -7/5, R0 = subst(Re, t, t0), lam = polcoef(R0, 1, M1_2)/polcoef(Phi, 1, M1_2), D = simplify(R0 - lam*Phi));
  my(lamC = if(lin == 0, 59/110, 27/50), DC = if(lin == 0, -(2/55)*M1_0^2 - (788/125)*M1_0 - (2/55)*M2_0 + 39391023/1100000,
                                                    -(1/25)*M1_0^2 - (678/125)*M1_0 - (1/25)*M2_0 + 1294023/50000));
  print("(E) Re(-7/5) = ", R0);
  print("    lambda (A2-cancelling) = ", lam, "  claimed ", lamC, ";  D - claimed = ", simplify(D - DC), ";  D = ", D);
  \\ direct check: Re(-7/5) = (Ve(-7/5) - 2^-8 (7/5)^8 - Z0 - Z1 (t0+1))/(t0+1)^2 as moment identity
  print("    direct (Q(t0) - Z0 - Z1(t0+1))/(t0+1)^2 - Re(t0) = ", simplify((subst(Q, t, t0) - Z0 - Z1*(t0 + 1))/(t0 + 1)^2 - R0));
  my(PhiBar = if(lin == 0, 2*binomial(9, 2)*18^2/2^10, 2*(binomial(8, 2)*18^2 + 2*8*18)/2^10));
  my(q2 = polcoef(D, 2, M1_0), q1 = polcoef(polcoef(D, 1, M1_0), 0, M2_0), qB = polcoef(D, 1, M2_0), q0 = polcoef(polcoef(D, 0, M1_0), 0, M2_0));
  print("    PhiBar = ", PhiBar, "; q2 = ", q2, " qB = ", qB, " q1 = ", q1, " q0 = ", q0);
  my(B = q2*'y^2 + q1*'y + q0 + lam*PhiBar, rts = polrootsreal(B));
  print("    bound poly ", B, " roots ", rts, "  => S <= ", vecmax(rts));
  \\ exact rational test of the bound at the cover ends
  my(Send = if(lin == 0, (57/8 + 1/2)^2/8 + 8/32, (13/2 + 1/2)^2/8 + 7/32));
  print("    cover end S = ", Send, " = ", Send*1., ";  bound poly at cover end = ", subst(B, 'y, Send), " (must be < 0)");
  \\ scan t0 for robustness: best bound over t0 in [-10, -1.025]
  my(best = [0, 10^9]);
  forstep(tn = -400, -41, 1, my(tt = tn/40, RR = subst(Re, t, tt), la = polcoef(RR, 1, M1_2)/polcoef(Phi, 1, M1_2));
    if(la <= 0, next); my(DD = RR - la*Phi, a2 = polcoef(DD, 2, M1_0), aB = polcoef(DD, 1, M2_0));
    if(polcoef(DD, 1, M1_2) != 0, error("A2 not cancelled"));
    if(a2 >= 0 || aB > 0, next);
    my(BB = a2*'y^2 + polcoef(polcoef(DD, 1, M1_0), 0, M2_0)*'y + polcoef(polcoef(DD, 0, M1_0), 0, M2_0) + la*PhiBar, rr = polrootsreal(BB));
    if(#rr && vecmax(rr) < best[2], best = [tt, vecmax(rr)]));
  print("    scan t0 in [-10, -1.025] step 1/40: best t0 = ", best[1], " S <= ", best[2]);
  \\ (F) export moment forms for the relaxation rows: Z0, Z1, Phi, Te_n (n = 1..8) -- even parts
  my(Te = vector(8, n, simplify((polcoef(Ve, 8 - n, t)))));
  write(Str("aud13_A_forms_", lin, ".gp"), "FORMS", lin, " = [", Z0, ", ", Z1, ", ", Phi, ", ", Te, "];");
  print("    exported forms; Te = ", Te);
);
}
