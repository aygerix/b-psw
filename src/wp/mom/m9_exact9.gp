\\ m9 line, k = 9: EXACT feasible points of the four boundary equations (real cubic Williams families).
\\ Construction ("self-dual ansatz"):
\\   group A : nA asymmetric members (w_i, v_i), w_i rational != 0, v_i = lam d_i rational;
\\   pairs   : nP reflected pairs (+-w, v), rational;
\\   group B : self-dual members u = 1/2 (w = 0): nBf free rational ones + 3 whose v's are the roots of a cubic.
\\ Odd equations Z0-, Z1- are linear in S = M1_0 with coefficients multiple of M1_1 and do not involve M_{a,0} otherwise
\\ (checked below), so a combination Phi = sum_A phi(w_i, v_i) is S-free.  phi(w, .) is a cubic in v without constant
\\ term: the TANGENT line at v = 0 (sum w_i phi_1(w_i) d_i = 0) meets Phi = 0 again at a RATIONAL lam.
\\ Then Z0- = 0 gives S (rational), the even equations are linear in (M2_0, M3_0) given S, and group B must realise
\\ the power sums p_a = M_{a,0} - (known part), a = 1, 2, 3 (cubic with 3 distinct real roots > 1/32).
\\ Exact verification: the boundary equations eqs(m, lin) (brute force, cubic_moments.gp) vanish identically in the
\\ splitting algebra of the cubic (b3 = e1 - b1 - b2, reduce mod C2(b1, b2), then mod C(b1)).
read("m9_lib.gp");
Zs(m, lin) = {my(r = mform(m, lin), Z = r[1], Zr = [reflect(Z[1]), reflect(Z[2])]);
  [(Z[1] + Zr[1])/2, (Z[2] + Zr[2])/2, (Z[1] - Zr[1])/2, (Z[2] - Zr[2])/2]};
\\ substitute members (list of [w, v]) into moment polynomial; M_{a,0} gets + PP[a] (group-B power sums, symbolic)
msub(P, mem, PP) = {my(R = P); foreach(syms(P), z, my(st = strsplit(Str(z), "_"), a = eval(strjoin(Vec(st[1])[2..#Vec(st[1])])), b = eval(st[2]));
  R = subst(R, z, sum(i = 1, #mem, mem[i][2]^a * mem[i][1]^b) + if(b == 0 && a <= #PP, PP[a], 0))); R};
rr(lo, hi, den) = lo + random(floor((hi - lo) * den) + 1) / den;
\\ one attempt; returns 0 or [members (exact), cubic C]
try9(m, lin, ZZ, phi, nA, nP, nBf) = {
  my(w = vector(nA), d = vector(nA), lam, A = List(), mem, S, P1, P2, P3, e1, e2, e3, C, Q, sol);
  if(nA + 2*nP + nBf + 3 != m, error("count"));
  for(i = 1, nA, w[i] = rr(-5/2, 5/2, 8); while(w[i] == 0, w[i] = rr(-5/2, 5/2, 8)); d[i] = rr(1/16, 2, 16));
  \\ tangency at v = 0: sum w_i phi1(w_i) d_i = 0, solve for d[1]
  my(c1 = vector(nA, i, polcoef(subst(phi, ww, w[i]), 1, vv)), c2 = vector(nA, i, polcoef(subst(phi, ww, w[i]), 2, vv)),
     c3 = vector(nA, i, polcoef(subst(phi, ww, w[i]), 3, vv)));
  if(c1[1] == 0, return(0));
  d[1] = -sum(i = 2, nA, c1[i]*d[i]) / c1[1]; if(d[1] <= 0, return(0));
  my(q2 = sum(i = 1, nA, c2[i]*d[i]^2), q3 = sum(i = 1, nA, c3[i]*d[i]^3)); if(q3 == 0, return(0));
  lam = -q2/q3; if(lam <= 0, return(0));
  for(i = 1, nA, my(v = lam*d[i]); if((abs(w[i]) + 1/2)^2 >= 8*v, return(0)); listput(A, [w[i], v]));
  mem = Vec(A);
  for(j = 1, nP, my(v = rr(1/24, 2, 48), wmax = sqrtint(floor(64*(8*v)))/8 - 1/2, wp);
    if(wmax <= 1/16, return(0)); wp = rr(1/16, wmax, 16); if((wp + 1/2)^2 >= 8*v, return(0));
    mem = concat(mem, [[wp, v], [-wp, v]]));
  for(j = 1, nBf, mem = concat(mem, [[0, rr(1/24, 2, 48)]]));
  \\ odd equation Z0- (linear in S through PP[1])
  my(z0m = msub(ZZ[3], mem, [p1, p2, p3]), z1m = msub(ZZ[4], mem, [p1, p2, p3]));
  if(poldegree(z0m, p1) != 1, return(0));
  P1 = simplify(-polcoef(z0m, 0, p1) / polcoef(z0m, 1, p1));
  if(type(P1) != "t_INT" && type(P1) != "t_FRAC", error("P1 not rational"));
  if(simplify(subst(z1m, p1, P1)) != 0, error("Phi not satisfied: construction bug"));
  my(z0p = subst(msub(ZZ[1], mem, [p1, p2, p3]), p1, P1), z1p = subst(msub(ZZ[2], mem, [p1, p2, p3]), p1, P1));
  my(Mx = simplify([polcoef(z0p, 1, p2), polcoef(z0p, 1, p3); polcoef(z1p, 1, p2), polcoef(z1p, 1, p3)]));
  if(poldegree(z0p, p2) > 1 || poldegree(z0p, p3) > 1 || poldegree(z1p, p2) > 1 || poldegree(z1p, p3) > 1, error("not linear"));
  if(matdet(Mx) == 0, return(0));
  sol = matsolve(Mx, -[subst(subst(z0p, p2, 0), p3, 0), subst(subst(z1p, p2, 0), p3, 0)]~);
  P2 = simplify(sol[1]); P3 = simplify(sol[2]);
  if(type(P2) == "t_POL" || type(P3) == "t_POL", error("P2/P3 not rational"));
  if(P1 <= 3/32, return(0));
  e1 = P1; e2 = (P1^2 - P2)/2; e3 = (P1^3 - 3*P1*P2 + 2*P3)/6;
  C = x^3 - e1*x^2 + e2*x - e3;
  if(poldisc(C) == 0 || subst(C, x, 1/32) == 0, return(0));
  if(polsturm(C, [1/32, +oo]) != 3, return(0));
  for(j = 1, #mem, if(mem[j][1] == 0 && subst(C, x, mem[j][2]) == 0, return(0)));
  [mem, C]};
\\ exact verification of a found configuration with the brute-force boundary equations
verify9(m, lin, mem, C) = {
  my(E = eqs(m, lin), nk = #mem, e1 = -polcoef(C, 2, x), Cb1 = subst(C, x, b1), C2, out = vector(4));
  C2 = simplify((subst(C, x, b2) - subst(C, x, b1)) / (b2 - b1));
  for(j = 1, 4, my(e = E[j]);
    for(i = 1, nk, e = subst(subst(e, eval(Str("u", i)), 1/2 + mem[i][1]), eval(Str("v", i)), mem[i][2]));
    for(i = 1, 3, e = subst(subst(e, eval(Str("u", nk + i)), 1/2), eval(Str("v", nk + i)), [b1, b2, b3][i]));
    e = subst(e, b3, e1 - b1 - b2);
    e = divrem(e, C2, b2)[2]; e = divrem(e, Cb1, b1)[2];
    out[j] = e);
  out};
\\ inequality report: domain margins, distinctness
report9(mem, C) = {my(rt = polrootsreal(C), allm = concat(mem, vector(3, j, [0, rt[j]])), mg = vector(#allm, i, sqrt(8*allm[i][2]) - 1/2 - abs(allm[i][1])));
  print("  margins sqrt(8v) - 1/2 - |w| (numerical): ", mg);
  print("  min margin ", vecmin(mg), "; min pairwise distance ", vecmin(concat(vector(#allm, i, vecmin(concat([10^9], vector(#allm - i, j, abs(allm[i][1] - allm[i+j][1]) + abs(allm[i][2] - allm[i+j][2]))))))));
  print("  cubic roots (group B v's): ", rt)};
