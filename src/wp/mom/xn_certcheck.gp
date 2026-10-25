\\ xn line: INDEPENDENT exact checker (PARI/GP) for the k = 9 Farkas certificates written by xn_cover9.py.
\\ Independent of the Python code in three ways:
\\  (1) the moment forms (Es0, Es1, Phi, T_n^e) are recomputed here from the Laurent expansion (xn_k9mom.gp machinery),
\\      not taken from the Python constants;
\\  (2) the rows are rebuilt here from the certificate's row identifiers (kind, exact parameter);
\\  (3) positivity of Psi is checked by a different method: for v = (r+1/2)^2/8, W = theta r^2, v^2 Psi is
\\      A(r) theta^2 + B(r) theta + C(r); we check C > 0 and A + B + C > 0 on [0, Rb] by Sturm counts (polsturm),
\\      and 4AC - B^2 >= 0 wherever the vertex -B/(2A) lies in (0, 1) with A > 0 (root isolation + exact samples).
\\ A certificate [lin, Ra, Rb, y, z] proves: no k = 9 real cubic Williams family (lin = 1: with the member 2t+1)
\\ has S = sum_i 1/a_i in [S(Ra), S(Rb)], S(R) = (R + 1/2)^2/8 + (m-1)/32.
\\ Usage: gp -q xn_certcheck.gp, then  checkfile("xn_cert9_l0.gp", 0, Rneed)
read("xn_k9mom.gp");
Sv; rr; th; W; v;
\\ ---- moment forms (exact, recomputed) ----
mforms(lin) = {
  my(k = 9, ell = ellv(w + 1/2, v, k), ellr = ellv(1/2 - w, v, k), T = Tpoly(k), D = Vders(k, T));
  my(Dm = vector(2), Tn = vector(2));
  for(refl = 0, 1,
    my(Lv = vector(k - 1, n, tomom(if(refl, ellr[n], ell[n])) + lin * ellin(k)[n]));
    Dm[refl + 1] = vector(3, r, my(e = D[r]); for(n = 1, k - 1, e = subst(e, eval(Str("L", n)), Lv[n])); e);
    Tn[refl + 1] = vector(k - 1, n, my(e = polcoef(T, n, s)); for(j = 1, k - 1, e = subst(e, eval(Str("L", j)), Lv[j])); e));
  [vector(3, r, (Dm[1][r] + Dm[2][r])/2), vector(k - 1, n, (Tn[1][n] + Tn[2][n])/2)]};
\\ decompose X(moments) with M1_0 -> Sv:  [const(Sv), [coef of M1_2, M2_0, M1_4, M2_2, M3_0] (affine in Sv), coef of M1_1^2]
LINMOM = [M1_2, M2_0, M1_4, M2_2, M3_0];
decomp(X) = {
  my(Y = simplify(subst(X, M1_0, Sv)), cf = vector(5, i, simplify(polcoef(Y, 1, LINMOM[i]))), q = simplify(polcoef(Y, 2, M1_1)), c0 = Y);
  for(i = 1, 5, c0 = subst(c0, LINMOM[i], 0)); c0 = simplify(subst(c0, M1_1, 0));
  my(rec = c0 + sum(i = 1, 5, cf[i]*LINMOM[i]) + q*M1_1^2);
  if(Y != rec, error("unexpected moment structure: ", Y - rec));
  for(i = 1, 5, if(poldegree(cf[i], Sv) > 1 || #setminus(Set(variables(cf[i])), Set([Sv])) > 0, error("coef not affine in S")));
  if(#setminus(Set(variables(c0)), Set([Sv])) > 0 || #variables(q) > 0, error("stray variables"));
  if(poldegree(c0, Sv) > 2, error("const not quadratic"));
  if(type(q) != "t_INT" && type(q) != "t_FRAC", error("M1_1^2 coefficient not constant"));
  [c0, cf, q]};
\\ per-member values of the linear moments, and b, om^2
PMV = [v*W, v^2, v*W^2, v^2*W, v^3];
colof(dec) = sum(i = 1, 5, dec[2][i]*PMV[i]);       \\ column (function of W, v, Sv)
\\ ---- rows: [column(W, v, Sv), P-coefficient, rhs(Sv)] ----
rows_eq(lin, F) = {
  my(m = 9 - lin, E0 = decomp(F[1][1]), E1 = decomp(F[1][2]));
  [[1, 0, m], [v, 0, Sv], [colof(E0), 2*E0[3], -E0[1]], [colof(E1), 2*E1[3], -E1[1]]]};
row_ineq(lin, F, kind, prm) = {
  my(m = 9 - lin, PH = decomp(F[1][3]), b = 2 + 1/(2*v), om2 = W/v^2, blin = 2*lin, Mtot = m + lin);
  if(kind == "CS", return([-Sv/2*v*W, 1, 0]));
  if(kind == "PHI", return([-colof(PH), -2*PH[3], PH[1]]));
  if(kind == "N3L", my(B = prm); return([-colof(PH) + 2^-10*(2*B*b - b^2 - om2), -2*PH[3], PH[1] - 2^-10*(2*B*blin - B^2 - 4*lin)]));
  if(kind == "N3U", my(Bbar = 18*m + 2*lin);
    return([colof(PH) - 2^-9*((Bbar*b - b^2)/2 + (Mtot - 1)/2*om2), 2*PH[3], -PH[1] + 2^-9*(Bbar*blin - 4*lin)/2]));
  if(kind == "POS", my(t0 = prm, col = 0, pc = 0, rhs = (1 - 2^-8)*t0^8);
    for(n = 1, 8, my(Dn = decomp(F[2][n])); col -= colof(Dn)*t0^(8-n); pc -= 2*Dn[3]*t0^(8-n); rhs += Dn[1]*t0^(8-n));
    return([col, pc, rhs]));
  error("unknown row kind ", kind)};
\\ ---- positivity of q(theta) = A th^2 + B th + C on [0,1] for r in [0, R] ----
posint(p, R) = {   \\ p > 0 on [0, R] (exact)
  if(p == 0, return(0));
  if(poldegree(p, rr) <= 0, return(p > 0));
  (subst(p, rr, 0) > 0) && polsturm(p, [0, R]) == 0};
quadpos(Q, R) = {
  my(A = polcoef(Q, 2, th), B = polcoef(Q, 1, th), C = polcoef(Q, 0, th));
  if(!posint(C, R), return("C edge"));
  if(!posint(A + B + C, R), return("theta = 1 edge"));
  my(Dq = 4*A*C - B^2);
  if(posint(Dq, R), return(1));
  \\ partition by real roots of A, B, B + 2A, Dq in (0, R)
  my(pts = [0, R]);
  foreach([A, B, B + 2*A, Dq], p, if(p != 0 && poldegree(p, rr) > 0,
    foreach(polrootsreal(p, [0, R]), z, pts = concat(pts, [bestappr(z, 10^30)]))));
  pts = vecsort(pts, , 8);
  for(i = 1, #pts - 1, my(x = (pts[i] + pts[i+1])/2, a = subst(A, rr, x), b = subst(B, rr, x), d = subst(Dq, rr, x));
    if(a > 0 && b < 0 && b + 2*a > 0 && d <= 0, return(Str("interior vertex negative near r = ", x*1.))));
  \\ isolated roots of Dq with interior vertex: Dq = 0 gives min q = 0 >= 0, acceptable for Farkas (Psi >= 0 suffices)
  1};
SofR(R, m) = (R + 1/2)^2/8 + (m - 1)/32;
checkcert(c, F) = {
  my(lin = c[1], Ra = c[2], Rb = c[3], y = c[4], zl = c[5], m = 9 - lin, Sa = SofR(Ra, m), Sb = SofR(Rb, m));
  my(E = rows_eq(lin, F), I = vector(#zl, j, row_ineq(lin, F, zl[j][1], zl[j][2])));
  if(#y != 4, error("y length"));
  for(j = 1, #zl, if(zl[j][3] < 0, return("negative z")));
  my(col = sum(i = 1, 4, y[i]*E[i][1]) + sum(j = 1, #zl, zl[j][3]*I[j][1]));
  my(pc = sum(i = 1, 4, y[i]*E[i][2]) + sum(j = 1, #zl, zl[j][3]*I[j][2]));
  my(kap = sum(i = 1, 4, y[i]*E[i][3]) + sum(j = 1, #zl, zl[j][3]*I[j][3]));
  if(pc < 0, return("P column negative"));
  \\ kappa(S) < 0 on [Sa, Sb]
  kap = subst(kap, Sv, 'X);
  my(cand = [Sa, Sb]); if(poldegree(kap, 'X) == 2, my(xv = -polcoef(kap, 1, 'X)/(2*polcoef(kap, 2, 'X))); if(Sa < xv && xv < Sb, cand = concat(cand, [xv])));
  for(i = 1, #cand, if(subst(kap, 'X, cand[i]) >= 0, return("kappa not negative")));
  \\ Psi >= 0 on D_{Sb} at S = Sa and S = Sb (Psi affine in S)
  if(poldegree(col, Sv) > 1, return("column not affine in S"));
  foreach([Sa, Sb], S0,
    my(P0 = subst(col, Sv, S0), vv = (rr + 1/2)^2/8);
    my(Q = subst(subst(P0, W, th*rr^2), v, vv)*vv^2);
    Q = simplify(Q);
    if(type(Q) == "t_RFRAC", return("Psi*v^2 not polynomial"));
    my(res = quadpos(Q, Rb));
    if(res != 1, return(Str("positivity failed at S = ", S0*1., ": ", res))));
  1};
checkfile(fname, lin, Rneed) = {
  my(F = mforms(lin), L = readvec(fname), bad = 0, iv = List());
  for(i = 1, #L, my(c = L[i]); if(c[1] != lin, error("lin mismatch"));
    my(r = checkcert(c, F)); if(r != 1, bad++; print("certificate ", i, " R = [", c[2], ", ", c[3], "] FAILED: ", r), listput(iv, [c[2], c[3]])));
  \\ coverage of [0, Rneed]
  my(V = vecsort(Vec(iv), 1), reach = 0);
  for(i = 1, #V, if(V[i][1] > reach, break); reach = max(reach, V[i][2]));
  print(fname, ": ", #L, " certificates, ", bad, " failed; verified intervals cover R in [0, ", reach, "] (need ", Rneed, "): ",
        if(reach >= Rneed && bad == 0, "COVERED", "NOT covered"));
  print("   S range covered: (", (9-lin)/32, ", ", SofR(reach, 9 - lin)*1., "]")};
