\\ vk_indep.gp -- INDEPENDENT verifier (line vk) for the m9 feasibility certificates.
\\ (1) Boundary equations by a different route: V = polynomial part of (F-1)/(t prod g_i) = Euclidean QUOTIENT
\\     of an exact polynomial division by a MONIC polynomial in x (= t), coefficients in Z[u, v]:
\\       v f_i = n_i = v(2x+1) + x(x+1)(x+u),   v g_i = q_i = x^2 + (1+u)x + u + 2v,
\\       no linear member:  V = (prod n_i - prod v_i)            \ (x prod q_i)
\\       linear 2t+1:       V = ((2x+1) prod n_i - prod v_i)/2   \ (x prod q_i)
\\     E = [V(-1) - 2^(1-k), V'(-1) - (1-k)2^(1-k)] for the family and for its reflection u -> 1 - u.
\\ (2) Exact comparison with eqs(m, lin) of cubic_moments.gp at random rational points.
\\ (3) Own Krawczyk existence test on an own slice (own denominators, own rho), with the Jacobian over the box
\\     enclosed by an exact Taylor-shift bound (sum of |coefficients| of J(xc + d) - J(xc) at |d_i| = rho),
\\     not by interval Horner.  Everything that enters the verdict is an exact rational.
\p 200
bnd(UU, VV, lin) = {
  my(m = #UU, N = 1, D = x, pv = 1, k = m + lin, Q);
  for(i = 1, m, N *= VV[i]*(2*x+1) + x*(x+1)*(x+UU[i]); D *= x^2 + (1+UU[i])*x + UU[i] + 2*VV[i]; pv *= VV[i]);
  if(lin, N = ((2*x+1)*N - pv)/2, N = N - pv);
  Q = divrem(N, D, x)[1];
  if(poldegree(Q, x) != k-1 || polcoef(Q, k-1, x) != 1, error("V not monic of degree k-1"));
  [subst(Q, x, -1) - 2^(1-k), subst(deriv(Q, x), x, -1) - (1-k)*2^(1-k)]};
allbnd(UU, VV, lin) = concat(bnd(UU, VV, lin), bnd(vector(#UU, i, 1 - UU[i]), VV, lin));
\\ absolute-coefficient polynomial (recursive over all variables)
absP(P) = if(type(P) != "t_POL", abs(P), my(z = variable(P)); sum(j = 0, poldegree(P), absP(polcoef(P, j)) * z^j));
rat(z, Dn) = round(z * Dn) / Dn;
\\ numerical Jacobian of allbnd w.r.t. all 2m coordinates (central differences, h = 1e-60)
numjac(U0, V0, lin) = {
  my(m = #U0, h = 1e-60, X0 = concat(U0, V0), J = matrix(4, 2*m));
  for(c = 1, 2*m, my(Xp = X0, Xm = X0); Xp[c] += h; Xm[c] -= h;
    my(Fp = allbnd(Xp[1..m], Xp[m+1..2*m], lin), Fm = allbnd(Xm[1..m], Xm[m+1..2*m], lin));
    for(r = 1, 4, J[r, c] = (Fp[r] - Fm[r]) / (2*h)));
  J};
\\ independent Krawczyk on the slice 'sel' (4 indices into concat(U, V)), fixed coords rounded with denominator Dn
vkraw(U0, V0, lin, sel, Dn, rho) = {
  my(m = #U0, n = 2*m, X0 = concat(U0, V0), yv = [y1, y2, y3, y4], dv = [d1, d2, d3, d4], Xs = vector(n), F, J, xc,
     Y, Fxc, YF, C, R, ok = 1, rows = vector(4));
  for(c = 1, n, Xs[c] = rat(X0[c], Dn));
  for(j = 1, 4, Xs[sel[j]] = yv[j]);
  F = allbnd(Xs[1..m], Xs[m+1..n], lin);
  J = matrix(4, 4, r, c, deriv(F[r], yv[c]));
  xc = vector(4, j, X0[sel[j]]);
  for(it = 1, 30, my(Fx = vector(4, r, substvec(F[r], yv, xc))~, Jx = matrix(4, 4, r, c, substvec(J[r, c], yv, xc)));
    xc = xc - (matsolve(Jx, Fx))~);
  print("   Newton float residual: ", normlp(vector(4, r, substvec(F[r], yv, xc)), oo), "   moved from start: ",
        normlp(xc - vector(4, j, X0[sel[j]]), oo));
  xc = vector(4, j, rat(xc[j], 10^130));
  Y = matrix(4, 4, r, c, substvec(J[r, c], yv, 1.0 * xc))^(-1); Y = matrix(4, 4, r, c, rat(Y[r, c], 10^50));
  Fxc = vector(4, r, substvec(F[r], yv, xc))~;  YF = Y * Fxc;
  my(Jc = matrix(4, 4, r, c, substvec(J[r, c], yv, xc)), dJ = matrix(4, 4));
  for(r = 1, 4, for(c = 1, 4,
    my(Ps = substvec(J[r, c], yv, vector(4, j, xc[j] + dv[j])), A = absP(Ps));
    dJ[r, c] = substvec(A, dv, vector(4, j, rho)) - abs(substvec(Ps, dv, [0, 0, 0, 0]))));
  C = matid(4) - Y * Jc;
  for(r = 1, 4, rows[r] = sum(c = 1, 4, abs(C[r, c]) + sum(l = 1, 4, abs(Y[r, l]) * dJ[l, c]));
    if(!(abs(YF[r]) + rho * rows[r] < rho), ok = 0));
  print("   max |Y F(xc)| = ", 1.0 * vecmax(abs(YF~)), "   max dJ = ", 1.0 * vecmax(abs(concat(Vec(dJ)))),
        "   row sums |I - Y J(B)| <= ", 1.0 * rows);
  \\ domain + distinctness over the box
  my(IU = vector(m), IW = vector(m), dom = 1, dist = 1, marg = oo);
  for(i = 1, m,
    my(ju = 0, jw = 0); for(j = 1, 4, if(sel[j] == i, ju = j); if(sel[j] == m + i, jw = j));
    IU[i] = if(ju, [xc[ju] - rho, xc[ju] + rho], [Xs[i], Xs[i]]);
    IW[i] = if(jw, [xc[jw] - rho, xc[jw] + rho], [Xs[m+i], Xs[m+i]]);
    my(lo = IU[i][1], hi = IU[i][2], top = max(max(lo^2, hi^2), max((1-lo)^2, (1-hi)^2)));
    if(!(IW[i][1] > 0 && top < 8 * IW[i][1]), dom = 0);
    marg = min(marg, 8 * IW[i][1] - top));
  for(i = 1, m, for(j = i+1, m,
    if(!(IU[i][2] < IU[j][1] || IU[j][2] < IU[i][1] || IW[i][2] < IW[j][1] || IW[j][2] < IW[i][1]), dist = 0)));
  print("   KRAWCZYK ", if(ok, "PASS", "FAIL"), ";  domain ", if(dom, "HOLDS", "FAILS"), " (margin ", 1.0*marg,
        ");  distinct ", if(dist, "YES", "NO"));
  \\ moment sanity at the centre (float, informational only)
  my(uu = vector(m, i, (IU[i][1] + IU[i][2])/2), ww = vector(m, i, (IW[i][1] + IW[i][2])/2), S = vecsum(ww),
     K = sum(i = 1, m, ww[i]*(uu[i] - 1/2)^2), M = sum(i = 1, m, ww[i]^2), Mmax = (S - (m-1)/32)^2 + (m-1)/1024);
  print("   centre: S = ", 1.0*S, "  K = ", 1.0*K, "  Mmax - M = ", 1.0*(Mmax - M), "  min v = ", 1.0*vecmin(ww),
        "  exact F at centre (float of exact rational): ", 1.0 * vecmax(abs(Fxc~)));
  [ok && dom && dist, xc, Xs]};
