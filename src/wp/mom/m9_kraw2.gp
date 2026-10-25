\\ m9 line: Krawczyk certificate for LEVEL 2 of the Taylor hierarchy (orders 0, 1, 2 at t = -1, both halves):
\\   E1..E4 = eqs(m, lin) (boundary equations),
\\   E5 = [e_2(s) - (2^k V''(-1) - 2(k-1)(k-2))] * prod v_i,   s_i = 2 + (1 - u_i)/v_i  (linear member: s = 2),
\\   E6 = the same for the reflected family (u -> 1 - u), i.e. with alpha_i = 2 + u_i/v_i and V~.
\\ (E5 = 0 <=> the z^2 Taylor coefficient of H - V at t = -1 vanishes:  r_2 = 2^(1-k)[(k-1)(k-2)/2 + e_2(s)/4] - V''(-1)/2.)
\\ Same exact rational interval Krawczyk test as m9_kraw.gp, with 6 unknowns.
read("m9_kraw.gp");
lev2eqs(m, lin) = {
  my(k = m + lin, U = vector(m, i, eval(Str("u", i))), W = vector(m, i, eval(Str("v", i))), E = eqs(m, lin), out = vector(6));
  for(j = 1, 4, out[j] = E[j]);
  for(refl = 0, 1, my(Uq = if(refl, vector(m, i, 1 - U[i]), U), V = Vpoly(m, lin, Uq, W), V2 = subst(deriv(deriv(V, t), t), t, -1),
       PV = prod(i = 1, m, W[i]), sn = vector(m, i, 2*W[i] + 1 - Uq[i]), e2n = 0);
    \\ e_2(s) * prod v:  pairs of cubics, plus (if lin) the linear member s_0 = 2 with each cubic
    for(i = 1, m, for(j = i + 1, m, e2n += sn[i] * sn[j] * prod(l = 1, m, if(l == i || l == j, 1, W[l]))));
    if(lin, for(i = 1, m, e2n += 2 * sn[i] * prod(l = 1, m, if(l == i, 1, W[l]))));
    out[5 + refl] = e2n - (2^k * V2 - 2*(k-1)*(k-2)) * PV);
  out};
krawczykN(E, ALL, X0, D, ne) = {
  my(n = #ALL, Jn, best = 0, bsel, sel, fix, F, J, xt, Y, rho, B, JB);
  Jn = matrix(ne, n, r, c, substvec(deriv(E[r], ALL[c]), ALL, X0));
  forsubset([n, ne], s, my(dd = abs(matdet(vecextract(Jn, Str("1..", ne), Vec(s))))); if(dd > best, best = dd; bsel = Vec(s)));
  sel = bsel; print("  unknown coordinates: ", vector(ne, j, ALL[sel[j]]), "  |det| of numerical minor = ", precision(best, 10));
  fix = setminus([1..n], Set(sel));
  my(Fv = vector(#fix, j, ALL[fix[j]]), Fq = vector(#fix, j, rat(X0[fix[j]], D)), XV = vector(ne, j, ALL[sel[j]]));
  F = vector(ne, r, substvec(E[r], Fv, Fq));
  J = matrix(ne, ne, r, c, deriv(F[r], XV[c]));
  xt = vector(ne, j, X0[sel[j]]);
  for(it = 1, 50, my(Fx = vector(ne, r, substvec(F[r], XV, xt))~, Jx = matrix(ne, ne, r, c, substvec(J[r, c], XV, xt)));
    xt = xt - (matsolve(Jx, Fx))~);
  print("  Newton residual (float): ", precision(normlp(vector(ne, r, substvec(F[r], XV, xt)), oo), 10),
        ";  correction from the numerical start: ", precision(normlp(xt - vector(ne, j, X0[sel[j]]), oo), 10));
  xt = vector(ne, j, rat(xt[j], 10^110));
  Y = matrix(ne, ne, r, c, substvec(J[r, c], XV, 1.0 * xt)); Y = Y^(-1); Y = matrix(ne, ne, r, c, rat(Y[r, c], 10^40));
  rho = 1/10^60;
  B = vector(ne, j, [xt[j] - rho, xt[j] + rho]);
  my(Fxt = vector(ne, r, substvec(F[r], XV, xt))~, YF = Y * Fxt, okk = 1, rowb = vector(ne));
  JB = matrix(ne, ne, r, c, ieval(J[r, c], XV, B));
  for(r = 1, ne, my(sb = 0);
    for(c = 1, ne, my(acc = [(r == c), (r == c)]);
      for(l = 1, ne, my(y = Y[r, l]); acc = iadd(acc, imul([-y, -y], JB[l, c])));
      sb += max(abs(acc[1]), abs(acc[2])));
    rowb[r] = sb; if(abs(YF[r]) + rho * sb >= rho, okk = 0));
  print("  Krawczyk: max |Y F(x~)| = ", precision(1.0 * vecmax(abs(YF~)), 6), ", rho = 1e-60, max row sum |I - Y J(B)| = ", precision(1.0 * vecmax(rowb), 6));
  print("  KRAWCZYK ", if(okk, "PASSED", "FAILED"));
  [okk, sel, xt, Fq, fix, B]};
\\ domain/distinctness over the box for given selection
domcheck(m, sel, B, U0, V0, D) = {
  my(IU = vector(m), IVv = vector(m), dom = 1, dist = 1, minmarg = oo);
  for(i = 1, m, my(iu = 0, iv = 0);
    for(j = 1, #sel, if(sel[j] == i, iu = j); if(sel[j] == m + i, iv = j));
    IU[i] = if(iu, B[iu], my(q = rat(U0[i], D)); [q, q]);
    IVv[i] = if(iv, B[iv], my(q = rat(V0[i], D)); [q, q]);
    my(u2 = max(IU[i][1]^2, IU[i][2]^2), w2 = max((1 - IU[i][1])^2, (1 - IU[i][2])^2), v8 = 8 * IVv[i][1]);
    if(!(u2 < v8 && w2 < v8 && IVv[i][1] > 0), dom = 0);
    minmarg = min(minmarg, v8 - max(u2, w2)));
  for(i = 1, m, for(j = i + 1, m,
    if(!(IU[i][2] < IU[j][1] || IU[j][2] < IU[i][1] || IVv[i][2] < IVv[j][1] || IVv[j][2] < IVv[i][1]), dist = 0)));
  print("  domain inequalities on the box: ", if(dom, "HOLD", "FAIL"), " (min 8v - max(u^2,(1-u)^2) = ", precision(1.0 * minmarg, 6), "); distinct: ", if(dist, "YES", "NO"));
  dom && dist};
