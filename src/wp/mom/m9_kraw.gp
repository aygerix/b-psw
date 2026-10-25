\\ m9 line: RIGOROUS existence certificate (Krawczyk test in exact rational interval arithmetic) for a real point
\\ satisfying the four boundary equations eqs(m, lin) (brute force, cubic_moments.gp) inside the open admissible domain
\\   u_i^2 < 8 v_i,  (1 - u_i)^2 < 8 v_i  (<=> (|u_i - 1/2| + 1/2)^2 < 8 v_i, which implies v_i > 1/32),
\\ with pairwise distinct members.
\\ Method: 2m - 4 coordinates are fixed at RATIONAL values (rounded from a numerical solution); the remaining 4 unknowns
\\ X are chosen to maximise |det| of the 4 x 4 Jacobian minor.  x~ = rational Newton approximation, Y ~ J(x~)^(-1)
\\ rational, box B = x~ + [-rho, rho]^4.  Krawczyk: K = x~ - Y F(x~) + (I - Y J(B)) (B - x~) inside int(B) =>
\\ F has a unique zero in B (and J is nonsingular on B since ||I - Y J(B)||_inf < 1).  All interval operations are
\\ exact rational arithmetic, so the certificate is a proof (no floating point enters the verdict).
read("cubic_moments.gp");
\\ ---- exact rational interval arithmetic
iadd(a, b) = [a[1] + b[1], a[2] + b[2]];
imul(a, b) = {my(p = [a[1]*b[1], a[1]*b[2], a[2]*b[1], a[2]*b[2]]); [vecmin(p), vecmax(p)]};
\\ interval Horner evaluation of a polynomial whose variables are among VL (vector of variables), IV their intervals
ieval(P, VL, IV) = {
  if(type(P) != "t_POL", return([P, P]));
  if(P == 0, return([0, 0]));
  my(z = variable(P), idx = 0, n = poldegree(P), acc);
  for(j = 1, #VL, if(VL[j] == z, idx = j));
  if(idx == 0, error("ieval: unknown variable ", z));
  acc = ieval(polcoef(P, n), VL, IV);
  forstep(j = n - 1, 0, -1, acc = iadd(imul(acc, IV[idx]), ieval(polcoef(P, j), VL, IV)));
  acc};
rat(x, D) = round(x * D) / D;
\\ krawczyk(m, lin, U0, V0, D): U0, V0 numerical solution (real vectors), D = denominator for the fixed coordinates
krawczyk(m, lin, U0, V0, D) = {
  my(k = m + lin, E = eqs(m, lin), U = vector(m, i, eval(Str("u", i))), W = vector(m, i, eval(Str("v", i))),
     ALL = concat(U, W), X0 = concat(U0, V0), n = 2*m, Jn, best = 0, bsel, sel, fix, F, J, xt, Y, rho, B, JB, M, K, ok);
  \\ numerical Jacobian of the 4 equations w.r.t. all 2m coordinates
  Jn = matrix(4, n, r, c, substvec(deriv(E[r], ALL[c]), ALL, X0));
  forsubset([n, 4], s, my(dd = abs(matdet(vecextract(Jn, "1..4", Vec(s))))); if(dd > best, best = dd; bsel = Vec(s)));
  sel = bsel; print("  unknown coordinates: ", vector(4, j, ALL[sel[j]]), "  |det| of numerical minor = ", best);
  fix = setminus([1..n], Set(sel));
  \\ substitute rational fixed coordinates
  my(Fv = vector(#fix, j, ALL[fix[j]]), Fq = vector(#fix, j, rat(X0[fix[j]], D)), XV = vector(4, j, ALL[sel[j]]));
  F = vector(4, r, substvec(E[r], Fv, Fq));
  J = matrix(4, 4, r, c, deriv(F[r], XV[c]));
  \\ Newton in high precision
  xt = vector(4, j, X0[sel[j]]);
  for(it = 1, 40, my(Fx = vector(4, r, substvec(F[r], XV, xt))~, Jx = matrix(4, 4, r, c, substvec(J[r, c], XV, xt)));
    xt = xt - (matsolve(Jx, Fx))~);
  my(Fx = vector(4, r, substvec(F[r], XV, xt))); print("  Newton residual (float): ", normlp(Fx, oo));
  print("  Newton correction of the unknowns from the numerical start: ", normlp(xt - vector(4, j, X0[sel[j]]), oo));
  xt = vector(4, j, rat(xt[j], 10^110));
  Y = matrix(4, 4, r, c, substvec(J[r, c], XV, 1.0 * xt)); Y = Y^(-1); Y = matrix(4, 4, r, c, rat(Y[r, c], 10^40));
  rho = 1/10^60;
  B = vector(4, j, [xt[j] - rho, xt[j] + rho]);
  my(Fxt = vector(4, r, substvec(F[r], XV, xt))~);     \\ exact rational
  my(YF = Y * Fxt);
  JB = matrix(4, 4, r, c, ieval(J[r, c], XV, B));
  \\ M = I - Y J(B) as interval matrix; bound row sums of |M|
  my(rowb = vector(4), okk = 1);
  for(r = 1, 4, my(sb = 0);
    for(c = 1, 4, my(acc = [(r == c), (r == c)]);
      for(l = 1, 4, my(y = Y[r, l]); acc = iadd(acc, imul([-y, -y], JB[l, c])));
      sb += max(abs(acc[1]), abs(acc[2])));
    rowb[r] = sb;
    if(abs(YF[r]) + rho * sb >= rho, okk = 0));
  print("  Krawczyk: max |Y F(x~)| = ", 1.0 * vecmax(abs(YF~)), ", rho = ", 1.0 * rho, ", row sums of |I - Y J(B)|: ", 1.0 * rowb);
  print("  KRAWCZYK ", if(okk, "PASSED: unique zero of the 4 x 4 system in the box", "FAILED"));
  \\ domain + distinctness checks over the box (exact rational intervals)
  my(IU = vector(m), IVv = vector(m), dom = 1, dist = 1, minmarg = oo);
  for(i = 1, m, my(iu = 0, iv = 0);
    for(j = 1, 4, if(sel[j] == i, iu = j); if(sel[j] == m + i, iv = j));
    IU[i] = if(iu, B[iu], my(q = rat(U0[i], D)); [q, q]);
    IVv[i] = if(iv, B[iv], my(q = rat(V0[i], D)); [q, q]);
    my(u2 = max(IU[i][1]^2, IU[i][2]^2), w2 = max((1 - IU[i][1])^2, (1 - IU[i][2])^2), v8 = 8 * IVv[i][1]);
    if(!(u2 < v8 && w2 < v8 && IVv[i][1] > 0), dom = 0);
    minmarg = min(minmarg, v8 - max(u2, w2)));
  for(i = 1, m, for(j = i + 1, m,
    if(!(IU[i][2] < IU[j][1] || IU[j][2] < IU[i][1] || IVv[i][2] < IVv[j][1] || IVv[j][2] < IVv[i][1]), dist = 0)));
  print("  domain inequalities u^2 < 8v, (1-u)^2 < 8v on the whole box: ", if(dom, "HOLD", "FAIL"), "  (min of 8v - max(u^2,(1-u)^2): ", 1.0 * minmarg, ")");
  print("  members pairwise distinct on the box: ", if(dist, "YES", "NO"));
  [okk && dom && dist, sel, xt, Fq, fix]};
