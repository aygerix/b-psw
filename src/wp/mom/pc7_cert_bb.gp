\\ pc7_cert_bb.gp -- EXACT branch-and-bound lower bound for E_alpha = sum_i alpha_i D_i on [LO, HI]^7 at a rational shape.
\\ E = kap + sum_j a(w_j) + sum_{i<j} b(w_i, w_j) + lam3 e3 (3-body Hankel form, D-basis).  [LO, HI] is split into N intervals
\\ (geometric spacing, rational end points); rigorous lower bounds of a on each interval and of b on each pair of intervals come
\\ from Bernstein coefficients (exact rationals).  A depth-first search over sorted 7-tuples of intervals adds these bounds and
\\ prunes a subtree as soon as a rigorous lower bound for all its completions is > 0.  Returns the minimum lower bound found
\\ (> 0 proves E_alpha > 0 on [LO, HI]^7).  Inputs (globals): SHAPE = [r1, r2], ALPHA (8 rationals), LO, HI, N.
read("pc7_cert_shape.gp");
s; s2;
bern1(P, x0, x1) = {my(Q = subst(P, w, x0 + (x1 - x0)*s), d = poldegree(Q, s)); if(d < 0, return(0));
  vecmin(vector(d + 1, k, sum(i = 0, k - 1, binomial(k - 1, i) / binomial(d, i) * polcoef(Q, i, s))))};
bern2(P, x0, x1, y0, y1) = {my(Q = substvec(P, [w, v], [x0 + (x1 - x0)*s, y0 + (y1 - y0)*s2]));
  my(d1 = poldegree(Q, s), d2 = poldegree(Q, s2), co = matrix(d1 + 1, d2 + 1, i, j, polcoef(polcoef(Q, i - 1, s), j - 1, s2)), m = oo);
  for(k = 0, d1, for(l = 0, d2, my(bkl = sum(i = 0, k, sum(j = 0, l, binomial(k, i) * binomial(l, j) / (binomial(d1, i) * binomial(d2, j)) * co[i + 1, j + 1])));
    if(bkl < m, m = bkl))); m};
BBMIN = oo; BBNODES = 0;
bbrec(m, st, part, ch, LA, LB, ex, lam3) = {
  BBNODES++;
  if(m == 7, my(t = part + lam3 * sum(i = 1, 7, sum(j = i + 1, 7, sum(l = j + 1, 7, ex[ch[i]] * ex[ch[j]] * ex[ch[l]]))));
     if(t < BBMIN, BBMIN = t); return);
  my(r = 7 - m, nN = #LA);
  for(k = st, nN,
    my(np = part + LA[k] + sum(i = 1, m, LB[ch[i], k]));
    \\ rigorous bound for the r - 1 members still to come (all from intervals >= k) and their pairs
    my(rest = 0);
    if(r > 1, my(mn1 = oo, mn2 = oo);
      for(q = k, nN, my(c1 = LA[q] + sum(i = 1, m, LB[ch[i], q]) + LB[k, q]); if(c1 < mn1, mn1 = c1));
      for(q = k, nN, for(q2 = q, nN, if(LB[q, q2] < mn2, mn2 = LB[q, q2])));
      rest = (r - 1) * mn1 + binomial(r - 1, 2) * min(mn2, 0));
    \\ e3 term: >= 0 if lam3 >= 0; if lam3 < 0, >= lam3 * e3(upper end points, remaining members at HI)
    my(e3lb = if(lam3 >= 0, 0, my(z = concat(vector(m, i, ex[ch[i]]), concat([ex[k]], vector(r - 1, i, ex[nN])))); lam3 * polcoef(prod(i = 1, 7, 1 + z[i]*'T), 3, 'T)));
    if(np + rest + e3lb > 0, next);  \\ whole subtree certified
    bbrec(m + 1, k, np, concat(ch, [k]), LA, LB, ex, lam3))};
bb() = {
  my(S = shdata(SHAPE[1], SHAPE[2]), K = S[1], A = S[2], C = S[3]);
  my(kap = sum(i = 1, 8, ALPHA[i] * K[i]), lam3 = ALPHA[1] - ALPHA[5]);
  my(aP = sum(m = 1, 9, sum(i = 1, 8, ALPHA[i] * A[i, m]) * w^m));
  my(bP = sum(n = 2, 6, sum(i = 1, 8, ALPHA[i] * C[i, n - 1]) * sum(q = 1, n - 1, w^q * v^(n - q))));
  my(X = vector(N + 1, k, if(k == 1, LO, if(k == N + 1, HI, bestappr(LO * (HI/LO)^((k - 1)/N) * 1., 10^6)))));
  my(LA = vector(N, k, bern1(aP, X[k], X[k + 1])));
  my(LB = matrix(N, N)); for(k = 1, N, for(l = k, N, LB[k, l] = bern2(bP, X[k], X[k + 1], X[l], X[l + 1]); LB[l, k] = LB[k, l]));
  my(ex = vector(N, k, if(lam3 >= 0, X[k], X[k + 1])));
  BBMIN = oo; BBNODES = 0;
  bbrec(0, 1, kap, [], LA, LB, ex, lam3);
  if(BBMIN == oo, BBMIN = 1);   \\ every subtree pruned with a positive bound
  print("N = ", N, "  kap_alpha = ", kap * 1., "  lam3 = ", lam3 * 1., "  nodes ", BBNODES, "  min leaf LB = ", BBMIN * 1., "  (sign ", sign(BBMIN), ")");
  BBMIN};
