\\ m9vmf_indep.gp -- ADVERSARIAL, INDEPENDENT check of the m9 moment form printed in m9_mform.out (k = 9, 11).
\\ Shares no code with cubic_moments.gp / m9_moments.gp / m9_lib.gp (no xser, trunc, ellv, Tform, evalmom, reflect).
\\ E-side computed directly: V = Euclidean quotient of F - 1 by t*G, F = prod f_i, G = prod (f_i - 1)/t,
\\ f_i = 2t + 1 + t(t+1)(a_i t + c_i), linear member 2t + 1 (g = 2); reflection computed as -f(-1-t) directly.
\\ (A) exact rational random points, ALL m members, reflected family read off from the reflected polynomials;
\\ (B) exact SYMBOLIC identity with r active members (u_i, v_i symbolic), V via monic division:
\\     t^(k-1) prod(1 + x_i) (1 + s/2)^lin = t^(k-1-r-lin) (t + 1/2)^lin prod(v_i + t gh_i) / prod gh_i,
\\     gh_i = (f_i - 1)/(t a_i) monic.  Combined with the structural count (each monomial of E involves at most
\\     floor((k-1)/3) members, each monomial of the moment form at most deg_M members) this proves the identity.
tvar = t;
pre(s, p) = #s >= #p && strjoin(Vec(s)[1..#p]) == p;
post(s, n) = strjoin(Vec(s)[n+1..#s]);
parseout(fn) = {
  my(L = readstr(fn), cs = List(), cur = 0, pr = ["Z0+ = ", "Z1+ = ", "Z0- = ", "Z1- = "]);
  for(i = 1, #L, my(s = L[i]);
    if(pre(s, "===="),
      my(tok = strsplit(s, " "));
      if(cur, listput(cs, cur));
      cur = [eval(strjoin(Vec(tok[4])[1..#tok[4]-1])), eval(strjoin(Vec(tok[9])[1..#tok[9]-1])), eval(tok[12]), 0, 0, 0, 0];
      next);
    for(j = 1, 4, if(pre(s, pr[j]), cur[3 + j] = eval(post(s, #pr[j])))));
  if(cur, listput(cs, cur)); Vec(cs)};
\\ moment symbol M{a}_{b} -> [a, b] (or 0 if not a moment symbol)
mab(z) = {my(nm = Str(z), c = Vec(nm)); if(c[1] != "M", return(0)); my(sp = strsplit(nm, "_"));
  [eval(strjoin(Vec(sp[1])[2..#sp[1]])), eval(sp[2])]};
momval(P, U, W) = {my(R = P); foreach(variables(P), z, my(ab = mab(z)); if(ab,
  R = subst(R, z, sum(j = 1, #U, W[j]^ab[1] * (U[j] - 1/2)^ab[2])))); R};
flipodd(P) = {my(R = P); foreach(variables(P), z, my(ab = mab(z)); if(ab && ab[2] % 2, R = subst(R, z, -z))); R};
mdeg(P) = {my(R = P); foreach(variables(P), z, if(mab(z), R = subst(R, z, z*eps))); poldegree(R, eps)};
\\ direct boundary quantities of a list of member polynomials fs (k = #fs, odd)
Edir(fs) = {
  my(k = #fs, c0 = 2^(1-k), c1 = (1-k)*2^(1-k), F = prod(i = 1, k, fs[i]), G = 1);
  for(i = 1, k, my(g = (fs[i] - 1)/t); if(type(g) == "t_RFRAC", error("t does not divide f - 1")); G *= g);
  my(qr = divrem(F - 1, t*G, t), V = qr[1]);
  if(F - 1 != V*t*G + qr[2] || poldegree(qr[2], t) >= poldegree(t*G, t), error("bad division"));
  if(poldegree(V, t) != k - 1 || polcoef(V, k - 1, t) != 1, error("V not monic of degree k-1"));
  \\ sanity: the forced values hold for the rational function (F-1)/(tG) identically
  my(Rf = (F - 1)/(t*G));
  if(subst(Rf, t, -1) != c0 || subst(deriv(Rf, t), t, -1) != c1, error("forced values fail for (F-1)/(tG)"));
  [subst(V, t, -1) - c0, subst(deriv(V, t), t, -1) - c1]};
\\ read (a, c) back from a normal-form cubic a t^3 + (a + c) t^2 + (c + 2) t + 1
readac(f) = {my(a = polcoef(f, 3, t), c = polcoef(f, 2, t) - a);
  if(poldegree(f, t) != 3 || polcoef(f, 0, t) != 1 || polcoef(f, 1, t) != c + 2, error("not normal form")); [a, c]};
rr() = (random(4001) - 2000) / (random(53) + 1);
partA(cs, npts) = {
  my(m = cs[2], lin = cs[3], k = cs[1], Z0 = cs[4] + cs[6], Z1 = cs[5] + cs[7], bad = 0);
  for(it = 1, npts,
    my(A = vector(m, i, my(a = 0); while(a == 0, a = rr()); a), C = vector(m, i, rr()));
    my(fs = vector(m, i, 2*t + 1 + t*(t+1)*(A[i]*t + C[i])));
    if(lin, fs = concat(fs, [2*t + 1]));
    my(fr = vector(#fs, i, -subst(fs[i], t, -1 - t)));
    if(lin && fr[#fr] != 2*t + 1, error("linear member not fixed by reflection"));
    my(E = Edir(fs), Er = Edir(fr));
    my(U = vector(m, i, C[i]/A[i]), W = vector(m, i, 1/A[i]));
    my(acr = vector(m, i, readac(fr[i])), Ur = vector(m, i, acr[i][2]/acr[i][1]), Wr = vector(m, i, 1/acr[i][1]));
    if(Ur != vector(m, i, 1 - U[i]) || Wr != W, error("reflection is not u -> 1-u, v -> v"));
    my(z0 = momval(Z0, U, W), z1 = momval(Z1, U, W), z0r = momval(Z0, Ur, Wr), z1r = momval(Z1, Ur, Wr));
    my(res = [E[1] - z0, E[2] - (-(k-1)*z0 - z1), Er[1] - z0r, Er[2] - (-(k-1)*z0r - z1r)]);
    if(res != [0,0,0,0], bad++; print("   MISMATCH at point ", it, ": ", res)));
  bad};
partB(cs, r) = {
  my(m = cs[2], lin = cs[3], k = cs[1], Z0 = cs[4] + cs[6], Z1 = cs[5] + cs[7], c0 = 2^(1-k), c1 = (1-k)*2^(1-k));
  if(k - 1 - r - lin < 0, error("too many active members"));
  my(U = vector(r, i, eval(Str("p", i))), W = vector(r, i, eval(Str("q", i))), res = vector(4));
  for(rf = 0, 1,
    \\ gh_i from f_i = 1 + t gh_i / v_i; reflected gh from -f(-1-t) = 1 + t gh^r / v  (computed, not assumed)
    my(gh = vector(r, i, t^2 + (1 + U[i])*t + U[i] + 2*W[i]));
    if(rf, gh = vector(r, i, my(fi = 1 + t*gh[i]/W[i], frf = -subst(fi, t, -1 - t), g2 = (frf - 1)*W[i]/t);
      if(type(g2) == "t_RFRAC" || poldegree(g2, t) != 2 || polcoef(g2, 2, t) != 1, error("reflected gh not monic quadratic")); g2));
    my(N = t^(k-1-r-lin) * (t + 1/2)^lin * prod(i = 1, r, W[i] + t*gh[i]), D = prod(i = 1, r, gh[i]));
    my(qr = divrem(N, D, t), V = qr[1]);
    if(N != V*D + qr[2] || poldegree(qr[2], t) >= poldegree(D, t), error("bad symbolic division"));
    my(z0 = momval(Z0, U, W), z1 = momval(Z1, U, W));
    if(rf, z0 = momval(flipodd(Z0), U, W); z1 = momval(flipodd(Z1), U, W));
    res[2*rf + 1] = subst(V, t, -1) - c0 - z0;
    res[2*rf + 2] = subst(deriv(V, t), t, -1) - c1 - (-(k-1)*z0 - z1));
  res};
{
  setrand(20261025);
  my(CS = parseout("m9_mform.out"));
  print("parsed cases: ", vector(#CS, i, CS[i][1..3]));
  foreach(CS, cs, my(k = cs[1], m = cs[2], lin = cs[3]);
    print("==== k = ", k, ", m = ", m, ", lin = ", lin);
    \\ parity of the split: Z+ even, Z- odd under M(a,b) -> (-1)^b M(a,b)
    print("  parity: Z0+ even ", flipodd(cs[4]) == cs[4], ", Z1+ even ", flipodd(cs[5]) == cs[5],
          ", Z0- odd ", flipodd(cs[6]) == -cs[6], ", Z1- odd ", flipodd(cs[7]) == -cs[7]);
    my(syms = select(z -> mab(z) != 0, variables(cs[4] + cs[5] + cs[6] + cs[7])));
    print("  #moment symbols: ", #syms, "  max (a,b): a <= ", vecmax(vector(#syms, i, mab(syms[i])[1])),
          ", b <= ", vecmax(vector(#syms, i, mab(syms[i])[2])), ", max a+b = ", vecmax(vector(#syms, i, vecsum(mab(syms[i])))));
    my(dM = vecmax(vector(4, j, mdeg(cs[3 + j]))), dE = (k - 1) \ 3, R = max(dM, dE));
    print("  max product degree in moments = ", dM, ", members per E-monomial <= ", dE, " => R = ", R);
    my(t0 = getabstime(), b = partA(cs, 40));
    print("  (A) 40 random exact points, all ", m, " members, direct V and direct reflection: mismatches = ", b, "  (", getabstime() - t0, " ms)");
    for(r = R, R + 2, if(k - 1 - r - lin < 0, next); t0 = getabstime();
      print("  (B) symbolic, ", r, " active members: residuals ", partB(cs, r), "  (", getabstime() - t0, " ms)")));
}
