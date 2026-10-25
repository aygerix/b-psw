\\ src/cs/e3/scan_core.gp -- E3 scanning core shared by rho_scan.gp and p23_scan.gp.
\\ scanpolys(f, lab, d, Xmax, Mmax, fh): f = [f1,f2,f3] in Q[x]; scans x = X/d, |X| <= Xmax, over the residues mod M
\\ (M = lcm of the content denominators, skipped if > Mmax) where all three values are integers; keeps values > 3,
\\ distinct, all (pseudo)prime.  Per tuple: Selfridge D, Jacobi symbols, sprp2, slprp, Korselt (p_j - 1 | n - 1),
\\ Fermat defect dF = lcm odd(ord_p 2)/gcd(., n-1), Lucas defects dL_j = odd(z_j)/gcd(odd(z_j), n+1), dLt (lcm version).
read("../../fam/bpsw.gp");
oddp(m) = m >> valuation(m, 2);
zrank(p, Q, D) = {
  my(N = p - kronecker(D, p), fa = factor(N), z = N);
  for (i = 1, #fa~, my(l = fa[i,1]);
    for (e = 1, fa[i,2], if (lucUV(z / l, Q, D, p)[1] == 0, z /= l, break)));
  z };
STATS = vector(12); STATS[10] = 1.0;
\\ STATS = [tuples, allinert, sprp2, slprp, bpsw, korselt_all, dF==1, dL1==1, sum lg dLt/lg n, min, sum lg dL2/lg p2, sum lg dL3/lg p3]
scanpolys(f, lab, d, Xmax, Mmax, fh) = {
  my(g = vector(3, j, subst(f[j], x, x/d)), M = 1, res = List(), cnt = 0);
  for (j = 1, 3, M = lcm(M, denominator(content(g[j]))));
  if (M > Mmax, return(-1));
  for (X = 0, M - 1, my(ok = 1); for (j = 1, 3, if (denominator(subst(g[j], x, X)) != 1, ok = 0; break)); if (ok, listput(res, X)));
  if (#res == 0, return(0));
  forstep (X0 = M * floor(-Xmax / M), Xmax, M,
    for (i = 1, #res, my(X = X0 + res[i], v = vector(3, j, subst(g[j], x, X)), ok = 1);
      for (j = 1, 3, if (v[j] < 5, ok = 0; break)); if (!ok, next);
      if (#Set(v) < 3, next);
      for (j = 1, 3, if (!ispseudoprime(v[j]), ok = 0; break)); if (!ok, next);
      cnt++;
      my(n = v[1]*v[2]*v[3], D = selD(n), Q = (1 - D)/4, kr = vector(3, j, kronecker(D, v[j])),
         s2 = sprp2(n), sl = slprp(n), kor = 1, AF = 1, AL = 1, dL = vector(3));
      for (j = 1, 3, if ((n - 1) % (v[j] - 1), kor = 0));
      for (j = 1, 3, AF = lcm(AF, oddp(znorder(Mod(2, v[j])))));
      for (j = 1, 3, if (gcd(v[j], 2*Q*D) == 1, my(z = oddp(zrank(v[j], Q, D))); AL = lcm(AL, z); dL[j] = z / gcd(z, n + 1), dL[j] = -1));
      my(dF = AF / gcd(AF, n - 1), dLt = AL / gcd(AL, n + 1), lr = if (dLt > 1, log(dLt) / log(n), 0.));
      STATS[1]++; STATS[2] += (kr == [-1,-1,-1]); STATS[3] += s2; STATS[4] += sl; STATS[5] += (s2 && sl);
      STATS[6] += kor; STATS[7] += (dF == 1); STATS[8] += (dL[1] == 1); STATS[9] += lr; STATS[10] = min(STATS[10], lr);
      if (dL[2] > 0, STATS[11] += log(max(dL[2],1)) / log(v[2])); if (dL[3] > 0, STATS[12] += log(max(dL[3],1)) / log(v[3]));
      filewrite(fh, Str(n, " ", v[1], " ", v[2], " ", v[3], " ", lab, " ", d, " ", X, " ", D, " ", s2, " ", sl, " ", kor, " ", dF, " ", dL[1], " ", dL[2], " ", dL[3], " ", dLt));
      if (s2 && sl, print("BPSW HIT ", lab, " d=", d, " X=", X, " n=", n, " factors=", v))));
  cnt };
printstats() = {
  print("prime tuples: ", STATS[1], "  all inert: ", STATS[2], "  sprp2: ", STATS[3], "  slprp: ", STATS[4], "  BPSW: ", STATS[5]);
  print("Korselt p_j-1 | n-1 for all j: ", STATS[6], "  Fermat defect 1: ", STATS[7], "  member-1 Lucas defect 1: ", STATS[8]);
  if (STATS[1], printf("Lucas defect log dL/log n: mean %.4f  min %.4f   mean log dL2/log p2 %.4f   mean log dL3/log p3 %.4f\n",
     STATS[9]/STATS[1], STATS[10], STATS[11]/STATS[1], STATS[12]/STATS[1])) };
