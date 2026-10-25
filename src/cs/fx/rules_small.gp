\\ FX (1): exact 2-adic class rules for the Selfridge Lucas sequences U(P = 1, Q = (1 - D)/4), checked on all primes.
\\ For p not dividing 2 D Q: eps = (D/p), m = p - eps, v = v2(m); z = rank of apparition (z | m); j2 = v2(z);
\\ j1 = v2(ord_p 2), w = v2(p - 1).
\\ Rules tested (each a PROVED statement, the check is a guard against slips):
\\  R0  (Fermat)   j1 = w  iff (2/p) = -1.
\\  R1  (any D)    j2 = v  iff (Q/p) = -1   [theta = alpha/beta = alpha^2/Q is a non-square in the cyclic group G_p of order m
\\                 iff Q is a non-residue; inert case via x -> x^(1-p): F_{p^2}^* -> T, kernel F_p^* inside the squares].
\\  R2  (Q = c^2)  theta = x^2 with x = alpha/c = (alpha + c)/(beta + c), N(alpha + c) = c(2c + 1); x is a non-square in G_p iff
\\                 (c(2c+1)/p) = -1, and then j2 = v - 1; if (c(2c+1)/p) = +1 then j2 <= max(v - 2, 0).  D = -15: c = 2, 10.
\\  R3  (Q = -c^2, p == 3 mod 4)  j2 = v (special case of R1);  (p == 1 mod 4, inert) v = 1 so j2 = 0.
\\ Output: violations (must be 0) and the class table (j1, j2) for inert / split primes, per D.
XMAX = 100000;
DL = [5, -7, -11, 13, -15, 17, -19, -23, 29, -31];
U(k, Q, p) = (Mod([1, -Q; 1, 0], p)^k)[2, 1];
j2of(p, D) =
{
  my(Q = (1 - D) / 4, e = kronecker(D, p), m = p - e, v = valuation(m, 2), t = m >> v, j = 0);
  if (U(m, Q, p) != 0, error("z does not divide m at p = ", p, " D = ", D));
  while (U(t << j, Q, p) != 0, j++);
  [j, v, e];
}
j1of(p) = valuation(znorder(Mod(2, p)), 2);
check(X) =
{
  my(viol = vector(4), nch = 0, tab = Map());
  forprime(p = 7, X,
    my(j1 = j1of(p), w = valuation(p - 1, 2));
    if ((j1 == w) != (kronecker(2, p) == -1), viol[1]++);
    for (k = 1, #DL,
      my(D = DL[k], Q = (1 - D) / 4, r, j2, v, e);
      if (D % p == 0 || Q % p == 0, next);
      r = j2of(p, D); j2 = r[1]; v = r[2]; e = r[3]; nch++;
      if ((j2 == v) != (kronecker(Q, p) == -1), viol[2]++);
      if (D == -15,
        if (kronecker(10, p) == -1, if (j2 != v - 1, viol[3]++), if (j2 > max(v - 2, 0), viol[3]++)));
      if (D == 5 || D == 17,
        if (p % 4 == 3 && j2 != v, viol[4]++);
        if (p % 4 == 1 && e == -1 && j2 != 0, viol[4]++));
      my(key = [D, e, j1, j2], c = 0);
      if (mapisdefined(tab, key, &c), mapput(tab, key, c + 1), mapput(tab, key, 1))));
  print("primes 7..", X, ": (p, D) checks ", nch, "; violations R0 R1 R2(-15) R3(5,17) = ", viol);
  tab;
}
report(tab) =
{
  my(M = Mat(tab));
  for (s = 1, 2, my(D = [5, -15][s]);
    foreach([-1, 1], e,
      print("D = ", D, if (e == -1, " inert", " split"), ":");
      my(tot = 0); for (i = 1, #M~, if (M[i, 1][1] == D && M[i, 1][2] == e, tot += M[i, 2]));
      for (i = 1, #M~, if (M[i, 1][1] == D && M[i, 1][2] == e && M[i, 2] >= tot / 400,
        print("   (", M[i, 1][3], ",", M[i, 1][4], ")  ", M[i, 2], "  ", Strprintf("%.4f", M[i, 2] / tot))))));
  print("inert primes with j1 = 2, by j2 (all ten D):");
  for (k = 1, #DL, my(D = DL[k], c = vector(4));
    for (i = 1, #M~, my(K = M[i, 1]); if (K[1] == D && K[2] == -1 && K[3] == 2, c[min(K[4], 3) + 1] += M[i, 2]));
    print("   D = ", D, "  Q = ", (1 - D) / 4, "  j2 = 0,1,2,>=3: ", c, "  frac(j2=0) = ", Strprintf("%.4f", c[1] / vecsum(c))));
}
report(check(XMAX));
