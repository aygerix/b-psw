\\ Independent audit census (PARI): for every prime 3 <= p in [a, b), exact ord_p 2 (znorder) and exact Lucas rank
\\ z_D(p) = least d | m = p - (D/p) with U_d(1, Q) == 0 mod p, U_d from 2x2 matrix powers (not the doubling ladder of
\\ census3.c, not the gamma-order of chm_pool.gp).  Keep (D, p) when odd(ord) and odd(z) are both 200-smooth, and for
\\ split p also gcd = 1.  Same record format as census3.c: "D p t j1 j2 o1 o2".
DL = [5, -7, -11, 13, -15, 17, -19, -23, 29, -31];
odd(x) = x >> valuation(x, 2);
smooth(x) = if (x == 1, 1, my(f = factor(x)); f[#f~, 1] <= 200);
Uk(k, Q, p) = lift((Mod([1, -Q; 1, 0], p)^k)[2, 1]);
zrank(p, Q, m) =
{
  my(z = m, f = factor(m)[, 1]);
  if (Uk(m, Q, p) != 0, error("U_m != 0 at p = ", p));
  for (i = 1, #f, my(l = f[i]); while (z % l == 0 && Uk(z / l, Q, p) == 0, z /= l));
  z;
}
run(a, b, outfile) =
{
  my(fh = fileopen(outfile, "w"), np = 0, nf = 0, nr = 0, t0 = getwalltime());
  forprime(p = max(a, 3), b - 1,
    np++;
    my(o = znorder(Mod(2, p)), o1 = odd(o));
    if (!smooth(o1), next);
    nf++;
    for (k = 1, #DL,
      my(D = DL[k], Q = (1 - D) / 4, js, m, z, o2);
      if (D % p == 0 || Q % p == 0, next);
      js = kronecker(D, p); if (js == 0, next);
      m = p - js; z = zrank(p, Q, m); o2 = odd(z);
      if (!smooth(o2), next);
      if (js == 1 && gcd(o1, o2) != 1, next);
      nr++;
      filewrite(fh, Str(D, " ", p, " ", (js == 1) * 1, " ", valuation(o, 2), " ", valuation(z, 2), " ", o1, " ", o2))));
  fileclose(fh);
  print("range [", a, ", ", b, "): primes ", np, ", Fermat-smooth ", nf, ", records ", nr, ", ", (getwalltime() - t0) / 1000., " s");
}
