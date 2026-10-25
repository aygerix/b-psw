\\ E1 pool from the CHM twin-smooth list: for every n in the list with p = 2n + 1 prime (so p - 1 and p + 1 are both
\\ y-smooth), compute exact ord_p 2 and exact z_D(p) = order of gamma = alpha / alphabar, alpha = (1 + sqrt D)/2,
\\ in F_p (split) or F_{p^2} (inert), for the ten Selfridge D (P = 1, Q = (1 - D)/4).  U_k == 0 mod p iff gamma^k = 1.
\\ Output lines "D p t j1 j2 o1 o2" (t = 0 inert, 1 split; split kept only if gcd(o1, o2) = 1), same format as census3.c.
\\ Usage: gp -q chm_pool.gp  with  INFILE, OUTFILE set by the caller (default data/sn_200_aug.txt -> chm200_pool.txt).
DL = [5, -7, -11, 13, -15, 17, -19, -23, 29, -31];
odd(x) = x >> valuation(x, 2);
zrank(p, D) =
{
  my(js = kronecker(D, p), m = p - js, g, s, t);
  if (js == 1,
    s = sqrt(Mod(D, p)); g = (1 + s) / (1 - s); return(znorder(g, m)),
    t = ffgen(Mod(1, p) * (x^2 - D), 't); g = (1 + t) / (1 - t); return(fforder(g, [m, factor(m)])));
}
run(infile, outfile) =
{
  my(v = readvec(infile), fh = fileopen(outfile, "w"), np = 0, nr = vector(#DL));
  for (i = 1, #v,
    my(n = v[i], p = 2 * n + 1);
    if (p < 3 || !isprime(p), next);
    np++;
    my(o = znorder(Mod(2, p), [p - 1, factor(p - 1)]), o1 = odd(o), j1 = valuation(o, 2));
    for (k = 1, #DL,
      my(D = DL[k], Q = (1 - D) / 4, js, z, o2);
      if (D % p == 0 || Q % p == 0, next);
      js = kronecker(D, p); if (js == 0, next);
      z = zrank(p, D); o2 = odd(z);
      if (js == 1 && gcd(o1, o2) != 1, next);
      nr[k]++;
      filewrite(fh, Str(D, " ", p, " ", (js == 1) * 1, " ", j1, " ", valuation(z, 2), " ", o1, " ", o2))));
  fileclose(fh);
  print("CHM pool: ", #v, " list entries, ", np, " primes p = 2n+1; records per D: ", nr);
}
