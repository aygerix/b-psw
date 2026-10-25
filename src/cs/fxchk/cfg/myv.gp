\\ Adversarial FX check: own verification of a primal configuration.  For each listed prime: proven primality, class
\\ (j1, j2) with j1 = v2(ord_p 2), j2 = v2(z_D(p)); z from x^k mod (x^2 - x + Q) over F_p (coefficient of x = U_k).
\\ Usable iff odd(ord) | E1 and odd(z) | E2.  Also recomputes the MINIMAL caps needed by the usable primes and its bits.
Ux(k, Q, p) = polcoeff(lift(lift(Mod(Mod(1, p) * x, Mod(1, p) * (x^2 - x + Q))^k)), 1);
zr(p, Q, m) = { my(z = m, f = factor(m)[, 1]); if (Ux(m, Q, p) != 0, error("U_m != 0 at ", p)); for (i = 1, #f, while (z % f[i] == 0 && Ux(z / f[i], Q, p) == 0, z /= f[i])); z; }
bitsof(caps) = sum(i = 1, #caps, log(eulerphi(caps[i][1]^caps[i][2])) / log(2));
myv(name, D, j1, j2, E1, E2, L) =
{
  my(Q = (1 - D) / 4, e1 = prod(i = 1, #E1, E1[i][1]^E1[i][2]), e2 = prod(i = 1, #E2, E2[i][1]^E2[i][2]), ok = 0, bad = 0,
     c1 = Map(), c2 = Map(), maxq = 0, ninert = 0, dup = #L - #Set(L));
  if (gcd(e1, e2) != 1, error("listed rows on both sides"));
  for (i = 1, #L, my(p = L[i], js, o, z, o1, o2, f, c);
    if (!isprime(p), bad++; next);
    js = kronecker(D, p); if (js == 0 || Q % p == 0, bad++; next);
    o = znorder(Mod(2, p)); z = zr(p, Q, p - js);
    if (valuation(o, 2) != j1 || valuation(z, 2) != j2, bad++; next);
    o1 = o >> j1; o2 = z >> j2;
    if (e1 % o1 || e2 % o2, bad++; next);
    ok++; ninert += (js == -1);
    f = factor(o1); for (k = 1, #f~, maxq = max(maxq, f[k, 1]); if (!mapisdefined(c1, f[k, 1], &c) || c < f[k, 2], mapput(c1, f[k, 1], f[k, 2])));
    f = factor(o2); for (k = 1, #f~, maxq = max(maxq, f[k, 1]); if (!mapisdefined(c2, f[k, 1], &c) || c < f[k, 2], mapput(c2, f[k, 1], f[k, 2]))));
  my(m1 = Mat(c1), m2 = Mat(c2), cm1 = if (#m1, vector(#m1~, i, [m1[i, 1], m1[i, 2]]), []), cm2 = if (#m2, vector(#m2~, i, [m2[i, 1], m2[i, 2]]), []));
  my(both = setintersect(Set(apply(v -> v[1], cm1)), Set(apply(v -> v[1], cm2))));
  my(bl = log(eulerphi(e1 * e2)) / log(2), bm = bitsof(cm1) + bitsof(cm2));
  printf("%s: D=%d (%d,%d) listed %d (dups %d) usable %d (inert %d) failed %d | listed-cap bits %.4f | minimal-cap bits %.4f, rows on both sides %d, max row %d | N-0.2b %.4f N-0.225b %.4f N-0.25b %.4f (minimal caps)\n",
         name, D, j1, j2, #L, dup, ok, ninert, bad, bl, bm, #both, maxq, ok - 0.2 * bm, ok - 0.225 * bm, ok - 0.25 * bm);
}
