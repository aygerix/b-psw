\\ List forced primes p > 2^33 that are usable at rows <= T (both orders' odd parts T-smooth), with full profile.
fibmod(m, p) = polcoeff(lift(Mod(Mod(1, p) * 'x, 'x^2 - 'x - 1)^m), 1);
oddpart(n) = n >> valuation(n, 2);
zent(p, N, fa) = { my(n = N); for (i = 1, #fa~, my(q = fa[i, 1]); for (j = 1, fa[i, 2], if (fibmod(n / q, p) == 0, n /= q, break))); n; }
gpf(n) = my(f = factor(oddpart(n))); if (#f~ == 0, 1, f[#f~, 1]);
prof(p, tag, ef) =
{
  my(o = znorder(Mod(2, p), [p - 1, factor(p - 1)]), z = zent(p, p + 1, factor(p + 1)));
  print(tag, " ", ef, " p=", p, " (", #binary(p), " bits, p mod 8 = ", p % 8, ") ord2=", factor(o), " z=", factor(z),
        " class (v2 ord, v2 z) = (", valuation(o, 2), ",", valuation(z, 2), ") gpf odd ord ", gpf(o), " gpf odd z ", gpf(z));
}
showL(f, T) =
{
  my(c = 1, fa); fordiv(f, d, c *= fibonacci(d)^moebius(f / d)); fa = alarm(60, factor(c)); if (type(fa) == "t_ERROR", return);
  for (i = 1, #fa~, my(p = fa[i, 1]); if (p < 2^33 || f % p == 0 || kronecker(5, p) != -1, next);
    my(o = znorder(Mod(2, p), [p - 1, factor(p - 1)])); if (gpf(o) <= T && gpf(f) <= T, prof(p, "LF", f)));
}
showF(e, T) =
{
  my(fa = alarm(60, factor(polcyclo(e, 2)))); if (type(fa) == "t_ERROR", return);
  for (i = 1, #fa~, my(p = fa[i, 1]); if (p < 2^33 || e % p == 0 || kronecker(5, p) != -1, next);
    my(z = zent(p, p + 1, factor(p + 1))); if (gpf(z) <= T && gpf(e) <= T, prof(p, "FF", e)));
}
