\\ Recompute the Lucas-forced primes p > 2^33 whose odd Fermat index is >= 10 (from full_L.out records).
fibmod(m, p) = polcoeff(lift(Mod(Mod(1, p) * 'x, 'x^2 - 'x - 1)^m), 1);
oddpart(n) = n >> valuation(n, 2);
showf(f) =
{
  my(c = 1, fa);
  fordiv(f, d, c *= fibonacci(d)^moebius(f / d));
  fa = factor(c);
  for (i = 1, #fa~, my(p = fa[i, 1], N, o, idx);
    if (p <= 5 || f % p == 0 || kronecker(5, p) != -1 || p < 2^33, next);
    N = p - 1; o = znorder(Mod(2, p), [N, factor(N)]); idx = oddpart(N) / oddpart(o);
    if (idx >= 10, print("f=", f, " factor(f)=", factor(f), " p=", p, " bits=", #binary(p),
        " p mod 8=", p % 8, " oddidx=", idx, " =", factor(idx), " cofactor(Phi*/p) bits=", #binary(c / p^fa[i,2]),
        " #prime factors of Phi*=", #fa~)));
}
