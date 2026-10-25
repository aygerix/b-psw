\\ Check the four 2-adic rules for D = 5 and the ord rule over all primes 7 <= p <= LIM (independent of census).
\\ usage: echo 'LIM=10^6; read("aud_2adic.gp")' | gp -q
zf(p) = my(kr = kronecker(5, p), a); if(kr == 1, my(r = polrootsmod(x^2 - x - 1, p)); znorder(Mod(r[1], p)/Mod(r[2], p)), a = ffgen(Mod(1,p)*(x^2 - x - 1), 'a); fforder(a/(1 - a)));
bad = vector(5); cnt = vector(5);
{
forprime(p = 7, LIM,
  my(kr = kronecker(5, p), z = zf(p), o = znorder(Mod(2, p)));
  if(kr == -1 && p % 4 == 1, cnt[1]++; if(z % 2 == 0, bad[1]++));
  if(kr == -1 && p % 4 == 3, cnt[2]++; if(valuation(z, 2) != valuation(p + 1, 2), bad[2]++));
  if(kr == 1 && p % 4 == 3, cnt[3]++; if(valuation(z, 2) != 1, bad[3]++));
  if(kr == 1 && p % 4 == 1, cnt[4]++; if(((p - 1)/2) % z != 0, bad[4]++));
  cnt[5]++; if(kronecker(2, p) == 1, if(((p - 1)/2) % o != 0, bad[5]++), if(valuation(o, 2) != valuation(p - 1, 2), bad[5]++)));
}
print("rules [inert 1mod4 z odd, inert 3mod4 v2 z = v2(p+1), split 3mod4 v2 z = 1, split 1mod4 z | (p-1)/2, ord rule]");
print("tested ", cnt, " violations ", bad);
