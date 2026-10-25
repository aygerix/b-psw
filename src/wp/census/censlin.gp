ok(F,f) = (F - f) % (f^3 - f) == 0;
crt2(a, A, b, B) =
{
  my(g = gcd(A, B), Ag, Bg, ge, u);
  if ((b - a) % g != 0, return(0));
  Ag = A / g; Bg = B / g;
  if (poldegree(Bg) == 0, return([a % A, A]));
  ge = gcdext(Ag, Bg); u = (((b - a) / g) * ge[1] / ge[3]) % Bg;
  return([(a + A*u) % (A*Bg), A*Bg]);
}
subsets_prod(S, m, R, start, acc, ~out) =
{
  if (m == 0, if (acc == R, listput(out, acc)); return(0));
  for (i = start, #S, if ((R % (acc*S[i])) == 0, subsets_prod(S, m-1, R, i+1, acc*S[i], ~out)));
  return(0);
}
handleF(F, f1, f2, f3, k, p, ~res) =
{
  my(o = Mod(1,p), R, S = List(), rts, out = List(), f);
  if (poldegree(F) != k || F % (f1*f2*f3) != 0, return(0));
  R = F / (f1*f2*f3);
  rts = polrootsmod(R, p);
  for (ir = 1, #rts, for (l = 1, p-1,
    f = l*(t - rts[ir])*o;
    if (f != f1 && f != f2 && f != f3 && ok(F, f), listput(S, f))));
  if (#S < k - 3, return(0));
  subsets_prod(Vec(S), k - 3, R, 1, o, ~out);
  if (#out > 0,
    res[1] += #out;
    if (F == pollead(F)*t^k, res[2] += #out, if (#res[3] < 4, listput(res[3], [f2, f3, F]))));
  return(0);
}
censlin(k, p) =
{
  my(o = Mod(1,p), f1 = o*t, A = f1^3 - f1, lins = List(), res = [0, 0, List()], c12, c, F0, L, f2, f3);
  for (a = 1, p-1, for (b = 0, p-1, my(f = o*(a*t + b)); if (f != f1, listput(lins, f))));
  for (i2 = 1, #lins,
    f2 = lins[i2]; c12 = crt2(f1, A, f2, f2^3 - f2);
    if (c12 != 0,
      for (i3 = i2+1, #lins,
        f3 = lins[i3]; c = crt2(c12[1], c12[2], f3, f3^3 - f3);
        if (c != 0,
          F0 = c[1]; L = c[2];
          if (poldegree(L) > k,
            handleF(F0, f1, f2, f3, k, p, ~res),
            if (poldegree(L) == k,
              for (s = 1, p-1, handleF(F0 + s*L, f1, f2, f3, k, p, ~res)),
              for (s1 = 1, p-1, for (s0 = 0, p-1, handleF(F0 + (s1*t + s0)*L, f1, f2, f3, k, p, ~res)))))))));
  return([k, p, res[1], res[2], Vec(res[3])]);
}
P = eval(getenv("PLIST")); K = eval(getenv("KLIST"));
for (ik = 1, #K, for (ip = 1, #P, my(t0 = getwalltime()); print(censlin(K[ik], P[ip]), "  ms=", getwalltime() - t0)));
