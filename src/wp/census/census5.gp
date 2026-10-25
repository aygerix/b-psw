\\ Complete F_p census of 5-member Williams families with members of degree <= 2
\\ (coprime AND non-coprime): F = prod f_j, F == f_j mod (f_j^3 - f_j) for all j.
ok(F,f) = (F - f) % (f^3 - f) == 0;
monicdivs(R, dmax) =
{
  my(fa = factor(R), L = [1], M, q);
  for (i = 1, #fa~,
    M = List();
    for (a = 1, #L,
      for (e = 0, fa[i,2],
        q = L[a] * fa[i,1]^e;
        if (poldegree(q) <= dmax, listput(M, q))));
    L = Vec(M));
  select(x -> poldegree(x) >= 1, L);
}
canon(v) = vecsort(apply(x -> Str(lift(x)), v));
ismono(v) =
{
  my(m = 1);
  for (a = 2, 5,
    if (poldegree(v[a]) != poldegree(v[1]), m = 0,
      if (v[a]*pollead(v[1]) != v[1]*pollead(v[a]), m = 0)));
  m;
}
census(p) =
{
  my(o = Mod(1,p), nr = 2, F1 = List(), fams = Map(), nmono = 0, nnon = 0, cand = 0, ex = List());
  my(d1, A, f2, B, g, Ag, Bg, uu, F0, L, Fs, R, S, n, v, key, ge, F);
  while (issquare(Mod(nr,p)), nr++);
  for (g0 = 0, p-1, listput(F1, o*t^2 + g0); listput(F1, nr*o*t^2 + g0));
  listput(F1, o*t);
  for (i1 = 1, #F1,
    my(f1 = F1[i1]);
    d1 = poldegree(f1); A = f1^3 - f1;
    forvec(c = vector(d1+1, i, [0, p-1]),
      f2 = o*Pol(c, t);
      if (poldegree(f2) >= 1 && f2 != f1,
        B = f2^3 - f2; g = gcd(A, B);
        if ((f2 - f1) % g == 0,
          Ag = A / g; Bg = B / g; L = A * Bg;
          if (poldegree(Bg) == 0, uu = 0,
            ge = gcdext(Ag, Bg); uu = (((f2 - f1) / g) * ge[1] / ge[3]) % Bg);
          F0 = (f1 + A * uu) % L;
          Fs = List();
          if (poldegree(L) > 10,
            if (poldegree(F0) <= 10, listput(Fs, F0)),
            for (s = 0, p-1, listput(Fs, F0 + s*L)));
          for (iF = 1, #Fs,
            F = Fs[iF];
            if (F != 0 && F % (f1*f2) == 0,
              R = F / (f1*f2);
              if (poldegree(R) >= 3,
                cand++;
                S = List();
                my(dv = monicdivs(R, d1));
                for (a = 1, #dv, for (l = 1, p-1,
                  my(f = l*dv[a]);
                  if (f != f1 && f != f2 && ok(F, f), listput(S, f))));
                n = #S;
                for (i = 1, n, for (j = i+1, n, for (k = j+1, n,
                  if (S[i]*S[j]*S[k] == R,
                    v = [f1, f2, S[i], S[j], S[k]]; key = canon(v);
                    if (!mapisdefined(fams, key),
                      mapput(fams, key, 1);
                      if (ismono(v), nmono++,
                        nnon++; if (#ex < 5, listput(ex, v)))))))))))))));
  [p, cand, nmono, nnon, Vec(ex)];
}
P = eval(getenv("PLIST")); for (ip = 1, #P, my(p = P[ip], t0 = getwalltime()); print(census(p), "  ms=", getwalltime() - t0));
