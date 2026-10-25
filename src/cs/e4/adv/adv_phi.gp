\\ Adversarial re-implementation of the E4 forced-prime enumeration (written independently of fic.gp).
\\ Globals (set in driver before read): MODE ("F" or "L"), XMIN, XMAX (exponent range), BLO, BHI (band on the
\\ size proxy: phi(e) for F, 0.6942*phi(f) for L; keep exponent iff BLO < proxy <= BHI), TL (seconds per factor),
\\ OUT (output file).
\\ Phi_e(2) is computed as prod_{d|e} (2^d-1)^mu(e/d)   (fic.gp uses polcyclo(e,2)).
\\ Phi*_f is computed as the homogenised cyclotomic value Phi_f(alpha,beta) in Z[t]/(t^2-t-1)   (fic.gp uses
\\ prod_{d|f} F_d^mu(f/d)); the Python checker recomputes it the other way and compares.
\\ Output lines:
\\   PHI <MODE> x <fac>|FAIL        full factorisation of Phi_x (all prime factors, exponents)
\\   REC <MODE> x p <fac(N)>|FAILN  every prime factor p of Phi_x (no filter here), N = p+1 (F) or p-1 (L)
\\ Filtering (inert, primitive, p > 5) is done ONLY in the Python checker.
fstr(fa) = if (#fa~ == 0, "1", strjoin(vector(#fa~, i, Str(fa[i, 1], "^", fa[i, 2])), "*"));
tfac(N) = my(r = iferr(alarm(TL, factor(N)), E, E)); if (type(r) == "t_ERROR", 0, r);
phiF(e) = my(c = 1); fordiv(e, d, c *= (2^d - 1)^moebius(e / d)); c;
phiL(f) =
{
  my(P = polcyclo(f, 'y), n = poldegree(P, 'y), m = 't^2 - 't - 1, a = Mod('t, m), b = Mod(1 - 't, m), s = Mod(0, m),
     ap = vector(n + 1), bp = vector(n + 1));
  ap[1] = Mod(1, m); bp[1] = Mod(1, m);
  for (k = 1, n, ap[k + 1] = ap[k] * a; bp[k + 1] = bp[k] * b);
  for (k = 0, n, s += polcoeff(P, k, 'y) * ap[k + 1] * bp[n - k + 1]);
  s = lift(s);
  if (poldegree(s, 't) > 0, error("Phi*_f irrational at f = ", f));
  abs(polcoeff(s, 0, 't));
}
proxy(x) = if (MODE == "F", eulerphi(x), 0.6942 * eulerphi(x));
{
  my(t0 = getabstime(), nx = 0, nfail = 0, nrec = 0, nfailn = 0);
  for (x = XMIN, XMAX,
    if (x % 5 == 0, next);
    if (MODE == "L" && x < 3, next);
    my(pr = proxy(x));
    if (pr <= BLO || pr > BHI, next);
    nx++;
    my(c = if (MODE == "F", phiF(x), phiL(x)), fa = tfac(c));
    if (fa == 0, nfail++; write(OUT, "PHI ", MODE, " ", x, " FAIL"); next);
    write(OUT, "PHI ", MODE, " ", x, " ", fstr(fa));
    for (i = 1, #fa~,
      my(p = fa[i, 1], N = if (MODE == "F", p + 1, p - 1), fN);
      if (p <= 5, next);
      \\ only factor N for candidates the checker could keep (cheap pre-filter; the checker re-derives the set)
      if (kronecker(5, p) != -1 || x % p == 0, next);
      fN = tfac(N);
      nrec++;
      if (fN == 0, nfailn++; write(OUT, "REC ", MODE, " ", x, " ", p, " FAILN"); next);
      write(OUT, "REC ", MODE, " ", x, " ", p, " ", fstr(fN))));
  write(OUT, "#END MODE ", MODE, " X ", XMIN, "..", XMAX, " band (", BLO, ",", BHI, "] TL ", TL, " exps ", nx,
        " phifail ", nfail, " recs ", nrec, " Nfail ", nfailn, " ms ", getabstime() - t0);
}
quit;
