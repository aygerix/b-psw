\\ E4: forced-index coupling (FIC) test.  See PREREG.txt.
\\ Globals to set before read(): EMAX, FMAX, KCTL, TLIM, SEED, OUT (record file), EMIN, FMIN (defaults 2, 3).
\\ Records (one line per prime), written to OUT:
\\  tag e_or_f bits(p) bits(No) logidx Elog Vlog hq ehq  o3 o10 o100 o1000  p3 p10 p100 p1000  U200 P200 U1500 P1500 U1e4 P1e4
\\  tag in {FF, FFc, LF, LFc}: FF = p | Phi_e(2) (Fermat-forced; other side = gamma in T, order p+1),
\\  LF = p | Phi*_f (Lucas-forced; other side = 2 in F_p^*, order p-1); 'c' = matched control.
\\ (set parisizemax in the driver: changing it inside read() aborts the read)
if (type(EMIN) == "t_POL", EMIN = 2);
if (type(FMIN) == "t_POL", FMIN = 3);
if (type(PHIBITS) == "t_POL", PHIBITS = 10^9);
if (type(PHILO) == "t_POL", PHILO = 0);
setrand(SEED);

fibmod(m, p) = polcoeff(lift(Mod(Mod(1, p) * 'x, 'x^2 - 'x - 1)^m), 1);
oddpart(n) = n >> valuation(n, 2);
sfac(N) =
{
  my(r = alarm(TLIM, factor(N)));
  if (type(r) == "t_ERROR", 0, r);
}

\\ Fibonacci entry point of p, given N = multiple of z(p) and fa = factor(N)
zent(p, N, fa) =
{
  my(n = N);
  for (i = 1, #fa~, my(q = fa[i, 1]);
    for (j = 1, fa[i, 2], if (fibmod(n / q, p) == 0, n /= q, break)));
  n;
}

\\ statistics of the odd index of an element of order o in a cyclic group of order N (fa = factor(N))
\\ returns [bits(No), logidx, Elog, Vlog, hq, ehq, obs tails (4), pred tails (4), U200, P200, U1500, P1500, U1e4, P1e4]
TT = [3, 10, 100, 1000];
TU = [200, 1500, 10000];
stats(N, fa, o) =
{
  my(No = oddpart(N), oo = oddpart(o), idx = No / oo, Elog = 0., Vlog = 0., hq = 0, ehq = 0., cst = 1., ot = vector(4), pt = vector(4), res);
  if (No % oo, error("order does not divide group order"));
  for (i = 1, #fa~, my(q = fa[i, 1], a = fa[i, 2], s = 0., s2 = 0., w);
    if (q == 2, next);
    for (j = 1, a, w = if (j < a, q^-j * (1 - 1. / q), 1. * q^-a); s += j * w; s2 += j^2 * w);
    Elog += s * log(q); Vlog += (s2 - s^2) * log(q)^2; cst *= (1 - 1. / q);
    if (q <= 50, hq += (idx % q == 0); ehq += 1. / q));
  \\ P(idx = d) = phi(No/d)/No = (1/d) * prod_{q | No} (1-1/q) * prod_{q: v_q(d) = v_q(No)} (1-1/q)^-1
  for (k = 1, 4, my(t = TT[k], s = 0.);
    forstep (d = 1, t - 1, 2,
      if (No % d == 0,
        my(c = cst / d);
        for (i = 1, #fa~, my(q = fa[i, 1]); if (q != 2 && valuation(d, q) == fa[i, 2], c /= (1 - 1. / q)));
        s += c));
    pt[k] = 1 - s; ot[k] = (idx >= t));
  res = [#binary(No), log(idx), Elog, Vlog, hq, ehq, ot[1], ot[2], ot[3], ot[4], pt[1], pt[2], pt[3], pt[4]];
  for (k = 1, 3, my(T = TU[k], sT = 1);
    for (i = 1, #fa~, my(q = fa[i, 1]); if (q != 2 && q <= T, sT *= q^fa[i, 2]));
    res = concat(res, [(sT % oo == 0), 1. * sT / No]));
  res;
}

emit(tag, ef, p, st) = write(OUT, tag, " ", ef, " ", #binary(p), " ", strjoin(apply(x -> if (type(x) == "t_REAL", strjoin(strsplit(Strprintf("%.6g", x), " "), ""), Str(x)), st), " "));

nF = 0; nFfail = 0; nFN = 0; nL = 0; nLfail = 0; nLN = 0; nC = 0; nCfail = 0; nCskip = 0;

\\ ---- Fermat-forced: inert primitive prime factors of Phi_e(2); other side = gamma, group order p+1
ctlF(p, e) =
{
  my(B = #binary(p), lo = 2^(B - 1), hi = 2^B - 1);
  if ((hi - lo) / lcm(e, 5) < 64 * B, nCskip++; return);    \\ class too sparse at this size: no control
  for (k = 1, KCTL,
    my(r5 = if (random(2), 2, 3), cl = chinese(Mod(1, e), Mod(r5, 5)), pc, N, fN, z);
    pc = iferr(randomprime([lo, hi], cl), E, 0);
    if (pc == 0, nCfail++; next);
    N = pc + 1; fN = sfac(N); if (fN == 0, nCfail++; next);
    z = zent(pc, N, fN); nC++;
    emit("FFc", e, pc, stats(N, fN, z)));
}
runF() =
{
  for (e = EMIN, EMAX,
    if (e % 5 == 0, next);                       \\ primitive p == 1 mod 5: split, never inert
    if (eulerphi(e) > PHIBITS || eulerphi(e) <= PHILO, next);   \\ size band: Phi_e(2) has ~phi(e) bits
    my(c = polcyclo(e, 2), fa = sfac(c));
    if (fa == 0, nFfail++; write(OUT, "#FAILPHI e ", e, " digits ", #Str(c)); next);
    for (i = 1, #fa~, my(p = fa[i, 1], N, fN, z);
      if (p <= 5 || e % p == 0 || kronecker(5, p) != -1, next);
      N = p + 1; fN = sfac(N); if (fN == 0, nFN++; write(OUT, "#FAILN FF ", e, " ", p); next);
      z = zent(p, N, fN); nF++;
      emit("FF", e, p, stats(N, fN, z));
      ctlF(p, e)));
}

\\ ---- Lucas-forced: inert primitive prime factors of Phi*_f; other side = 2, group order p-1
ctlL(p, f) =
{
  my(B = #binary(p), lo = 2^(B - 1), hi = 2^B - 1);
  if ((hi - lo) / lcm(f, 5) < 64 * B, nCskip++; return);    \\ class too sparse at this size: no control
  for (k = 1, KCTL,
    my(r5 = if (random(2), 2, 3), cl = chinese(Mod(-1, f), Mod(r5, 5)), pc, N, fN, o);
    pc = iferr(randomprime([lo, hi], cl), E, 0);
    if (pc == 0, nCfail++; next);
    N = pc - 1; fN = sfac(N); if (fN == 0, nCfail++; next);
    o = znorder(Mod(2, pc), [N, fN]); nC++;
    emit("LFc", f, pc, stats(N, fN, o)));
}
runL() =
{
  for (f = FMIN, FMAX,
    if (f % 5 == 0, next);                       \\ inert p has 5 !| p+1, so f = z(p) is prime to 5
    if (0.6942 * eulerphi(f) > PHIBITS || 0.6942 * eulerphi(f) <= PHILO, next);   \\ size band: Phi*_f ~0.694 phi(f) bits
    my(c = 1, fa);
    fordiv(f, d, c *= fibonacci(d)^moebius(f / d));
    if (denominator(c) != 1, error("Phi*_f not integral"));
    fa = sfac(c);
    if (fa == 0, nLfail++; write(OUT, "#FAILPHI f ", f, " digits ", #Str(c)); next);
    for (i = 1, #fa~, my(p = fa[i, 1], N, fN, o);
      if (p <= 5 || f % p == 0 || kronecker(5, p) != -1, next);
      if (fibmod(f, p) != 0, error("not a divisor of F_f"));
      N = p - 1; fN = sfac(N); if (fN == 0, nLN++; write(OUT, "#FAILN LF ", f, " ", p); next);
      o = znorder(Mod(2, p), [N, fN]); nL++;
      emit("LF", f, p, stats(N, fN, o));
      ctlL(p, f)));
}

t0 = getabstime();
runF(); t1 = getabstime();
runL(); t2 = getabstime();
{
write(OUT, "#DONE EMAX ", EMAX, " FMAX ", FMAX, " KCTL ", KCTL, " TLIM ", TLIM, " PHIBITS ", PHIBITS, " PHILO ", PHILO, " nF ", nF, " nL ", nL, " nC ", nC, " failPhiF ", nFfail, " failNF ", nFN, " failPhiL ", nLfail, " failNL ", nLN, " failC ", nCfail, " skipC ", nCskip, " msF ", t1 - t0, " msL ", t2 - t1);
}
quit;
