\\ Factor the other group order of each control prime. Input file IN (read via read()) defines CT = vector of
\\ [tagcode, x, p] with tagcode 1 = FFc (N = p+1), 2 = LFc (N = p-1). Output: CTL <tag> x p <fac(N)>|FAILN.
fstr(fa) = if (#fa~ == 0, "1", strjoin(vector(#fa~, i, Str(fa[i, 1], "^", fa[i, 2])), "*"));
tfac(N) = my(r = iferr(alarm(TL, factor(N)), E, E)); if (type(r) == "t_ERROR", 0, r);
read(IN);
{
  my(t0 = getabstime(), nf = 0);
  for (i = 1, #CT,
    my(c = CT[i][1], x = CT[i][2], p = CT[i][3], N = if (c == 1, p + 1, p - 1), fN = tfac(N), tag = if (c == 1, "FFc", "LFc"));
    if (fN == 0, nf++; write(OUT, "CTL ", tag, " ", x, " ", p, " FAILN"); next);
    write(OUT, "CTL ", tag, " ", x, " ", p, " ", fstr(fN)));
  write(OUT, "#END ctl n ", #CT, " fail ", nf, " ms ", getabstime() - t0);
}
quit;
