\r census_all.gp
\r liftq_fn.gp
P = eval(getenv("PLIST"));
{
for (ip = 1, #P, my(p = P[ip], r = census(p), fams = r[5], nl = 0, res, lifted = List(), lv = List());
  for (i = 1, #fams,
    res = lifttestq(apply(f -> subst(lift(f), t, x) * Mod(1,p), fams[i]), p, 4);
    listput(lv, if (res[1] == "no lift at level", res[2], 99));
    if (res[1] == "lifted to p^", nl++; listput(lifted, fams[i])));
  print([p, "cand", r[2], "mono", r[3], "nonmono hits", #fams, "lifting to p^4", nl, "fail levels", Vec(lv), Vec(lifted)]));
}
