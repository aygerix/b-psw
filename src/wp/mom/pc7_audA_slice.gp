\\ pc7_audA_slice: independent re-run of the slice claim (e1 = 0) with FOUR variables, to confirm that no monomial with
\\ >= 4 parts occurs (so 3 variables capture the full monomial-symmetric expansion for k = 7 members), and positivity on (0,1].
read("q22_moments.gp");
{
  my(k = 7, B = bnd(k, 0, ep - 1, 0, 2/ep), P = vector(2*k, a, x^a + y^a + z^a + v^a));
  my(Q = B[1]); for(a = 1, 2*k, Q = subst(Q, eval(Str("p", a)), P[a]));
  print("variables in B1: ", variables(B[1]));
  \\ total degree in the p's:
  my(cnt = 0, bad = 0, four = 0, maxdeg = 0);
  for(a = 0, 14, for(b = 0, a, for(cc = 0, b, for(dd = 0, cc,
    my(C = polcoef(polcoef(polcoef(polcoef(Q, a, x), b, y), cc, z), dd, v));
    if(C != 0, cnt++; maxdeg = max(maxdeg, a);
      if(dd > 0, four++);
      my(num = numerator(C), den = denominator(C));
      my(n0 = num / ep^valuation(num, ep), d0 = den / ep^valuation(den, ep));
      my(pos = sign(subst(n0, ep, 1/2)) * sign(subst(d0, ep, 1/2)) > 0 && polsturm(n0, [0, 1]) == 0 && polsturm(d0, [0, 1]) == 0);
      if(!pos, bad++);
      print("  m(", a, ",", b, ",", cc, ",", dd, "): ", factor(C), "   positive on (0,1]: ", pos))))));
  print("B1: ", cnt, " nonzero monomial coefficients, ", four, " with 4 parts, ", bad, " not positive, max single exponent ", maxdeg);
  \\ random numerical control: B1 at random positive w's (7 members) by direct power sums
  my(mn = oo);
  for(tr = 1, 2000, my(epv = random(1000)/1000. + 0.001, wv = vector(k, j, 10.^(random(7) - 3) * random(1000)/1000.));
    my(val = subst(B[1], ep, epv)); for(a = 1, 2*k, val = subst(val, eval(Str("p", a)), sum(j = 1, k, wv[j]^a)));
    mn = min(mn, val));
  print("min B1 over 2000 random (ep, positive w): ", mn);
}
