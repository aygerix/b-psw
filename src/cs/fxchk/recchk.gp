\\ Adversarial FX check: recompute sampled pool records "D p t j1 j2 o1 o2" from scratch.
\\ t = 0 inert / 1 split; j1 = v2(ord_p 2), o1 = odd(ord_p 2); j2 = v2(z_D(p)), o2 = odd(z_D(p)); z via x^k mod (x^2 - x + Q).
Ux(k, Q, p) = polcoeff(lift(lift(Mod(Mod(1, p) * x, Mod(1, p) * (x^2 - x + Q))^k)), 1);
zr(p, Q, m) = { my(z = m, f = factor(m)[, 1]); if (Ux(m, Q, p) != 0, error("U_m != 0 at ", p)); for (i = 1, #f, while (z % f[i] == 0 && Ux(z / f[i], Q, p) == 0, z /= f[i])); z; }
read("recsample.gp");
{
  my(bad = 0, nb = 0, maxp = 0);
  for (i = 1, #R, my(r = R[i], D = r[1], p = r[2], Q = (1 - D) / 4, js, o, z, t);
    if (!isprime(p), bad++; print("not prime ", r); next);
    js = kronecker(D, p); t = if (js == 1, 1, 0);
    o = znorder(Mod(2, p)); z = zr(p, Q, p - js); maxp = max(maxp, p);
    my(ok = (js != 0) && t == r[3] && valuation(o, 2) == r[4] && valuation(z, 2) == r[5] && (o >> valuation(o, 2)) == r[6] && (z >> valuation(z, 2)) == r[7]);
    if (t == 1 && gcd(r[6], r[7]) != 1, ok = 0);
    if (!ok, bad++; if (nb++ < 10, print("MISMATCH ", r, " recomputed j1 ", valuation(o, 2), " o1 ", o >> valuation(o, 2), " j2 ", valuation(z, 2), " o2 ", z >> valuation(z, 2)))));
  print("records checked ", #R, ", mismatches ", bad, ", max p ", maxp);
}
