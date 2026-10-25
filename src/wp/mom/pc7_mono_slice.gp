\\ symmetric slice e1 = 0 (lam = 0, mu = 2/ep, ep = 1 + e2 in (0, 1]): exact monomial expansion of B1 and B1 + B5 (k = 7).
\\ Claim: every coefficient (constant, single-member w^a, pair u^a v^b, triple uvx) is a polynomial in ep with positive value on (0, 1].
read("q22_moments.gp");
{
  my(k = 7, B = bnd(k, 0, ep - 1, 0, 2/ep), P = vector(2*k, a, x^a + y^a + z^a));
  my(chk(Q, nm) = for(a = 1, 2*k, Q = subst(Q, eval(Str("p", a)), P[a]));
    my(bad = 0, cnt = 0);
    for(a = 0, 9, for(b = 0, a, for(c = 0, b, my(C = polcoef(polcoef(polcoef(Q, a, x), b, y), c, z));
      if(C != 0, cnt++; my(num = numerator(C), den = denominator(C));
        \\ positivity on (0, 1]: no root in (0, 1] and positive at 1/2, for numerator / denominator
        my(n0 = num / ep^valuation(num, ep), d0 = den / ep^valuation(den, ep)); my(sg = sign(subst(n0, ep, 1/2)) * sign(subst(d0, ep, 1/2)), nr = polsturm(n0, [0, 1]) + polsturm(d0, [0, 1]));
        if(sg <= 0 || nr > 0, bad++; print("  ", nm, " (", a, ",", b, ",", c, "): ", C))))));
    print(nm, ": ", cnt, " nonzero monomial coefficients, ", bad, " not positive on (0,1]"));
  chk(B[1], "B1"); chk(B[1] + B[5], "B1+B5"); chk(-B[4], "-B4"); chk(-(B[4] + B[8]), "-(B4+B8)");
  \\ the q-form identities recorded in the notebook (q_a = sum (mu w)^a, p_a = (ep/2)^a q_a)
  my(sb(Q) = for(a = 1, 9, Q = subst(Q, eval(Str("p", a)), (ep/2)^a * eval(Str("q", a)))); Q);
  my(S15 = sb(B[1] + B[5]), S48 = sb(B[4] + B[8]));
  print("B4+B8 identity: ", S48 + ep^2*(20*q1*ep^2 + q1*q2 + 9*q3) == 0);
  print("B1+B5 identity: ", S15 - (63/32*ep^6 + (3/2*q1^2 + 17/2*q2)*ep^4 + (q1*q3/2 + q2^2/4 + 17/4*q4)*ep^2) == 0);
}
