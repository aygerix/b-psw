\\ HD test 1 (exact): linear common-root k = 5, f_i = 1 + a_i x, members 1..3 Lucas-good:
\\ E_j: prod_{i != j} (a_j - 2 a_i) - a_j^4 = 0. Normalise a1 = 1, fix a2 = s rational, solve for a3, a4, a5 by resultants.
E(v, j) = prod(i = 1, #v, if (i == j, 1, v[j] - 2*v[i])) - v[j]^(#v - 1);
{
  for (sn = -12, 12, my(s = sn/4); if (s == 0 || s == 1, next);
    my(v = [1, s, 'a3, 'a4, 'a5]);
    my(e1 = E(v, 1), e2 = E(v, 2), e3 = E(v, 3));
    \\ e1, e2 linear in a5? e1: (1-2s)(1-2a3)(1-2a4)(1-2a5) - 1 ; linear in a5. Eliminate a5 from e1.
    my(r12 = polresultant(e1, e2, 'a5), r13 = polresultant(e1, e3, 'a5));
    my(R = polresultant(r12, r13, 'a4));
    R = R / content(R);
    my(fa = factor(R), nreal = 0, sols = List());
    for (i = 1, #fa~, my(P = fa[i,1]); if (poldegree(P, 'a3) < 1, next);
      my(rts = polrootsreal(P));
      for (t = 1, #rts, my(a3 = rts[t]);
        \\ back-substitute: solve r12(a3, a4) = 0 and r13 = 0 common root a4
        my(g1 = subst(r12, 'a3, a3), g2 = subst(r13, 'a3, a3));
        my(ra = polroots(g1));
        for (u = 1, #ra, my(a4 = ra[u]);
          if (abs(imag(a4)) > 1e-20 || abs(subst(g2, 'a4, a4)) > 1e-15 * (1 + abs(subst(g2,'a4,a4+1))), next);
          a4 = real(a4);
          my(l5 = substvec(e1, ['a3, 'a4], [a3, a4]));
          if (poldegree(l5, 'a5) != 1, next);
          my(a5 = -polcoeff(l5, 0, 'a5)/polcoeff(l5, 1, 'a5));
          my(w = [1, s, a3, a4, a5]);
          my(md = vecmin(concat(vector(5, i, vector(5, j, if (i == j, 1e9, abs(w[i] - w[j])))))), mz = vecmin(apply(abs, w)));
          if (md < 1e-8 || mz < 1e-8, next);
          my(chk = vector(5, j, E(w, j) / w[j]^4));
          listput(sols, [w, chk]))));
    if (#sols, print("s = ", s, ": ", #sols, " real distinct solutions; first: a = ", sols[1][1], "  E_j/a_j^4 = ", sols[1][2]))
  );
}
