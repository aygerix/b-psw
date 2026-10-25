\\ The k = 5 two-dead curve: a1 = 1, a5 = M; E_2 = E_3 = E_4 = 0.  Find the component through the certified point, with M symbolic.
E(v, j) = prod(i = 1, #v, if (i == j, 1, v[j] - 2*v[i])) - v[j]^(#v - 1);
{
  my(v = [1, 'a2, 'a3, 'a4, 'M]);
  my(e2 = E(v, 2), e3 = E(v, 3), e4 = E(v, 4));
  my(r23 = polresultant(e2, e3, 'a4), r24 = polresultant(e2, e4, 'a4));
  print("deg r23 in a3: ", poldegree(r23, 'a3), "  deg r24 in a3: ", poldegree(r24, 'a3));
  my(R = polresultant(r23, r24, 'a3));
  print("eliminant: deg_a2 ", poldegree(R, 'a2), ", deg_M ", poldegree(R, 'M));
  my(fa = factor(R));
  for (i = 1, #fa~, my(P = fa[i, 1]);
    print("factor ", i, ": mult ", fa[i, 2], "  deg_a2 ", poldegree(P, 'a2), "  deg_M ", poldegree(P, 'M), "   value at (1.99924256189045, 135): ", subst(subst(P, 'M, 135), 'a2, 1.99924256189045)));
  for (i = 1, #fa~, my(P = fa[i, 1]); if (poldegree(P, 'a2) == 3 || poldegree(P, 'M) <= 3, print("factor ", i, " = ", P)));
}
