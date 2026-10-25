\\ aud13_A_rows.gp -- export my forms (aud13_A_forms_{0,1}.gp) decomposed for the relaxation rows.
\\ Each form X(moments): M1_0 -> S; X = c(S) + sum_L d_L(S) L + q M1_1^2, L in [M1_2, M2_0, M1_4, M2_2, M3_0];
\\ output line: name lin | c0 c1 c2 | d0 d1 (x5) | q      (exact rationals, reconstruction checked)
read("aud13_A_forms_0.gp"); read("aud13_A_forms_1.gp");
LIN = [M1_2, M2_0, M1_4, M2_2, M3_0];
dec(X) = {my(Y = simplify(subst(X, M1_0, 'SS)), d = vector(5, i, simplify(polcoef(Y, 1, LIN[i]))), q = simplify(polcoef(Y, 2, M1_1)), c = Y);
  for(i = 1, 5, c = subst(c, LIN[i], 0)); c = simplify(subst(c, M1_1, 0));
  if(Y != c + sum(i = 1, 5, d[i]*LIN[i]) + q*M1_1^2, error("structure"));
  if(#setminus(Set(variables(c)), Set(['SS])) || poldegree(c, 'SS) > 2, error("c"));
  for(i = 1, 5, if(#setminus(Set(variables(d[i])), Set(['SS])) || poldegree(d[i], 'SS) > 1, error("d")));
  if(#variables(q), error("q"));
  [vector(3, j, polcoef(c, j - 1, 'SS)), vector(5, i, [polcoef(d[i], 0, 'SS), polcoef(d[i], 1, 'SS)]), q]};
{
  my(fn = "aud13_A_rows.txt");
  foreach([0, 1], lin, my(FM = if(lin, FORMS1, FORMS0), names = ["Z0", "Z1", "PHI"]);
    my(L = concat([FM[1], FM[2], FM[3]], FM[4]), nm = concat(names, vector(8, n, Str("T", n))));
    for(i = 1, #L, my(r = dec(L[i]), s = Str(nm[i], " ", lin, " |"));
      for(j = 1, 3, s = Str(s, " ", r[1][j])); s = Str(s, " |");
      for(j = 1, 5, s = Str(s, " ", r[2][j][1], " ", r[2][j][2])); s = Str(s, " | ", r[3]);
      write(fn, s)));
  print("written ", fn);
}
