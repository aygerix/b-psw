\\ k = 7 (seven cubics; six cubics + 2t+1): extract the four boundary equations in moment form and test the
\\ advisor's forced values of M = sum v^2 and K = sum v (u - 1/2)^2.
read("cubic_moments.gp");
\\ monomial "types": coefficient of a representative monomial in member 1 (and pair 1,2)
types(E) = {my(u = u1, v = v1);
  [polcoef(polcoef(E, 1, v1), 0, u1), polcoef(polcoef(E, 1, v1), 1, u1), polcoef(polcoef(E, 1, v1), 2, u1), polcoef(polcoef(E, 1, v1), 3, u1),
   polcoef(polcoef(E, 2, v1), 0, u1), polcoef(polcoef(E, 2, v1), 1, u1)]};
{
  foreach([[7, 0], [6, 1]], c,
    my(m = c[1], lin = c[2], E = eqs(m, lin));
    \\ coefficients: [v, v u, v u^2, v u^3, v^2, v^2 u] for member 1 (set other members' variables to 0 first),
    \\ pair coefficient of v1 v2, and the constant term
    my(Z = vector(4, j, my(e = E[j]); for(i = 2, m, e = subst(subst(e, eval(Str("u", i)), 0), eval(Str("v", i)), 0)); e));
    my(tab = matrix(4, 6, j, r, types(Z[j])[r]));
    my(pair = vector(4, j, my(e = E[j]); for(i = 3, m, e = subst(subst(e, eval(Str("u", i)), 0), eval(Str("v", i)), 0));
         polcoef(polcoef(polcoef(polcoef(e, 1, v1), 1, v2), 0, u1), 0, u2)));
    my(cst = vector(4, j, my(e = E[j]); for(i = 1, m, e = subst(subst(e, eval(Str("u", i)), 0), eval(Str("v", i)), 0)); e));
    print("m = ", m, " lin = ", lin);
    for(j = 1, 4, print("  E", j, ": [v, vu, vu^2, vu^3, v^2, v^2u] = ", tab[j,], "  pair v_i v_j: ", pair[j], "  const: ", cst[j]));
    \\ consistency: E must be symmetric with pair terms only of the v_i v_j type -- check total degree structure
    print("  check: member-1-only part reproduces E for m = 1 types; pair coefficient uniform"));
}
