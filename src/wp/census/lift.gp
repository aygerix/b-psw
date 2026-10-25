\\ Hensel test: does an F_p Williams family of linear members lift to Z/p^N ?
\\ unknowns x = [a2,b2,...,a5,b5]; f1 = t; equations F((e - b_j)/a_j) = e, e = +-1, j = 1..5.
vars = [a2,b2,a3,b3,a4,b4,a5,b5];
fl = [t, a2*t+b2, a3*t+b3, a4*t+b4, a5*t+b5];
Fp = prod(j = 1, 5, fl[j]);
eqs = List();
{
  for (j = 1, 5, my(aj = pollead(fl[j], t), bj = polcoef(fl[j], 0, t));
    foreach([1, -1], e, listput(eqs, subst(Fp, t, (e - bj)/aj) - e)));
}
eqs = Vec(eqs);
J = matrix(#eqs, #vars, i, j, deriv(eqs[i], vars[j]));
ev(ex, x) = substvec(ex, vars, x);
lifttest(x0, p, N) =
{
  my(x = x0, q = p, E, Jm, r, y);
  for (n = 1, N - 1,
    E = vector(#eqs, i, Mod(ev(eqs[i], x), q*p));
    if (#select(z -> lift(z) % q != 0, E) > 0, return(["not a solution mod", q]));
    r = -vector(#eqs, i, lift(E[i]) / q)~ * Mod(1, p);
    Jm = matrix(#eqs, #vars, i, j, Mod(ev(J[i,j], x), p));
    y = iferr(matsolvemod(lift(Jm), p, lift(r)), err, 0);
    if (y == 0, return(["no lift at level", n + 1, "rankJ", matrank(Jm)]));
    x = x + q * lift(y)~;
    q *= p);
  ["lifted to p^", N, "rankJ", matrank(Jm)];
}
\\ F_17 family {t, t+6, 8t+7, 10t+13, 13t+5}
print(lifttest([1,6,8,7,10,13,13,5], 17, 6));
print(lifttest([1,11,8,10,10,4,13,12], 17, 6));
