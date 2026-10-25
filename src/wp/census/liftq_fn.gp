\\ Hensel test for an F_p Williams family (members of degree <= 2, fixed degree pattern).
\\ Unknowns: free coefficients of members (f1 = alpha x^2 + gamma with alpha fixed, gamma free;
\\ others all coefficients up to their degree) and cofactors H_j with F - f_j = (f_j^3 - f_j) H_j.
lifttestq(fam, p, N) =
{
  my(k = #fam, degs = vector(k, j, poldegree(fam[j])), D, V = List(), x0 = List(), mem = vector(k), H = vector(k), cnt = 0, v, F, eqs = List(), J, xs, q, E, r, y, Jm, alpha, hj, Hp, vv, nv);
  D = vecsum(degs);
  alpha = lift(pollead(fam[1]));
  \\ member 1
  cnt++; v = varlower(Str("u", cnt)); listput(V, v); listput(x0, lift(polcoef(fam[1], 0)));
  mem[1] = alpha*x^2 + v;
  for (j = 2, k,
    mem[j] = 0;
    for (i = 0, degs[j], cnt++; v = varlower(Str("u", cnt)); listput(V, v); listput(x0, lift(polcoef(fam[j], i))); mem[j] += v*x^i));
  F = prod(j = 1, k, mem[j]);
  for (j = 1, k,
    Hp = (prod(i = 1, k, fam[i]) - fam[j]) / (fam[j]^3 - fam[j]);
    H[j] = 0;
    for (i = 0, D - 3*degs[j], cnt++; v = varlower(Str("u", cnt)); listput(V, v); listput(x0, lift(polcoef(Hp, i))); H[j] += v*x^i);
    my(ex = F - mem[j] - (mem[j]^3 - mem[j])*H[j]);
    for (i = 0, D, listput(eqs, polcoef(ex, i, x))));
  V = Vec(V); xs = Vec(x0); eqs = Vec(eqs); nv = #V;
  J = matrix(#eqs, nv, i, j, deriv(eqs[i], V[j]));
  q = p;
  for (n = 1, N - 1,
    E = vector(#eqs, i, substvec(eqs[i], V, xs));
    if (#select(z -> z % q != 0, E) > 0, return(["not a solution mod", q]));
    r = -vector(#eqs, i, (E[i] / q) % p)~;
    Jm = matrix(#eqs, nv, i, j, substvec(J[i,j], V, xs) % p);
    y = matsolvemod(Jm, p, r);
    if (y == 0, return(["no lift at level", n + 1, "unknowns", nv, "rankJ", matrank(Jm * Mod(1,p))]));
    xs = xs + q * y~;
    q *= p);
  ["lifted to p^", N, "unknowns", nv, "rankJ", matrank(Jm * Mod(1,p))];
}
