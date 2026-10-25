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
o19 = Mod(1,19);
famA = [o19*x^2, o19*(3*x+13), o19*(13*x^2+2), o19*(18*x^2+17), o19*(16*x+13)];
print("F19 (2,2,2,1,1): ", lifttestq(famA, 19, 5));
o17 = Mod(1,17);
famB = [o17*x^2, o17*(x^2+6), o17*(8*x^2+7), o17*(10*x^2+13), o17*(13*x^2+5)];
print("F17 linear-in-x^2: ", lifttestq(famB, 17, 5));
o11 = Mod(1,11); z = Mod(3,11); \\ 3 has order 5 mod 11
famC = vector(5, j, o11*z^(j-1)*(x^2 + 4));
print("F11 monomial control: ", lifttestq(famC, 11, 5));
