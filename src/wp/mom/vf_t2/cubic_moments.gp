\\ Independent check of the advisor's Laurent/reflection argument for real cubic Williams families.
\\ Normal form f_i = 2t + 1 + t(t+1)(a_i t + c_i); u_i = c_i/a_i, v_i = 1/a_i.
\\ x_i = 1/(f_i - 1) = v_i s^3 / (1 + (1+u_i) s + (u_i + 2 v_i) s^2), s = 1/t.  Linear member 2t+1: x_0 = s/2.
\\ V = t^(k-1) [prod(1 + x_i) - prod x_i] is a polynomial; its value and derivative at t = -1 are forced:
\\   V(-1) = 2^(1-k), V'(-1) = (1-k) 2^(1-k).  Reflection: (a, c) -> (a, a - c), i.e. u -> 1 - u, v -> v.
\\ For m symbolic cubic members (+ optional linear one) compute E1 = V(-1) - 2^(1-k), E2 = V'(-1) - (1-k)2^(1-k)
\\ for the family and for its reflection, as polynomials in (u_i, v_i).
\\ truncated power series in s handled as polynomials (no t_SER: avoids variable-priority problems)
trunc(P, n) = sum(j = 0, n - 1, polcoef(P, j, s) * s^j);
xser(u, v, n) = {my(c = vector(n + 1)); c[1] = 1; if(n >= 1, c[2] = -(1 + u));
  for(j = 2, n, c[j+1] = -(1 + u)*c[j] - (u + 2*v)*c[j-1]);
  trunc(v * s^3 * sum(j = 0, n, c[j+1]*s^j), n + 1)};
Vpoly(m, lin, U, W) = {
  my(k = m + lin, G = 1);
  for(i = 1, m, G = trunc(G * (1 + xser(U[i], W[i], k)), k));
  if(lin, G = trunc(G * (1 + s/2), k));
  sum(j = 0, k - 1, polcoef(G, j, s) * t^(k-1-j))};
eqs(m, lin) = {
  my(k = m + lin, U = vector(m, i, eval(Str("u", i))), W = vector(m, i, eval(Str("v", i))), Ur = vector(m, i, 1 - U[i]),
     V = Vpoly(m, lin, U, W), Vr = Vpoly(m, lin, Ur, W), c0 = 2^(1-k), c1 = (1-k)*2^(1-k));
  [subst(V, t, -1) - c0, subst(deriv(V, t), t, -1) - c1, subst(Vr, t, -1) - c0, subst(deriv(Vr, t), t, -1) - c1]};
\\ moment helpers
mom(m, f) = sum(i = 1, m, f(eval(Str("u", i)), eval(Str("v", i))));
\\ find constants lam with sum lam_j E_j == target (exact); returns lam or 0
certify(E, target) = {
  my(lam = vector(4, j, eval(Str("L", j))), R = target - sum(j = 1, 4, lam[j]*E[j]), vars, M, rhs, cf, sol);
  \\ collect all coefficients of R (as polynomial in u_i, v_i) -- linear in L1..L4 -- and solve
  cf = List(); my(P = R);
  \\ flatten: substitute random rationals is NOT exact; instead extract coefficients via variables list
  my(vl = variables(P)); vl = select(z -> !(z == L1 || z == L2 || z == L3 || z == L4), vl);
  my(rec(Q, idx) = if(idx > #vl, listput(cf, Q), for(d = 0, poldegree(Q, vl[idx]), rec(polcoef(Q, d, vl[idx]), idx + 1))));
  rec(P, 1);
  M = matrix(#cf, 4, r, c, polcoef(cf[r], 1, lam[c])); rhs = vector(#cf, r, -subst(subst(subst(subst(cf[r], L1, 0), L2, 0), L3, 0), L4, 0))~;
  sol = iferr(matsolve(M[1..min(#cf, 40), ], rhs[1..min(#cf, 40)]), e, 0);
  \\ least-squares-free: test candidate on all rows
  my(ok = 0, best = 0);
  forsubset([#cf, 4], ss, my(Ms = matrix(4, 4, r, c, M[ss[r], c]));
    if(matdet(Ms) != 0, my(l = matsolve(Ms, vector(4, r, rhs[ss[r]])~)); if(M*l == rhs, ok = 1; best = l; break)));
  if(ok, best, 0)};
{
  foreach([[5, 0], [4, 1]], c,
    my(m = c[1], lin = c[2], E = eqs(m, lin), S = mom(m, (u, v) -> v), T = mom(m, (u, v) -> u*v));
    \\ k = 5: claim S = 15/4 (resp 9/4) and T = -105/16 (resp -47/16) for the family, and S - T for the reflection
    print("k = 5, cubic members ", m, ", linear ", lin, ":");
    print("   E1 = ", E[1]); print("   E2 = ", E[2]); print("   E1r = ", E[3]); print("   E2r = ", E[4]));
}
