\\ structure: rho-(c) = (X + Y)/(Dl^4 c^4), rho+(c) = (X - Y)/(Dl^4 c^4), X = Dl^4 c^3, Y = A2 c^2 + A1 c + A0.
\\ In w = 1/c: phi(w) = Y/X = (A2 w + A1 w^2 + A0 w^3)/Dl^4. Express A0, A1, A2 in e1, e2 and lam, mu.
default(parisizemax, 1000000000);
t; c; r1; r2; e1; e2;
reA(A, Q) = {my(R = A % Q, a = polcoef(Q, 2, t), b = polcoef(Q, 1, t)); polcoef(R, 1, t) * (-b/(2*a)) + polcoef(R, 0, t)};
sym(P) = {  \\ symmetric polynomial in r1, r2 -> polynomial in e1 = r1 + r2, e2 = r1 r2
  my(Q = P, out = 0);
  while(Q != 0, my(d1 = poldegree(Q, r1), cf = polcoef(Q, d1, r1), d2 = poldegree(cf, r2), lc = polcoef(cf, d2, r2));
    \\ leading monomial r1^d1 r2^d2 with d1 >= d2: e1^(d1-d2) e2^d2
    if(d1 < d2, error("not symmetric"));
    out += lc * e1^(d1 - d2) * e2^d2; Q -= lc * (r1 + r2)^(d1 - d2) * (r1*r2)^d2; Q = simplify(Q));
  out};
{
  my(E1 = r1 + r2, E2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2), lamN = 2*E1, muN = 2*(1 - E1^2 + E2));
  my(G = t^2 - 1, H = t^2 - E1*t + E2, LN = lamN*t + muN, LpN = lamN*t + muN + lamN*E1);
  my(Rm = reA(G*H, c*Dl*H + LN), Rp = reA(G*H, c*Dl*G + LpN));
  my(Nm = numerator(Rm), Np = numerator(Rp), X = polcoef(Nm, 3, c) * c^3, Y = Nm - X);
  print("X = Dl^4 c^3: ", X == Dl^4 * c^3, "   Np = X - Y: ", Np == X - Y, "   denominators Dl^4 c^4: ", denominator(Rm) == Dl^4*c^4, " ", denominator(Rp) == Dl^4*c^4);
  my(A2 = sym(polcoef(Y, 2, c)), A1 = sym(polcoef(Y, 1, c)), A0 = sym(polcoef(Y, 0, c)), DL = sym(Dl));
  print("Dl = ", DL);
  print("A0 = ", factor(A0)); print("A1 = ", factor(A1)); print("A2 = ", factor(A2));
}
