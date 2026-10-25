\\ Check of the A-reduction lemma  gcd(F+1, f_j+1) = gcd(f_j+1, L_j),  L_j = (c_j B_j)^(k-1) - prod_{i!=j}(c_j B_j - 2 c_i B_i),
\\ on (i) the rho-family (k = 3), (ii) Lemma 5 (k = 3, shared A = (Y-2)(Y+1)), (iii) random families f_j = 1 + c_j A B_j, k = 3, 5.
Lred(c, B, j) = (c[j]*B[j])^(#c - 1) - prod(i = 1, #c, if (i == j, 1, c[j]*B[j] - 2*c[i]*B[i]));
check(A, c, B, tag) = {
  my(k = #c, f = vector(k, j, 1 + c[j]*A*B[j]), F = prod(j = 1, k, f[j]), ok = 1, out = []);
  for (j = 1, k, my(g1 = gcd(F + 1, f[j] + 1), L = Lred(c, B, j), g2 = if (L == 0, f[j] + 1, gcd(f[j] + 1, L)));
    ok = ok && (poldegree(g1) == poldegree(g2)) && (g1/pollead(g1) == g2/pollead(g2)); out = concat(out, poldegree(g1)));
  print(tag, ": gcd degrees ", out, "   A-reduction identity holds: ", ok) };
{
  my(r = 3/7);                                                  \\ rho-family at rho = 3/7
  check(2*x^2 - (4*r + 1)*x + 2*r - 2, [1, 1/(2*r), 1/(2*r - 1)], [1, x, x - 1], "rho-family rho=3/7");
  check((x - 2)*(x + 1), [1, 1, 1], [1, x, x + 1], "Lemma 5 (Psi10, Psi14, Psi18)");   \\ Psi14 - 1 = Y(Y-2)(Y+1), Psi18 - 1 = (Y-2)(Y+1)^2
  setrand(1);
  for (t = 1, 6, my(k = if (t <= 3, 3, 5), a = random(3) + 1, A = x^a + sum(i = 0, a - 1, (random(7) - 3)*x^i));
    my(c = vector(k, j, random(9) - 4 + 1/(j+1)), B = vector(k, j, my(b = random(3)); x^b + sum(i = 0, b - 1, (random(7) - 3)*x^i)));
    check(A, c, B, Str("random k=", k, " a=", a)));
}
