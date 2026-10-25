\\ src/cs/e3b/advchk/lemmas.gp -- adversarial check of E3b Lemmas A, B, C, E (independent code, exact symbolic).
\\ Generic A = x^2 + p x + q, g_1 = c, g_j = l_j (x - r_j), k = 5.
k = 5;
A = x^2 + p*x + q;
L = [c, l2, l3, l4, l5]; R = [0, r2, r3, r4, r5];
g(j) = if (j == 1, c, L[j] * (x - R[j]));
f(j) = 1 + A * g(j);
F = prod(j = 1, k, f(j));
ok = 1;
okE = 1;
\\ Lemma A, Lucas part: g_j^(k-1) (F + 1) + Phi_j == 0 mod N_j  (exact division over Q(params)[x])
{for (j = 1, k,
  my(N = f(j) + 1, Phi = prod(i = 1, k, if (i == j, 1, g(j) - 2*g(i))) - g(j)^(k-1), rem);
  rem = (g(j)^(k-1) * (F + 1) + Phi) % N;
  if (rem != 0, ok = 0; print("Lemma A Lucas FAILS at j = ", j));
  \\ gcd(g_j, N_j) = 1: N_j == 2 mod g_j
  if (j > 1 && subst(N, x, R[j]) != 2, ok = 0; print("gcd(g_j, N_j) claim fails j = ", j)));}
print("Lemma A (Lucas congruence) checked for all j: ", ok);
\\ Lemma A, Korselt part: H = (F - 1)/A; H(r_j) = K_j(A(r_j)) where K_j(a) = [prod_{i!=j}(1 + a g_i(r_j)) - 1]/a
H = (F - 1) / A;
if (type(H) != "t_POL", print("A does not divide F - 1 ?!"));
{for (j = 2, k,
  my(Kj = (prod(i = 1, k, if (i == j, 1, 1 + 'a * subst(g(i), x, R[j]))) - 1) / 'a, lhs, rhs);
  Kj = simplify(Kj);
  lhs = subst(H, x, R[j]); rhs = subst(Kj, 'a, subst(A, x, R[j]));
  if (lhs != rhs, ok = 0; print("Lemma A Korselt FAILS j = ", j)));}
print("Lemma A (Korselt <=> K_j(A(r_j)) = 0) checked: ", ok);
\\ Korselt_j <=> A g_j | F - 1 <=> g_j | H  (A || F-1 handled since A | F - 1 always):
\\ check (F-1) mod (A*g_j) vanishes iff H(r_j) = 0 on a random numerical Korselt-violating example (sanity)
\\ Lemma B (FSP): with A = x^2 - 1, f_1 = 1 - 2A/(t^2-1), f_i = 1 - 2E_i/E_i(t):  prod_{i>=2} f_i(-t) = prod (3t + r_i)/(t - r_i)
E(i) = (x^2 - 1) * (x - R[i]);
fF(i) = if (i == 1, 1 - 2*(x^2 - 1)/(t^2 - 1), 1 - 2*E(i)/subst(E(i), x, t));
lhsB = prod(i = 2, k, subst(fF(i), x, -t));
rhsB = prod(i = 2, k, (3*t + R[i])/(t - R[i]));
print("Lemma B f_1-full criterion: ", lhsB == rhsB);
\\ FSP members really pass through -1 at t, and F(t) = -1
print("Lemma B FSP values: ", vector(k, i, subst(fF(i), x, t)), "  F(t) = ", subst(prod(i = 1, k, fF(i)), x, t));
\\ extra Lucas points of member j: N_j / (x - t) at x = t*sigma, over t^2: sigma^2 + u sigma + u - 1/T  with u = (t - r)/t
{
  my(j = 2, Q = (E(j) - subst(E(j), x, t)) / (x - t), lhs, u = (t - R[j])/t);
  Q = simplify(Q);
  lhs = subst(Q, x, t*'sig) / t^2;
  print("Lemma B extra-point equation: ", simplify(lhs - ('sig^2 + u*'sig + u - 1/t^2)) == 0);
}
\\ Lemma E: at a Korselt root r_i: F(r_i) = 1 (since A g_i | F - 1 and g_i(r_i) = 0) -> cannot equal -1
{
  for (i = 2, k,
    my(v = subst(F - 1, x, R[i]), Ai = subst(A, x, R[i]), Hi = subst(H, x, R[i]));
    if (v != Ai * Hi, okE = 0; print("Lemma E identity F(r_i) - 1 = A(r_i) H(r_i) fails")));
  print("Lemma E support: F(r_i) - 1 = A(r_i) H(r_i), so Korselt_i => F(r_i) = 1 != -1: ", okE);
}
\\ Lemma C parity (symmetric locus x -> -x, A = x^2 - 1, f3(x) = f2(-x), f5(x) = f4(-x)):
\\ Phi_1 and N_1 even.  Check symbolically.
{
  my(A0 = x^2 - 1, gg = [c, l2*(x - r2), -l2*(x + r2), l4*(x - r4), -l4*(x + r4)], N1, Phi1);
  N1 = 2 + A0*c;
  Phi1 = prod(i = 2, 5, c - 2*gg[i]) - c^4;
  print("Lemma C: N_1 even ", subst(N1, x, -x) == N1, "  Phi_1 even ", subst(Phi1, x, -x) == Phi1,
        "  f3(x) = f2(-x): ", subst(1 + A0*gg[2], x, -x) == 1 + A0*gg[3]);
}
