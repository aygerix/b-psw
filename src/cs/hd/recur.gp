\\ Two-dead classification check (any k): dead members d = 1/a in {dm, dm'} -> B(x) = x^2 + b1 x + b0; good members are the roots
\\ of P (monic, deg k-2) given by the recurrence  p_{n-2}(w_{n-2}-1) + p_{n-1}(w_{n-1} q1 - b1) + p_n(w_n q0 - b0) = 0  (n = 2..k-1),
\\ w_i = 2^(k-2-i), q0 = 2^(3-k) b0, plus the constraint at n = 1.  Compare with the k = 5 cubic (dead slopes 1 and M).
recP(k, b0, b1, q1) = {
  my(w = vector(k, i, 2^(k-2-(i-1))), q0 = 2^(3-k)*b0, p = vector(k+1));   \\ p[i+1] = p_i, i = 0..k
  p[k-1] = 1; p[k] = 0; p[k+1] = 0;  \\ p_{k-2} = 1, p_{k-1} = p_k = 0
  forstep (n = k-1, 2, -1,
    \\ p_{n-2}(w_{n-2} - 1) = -[p_{n-1}(w_{n-1} q1 - b1) + p_n(w_n q0 - b0)]
    my(t = p[n] * (w[n] * q1 - b1) + if (n <= k-2, p[n+1] * (w[n+1] * q0 - b0), 0));
    p[n-1] = -t / (w[n-1] - 1));
  my(con = p[1] * (w[1] * q1 - b1) + p[2] * (w[2] * q0 - b0));
  [Pol(vector(k-1, i, p[k-i])), con] };
{
  my(k = 5, M = 135, b1 = 1 + 1/M, b0 = 1/M);
  my(r = recP(k, b0, b1, 'q), P = r[1], con = r[2]);
  print("constraint in q1 (numerator) factors: ", factor(numerator(con)));
  my(qs = polroots(numerator(con)));
  for (i = 1, #qs, if (abs(imag(qs[i])) > 1e-20, next);
    my(Pq = subst(P, 'q, real(qs[i])), d = polroots(Pq));
    print("q1 = ", real(qs[i]), " -> good slopes a = 1/d: ", vector(#d, j, 1/d[j])));
  print("cubic roots from c5exact: ", polroots(subst((7*M^3 + 57*M^2 + 57*M + 7)*y^3 + (-98*M^3 - 308*M^2 - 98*M)*y^2 + (392*M^3 + 392*M^2)*y - 448*M^3, 'M, M)));
}
