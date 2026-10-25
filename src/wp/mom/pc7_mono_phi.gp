\\ pc7_mono_phi: the CONNECTION MAP V -> W without power sums. With s = 1/t, l(s) = L/H, sg(s) = 1/(G H) (both series in s),
\\ A(s) := Lambda'(l(s)), Lambda(z) = sum_j log(1 + w_j z):
\\   S := log(s^12 V(1/s) / (1 - s^2)^6) = sg A + sg^2 Lambda''/2 + sg^3 Lambda'''/6      (mod s^13)
\\   T := log(s^12 W(1/s) / (1 - e1 s + e2 s^2)^6) = -[sg A + (3/2) sg^2 Lambda'' + (7/6) sg^3 Lambda''']   (mod s^13)
\\ with Lambda''(l(s)) = A'(s)/l'(s), Lambda'''(l(s)) = (A'/l')'/l'.  Off the slice (lam != 0) A is recovered from S by iteration,
\\ so W is an explicit polynomial function of the coefficients of V (cubic), with 1/lam = il.
\\ General k odd: s^(2k-2) and exponent k - 1; N = 2k - 1 series terms.
phiVW(k, V, e1, e2, lam, mu, il) = {
  my(N = 2*k - 1, so = O(s^N), dg = 2*k - 2);
  my(Vs = sum(i = 0, dg, polcoef(V, i, t) * s^(dg - i)) + so);
  my(l = s*(lam + mu*s)/(1 - e1*s + e2*s^2) + so, sg = s^4/((1 - s^2)*(1 - e1*s + e2*s^2)) + so);
  \\ 1/l'(s) with lam^-1 = il: l' = lam + O(s); write l' = lam (1 + u), 1/l' = il (1 - u + u^2 - ...), u = (l' - lam) il
  my(lp = deriv(l, s), u = (lp - lam)*il, ilp = il * sum(m = 0, N, (-u)^m) + so);
  my(S = log(Vs / (1 - s^2)^(k - 1)));
  my(D1(X) = deriv(X, s) * ilp);   \\ d/dl
  my(A = S / sg);                   \\ first approximation, then two corrections
  for(it = 1, 3, A = (S - sg^2 * D1(A)/2 - sg^3 * D1(D1(A))/6) / sg);
  my(T = -(sg * A + 3/2 * sg^2 * D1(A) + 7/6 * sg^3 * D1(D1(A))));
  my(Ws = (1 - e1*s + e2*s^2)^(k - 1) * exp(T));
  sum(i = 0, dg, polcoef(Ws, dg - i, s) * t^i)};
