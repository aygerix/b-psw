\\ m9 line: explicit MONOMIAL-moment form of the four boundary equations, k = 9, 11.
\\ Moment symbols M(a,b) = sum_i v_i^a w_i^b, w_i = u_i - 1/2 (reflection u -> 1-u is w -> -w, M(a,b) -> (-1)^b M(a,b)).
\\ Equations:  Z0 := T(-1) - 2^(1-k) = 0,  Z1 := T'(-1) = 0  (family), and the same with M(a,b) -> (-1)^b M(a,b).
\\ Even/odd split: Z^+ = (Z + Z^r)/2 (even in the odd moments), Z^- = (Z - Z^r)/2.
read("m9_moments.gp");
msym(a, b) = eval(Str("M", a, "_", b));
\\ ell_n as polynomial in (w, v) -> linear form in moment symbols
lin2mom(P) = {my(R = 0, Q = P);   \\ P polynomial in w, v with v | every monomial
  for(a = 1, poldegree(Q, v), my(C = polcoef(Q, a, v)); for(b = 0, poldegree(C, w), R += polcoef(C, b, w) * msym(a, b)));
  if(polcoef(Q, 0, v) != 0, error("ell has a v-free term")); R};
mform(m, lin) = {
  my(k = m + lin, c0 = 2^(1-k), Ls = Lsym(k), TF = Tform(k, Ls), el = ellv(1/2 + w, v, k), Lm = vector(k - 1), Z);
  for(n = 1, k - 1, Lm[n] = lin * ellin(k)[n] + if(n >= 3, lin2mom(el[n]), 0));
  Z = [TF[1] - c0, TF[2]];
  for(j = 1, 2, for(n = 1, k - 1, Z[j] = subst(Z[j], Ls[n], Lm[n])));
  [Z, Lm]};
\\ symbols occurring in a polynomial
syms(P) = select(z -> Vec(Str(z))[1] == "M", variables(P));
reflect(P) = {my(R = P); foreach(syms(P), z, my(st = Str(z), b = eval(strsplit(st, "_")[2])); if(b % 2, R = subst(R, z, -z))); R};
evalmom(P, m) = {my(R = P); foreach(syms(P), z, my(st = strsplit(Str(z), "_"), a = eval(strjoin(Vec(st[1])[2..#Vec(st[1])])), b = eval(st[2]));
  R = subst(R, z, sum(i = 1, m, eval(Str("v", i))^a * (eval(Str("u", i)) - 1/2)^b))); R};
