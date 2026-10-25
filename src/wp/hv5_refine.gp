\\ Refine a candidate solution of the fixed-V k = 5 "-" half (see hv5_newton.py) at high precision; report the
\\ residual history and the eigenvalues of J^H J (a near-zero eigenvalue => positive-dimensional component).
\\ Usage: gp -q hv5_refine.gp   (candidate V0 hard-coded below)
default(realprecision, 60);
K = 5;
Tpoly(v) = {my(b = vector(K, i, v[i]), d = vector(K, i, v[K+i]), a, G, e, S);
  a = vector(K, i, 2 + ('y+1)*d[i] + ('y+1)^2*b[i]);
  G = prod(i = 1, K, 1 + 'Xv*a[i]); e = vector(K+1, r, polcoef(G, r-1, 'Xv));
  S = sum(r = 1, K-1, 'y^(r-1)*e[r+1]) - e[K+1]*(15/4*'y + 45/16);
  S \ ('y+1)^2};
Fv(v) = {my(T = Tpoly(v)); vector(10, m, polcoef(T, m-1, 'y))~};
Jac(v) = {my(h = 10.^-25, f = Fv(v), cols = vector(10)); for(c = 1, 10, my(v2 = v); v2[c] += h; cols[c] = (Fv(v2) - f)/h); matconcat(cols)};
refine(v0) = {my(v = v0, f, J);
  for(it = 1, 30, f = Fv(v); J = Jac(v); v = v - (matsolve(J, f))~;
    if(it % 5 == 0, print("it ", it, "  |F| = ", normlp(Fv(v)) * 1., "  |v| = ", normlp(v) * 1.)));
  v};
V0 = [2.678983-0.800074*I, 0.887023-0.415276*I, 0.996391+0.914024*I, 0.54138-0.007592*I, 7.571175+5.102097*I, -4.880569+4.580591*I, -2.768024+0.698006*I, -3.388064-1.499503*I, -1.763176-0.017745*I, -6.669141-4.045422*I];
w = refine(V0);
print("b = ", vector(K, i, round(w[i] * 10^8) / 10.^8));
print("d = ", vector(K, i, round(w[K+i] * 10^8) / 10.^8));
J = Jac(w); M = conj(J)~ * J;
print("eigenvalues of J^H J: ", vecsort(apply(t -> abs(t) * 1., polroots(charpoly(M)))));
