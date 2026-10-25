\\ AUDIT (numeric, HEURISTIC evidence): multi-start complex Newton for the two-dead system at k = 5, 7, 9, fixed M:
\\ unknowns a_2..a_{k-1}; equations prod_{i != j}(a_j - 2a_i) = a_j^(k-1) (polynomial form), j = 2..k-1.
\\ Collect distinct unordered solutions with pairwise distinct entries, all != 0, 1, M. Compare with the roots of G_k.
default(realprecision, 60);
eqs(s, k) = vector(k-2, r, my(j = r+1); prod(i=1, k, if(i==j, 1, s[j] - 2*s[i])) - s[j]^(k-1))~;
jac(s, k) = {matrix(k-2, k-2, r, c, my(j = r+1, cc = c+1);
   if (cc != j, -2*prod(i=1, k, if(i==j || i==cc, 1, s[j] - 2*s[i])),
      sum(l=1, k, if(l==j, 0, prod(i=1, k, if(i==j || i==l, 1, s[j] - 2*s[i])))) - (k-1)*s[j]^(k-2)))};
run(k, M, N) = {
  my(sols = List(), nconv = 0);
  for(t = 1, N,
    my(x = vector(k-2, j, (random(1.)*2 - 1)*M*1.2 + I*(random(1.)*2-1)*M*0.6), s, ok = 0);
    for(it = 1, 150, s = concat([1], concat(x, [M]));
      my(F = eqs(s, k)); if (normlp(F) < 1e-35 * max(1, normlp(x))^(k-1), ok = 1; break);
      my(J = jac(s, k)); if (abs(matdet(J)) < 1e-200, break);
      x = x - (matsolve(J, F))~);
    if (!ok, next);
    s = concat([1], concat(x, [M]));
    \\ non-degeneracy
    my(bad = 0);
    for(i = 1, k, for(j = i+1, k, if (abs(s[i]-s[j]) < 1e-12*(1+abs(s[i])), bad = 1)));
    for(i = 2, k-1, if (abs(s[i]) < 1e-6, bad = 1));
    if (bad, next);
    nconv++;
    my(key = vecsort(apply(z -> [round(real(z)*1e8), round(imag(z)*1e8)], x)));
    my(found = 0); for(u = 1, #sols, if (sols[u] == key, found = 1; break));
    if (!found, listput(sols, key)));
  [nconv, sols];
}
recP(k, b0, b1, q1) = {
  my(w = vector(k, i, 2^(k-2-(i-1))), q0 = 2^(3-k)*b0, p = vector(k+1));
  p[k-1] = 1; p[k] = 0; p[k+1] = 0;
  forstep (n = k-1, 2, -1,
    my(t = p[n] * (w[n] * q1 - b1) + if (n <= k-2, p[n+1] * (w[n+1] * q0 - b0), 0));
    p[n-1] = -t / (w[n-1] - 1));
  [Pol(vector(k-1, i, p[k-i]), 'x), 0] };
{
  setrand(7);
  foreach([[5, 7, 600], [5, 135, 600], [7, 7, 1500], [7, 1024, 1500], [9, 13, 1500]], c,
    my(k = c[1], M = c[2], r = run(k, M, c[3]));
    my(b1 = 1 + 1/'M, b0 = 1/'M, G = numerator(subst(recP(k, b0, b1, b1/2^(k-2))[1], 'x, -1/'y) * 'y^(k-2)));
    my(Gn = subst(G, 'M, M), gr = vecsort(apply(z -> [round(real(z)*1e8), round(imag(z)*1e8)], Vec(polroots(Gn)~))));
    my(close(u) = vecmax(vector(#gr, i, abs(u[i][1]-gr[i][1]) + abs(u[i][2]-gr[i][2]))) < 1000);
    printf("k=%d M=%d: converged non-degenerate %d, distinct unordered solutions %d; G_k root set among them: %d\n",
           k, M, r[1], #r[2], sum(u=1, #r[2], close(r[2][u])));
    print("   G_k roots: ", apply(v -> v/1e8, gr));
    for(u=1, #r[2], if (!close(r[2][u]), print("   OTHER: ", apply(v -> v/1e8, r[2][u])))));
}
