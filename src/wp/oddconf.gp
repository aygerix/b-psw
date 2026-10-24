\\ Linear Williams families in "odd configuration" (Lemma 22): F(t) = t^m0 * prod_{i<=r} (t^2 - e_i^2) / U, k = m0 + 2r odd.
\\ A linear family f_j = (t - s_j)/d_j exists iff for every root s_j: F(s_j + d_j) = 1, F(s_j - d_j) = -1, and prod d_j = U.
\\ Root 0 (symmetry centre): any d with F(d) = 1 works (k candidates D). Root +-e_i: d in X_i = {x : F(e_i+x) = 1, F(e_i-x) = -1}
\\ (the same set serves -e_i, by oddness). Product condition: prod(chosen m0 of the D's) * prod_i x_i x_i' = U.
\\ Since prod(all k D's) = U, this says: prod of the 2r OMITTED D's = prod_i x_i x_i'.
\\ This script: r = 1 exactly (all k <= kmax), reporting the smallest |log| mismatch over all choices.
default(realprecision, 60);
r1(k) = {
  my(m0 = k - 2, best = [1e100, 0, 0], g, xs, cnt = 0, hits = 0);
  g = (1 + 'x)^m0 * (2 + 'x);
  xs = polroots(substpol((g - subst(g, 'x, -'x)), 'x, 'x) / 'x);   \\ odd part / x: nonzero solutions
  for (a = 1, #xs,
    my(x = xs[a], U = x * subst(g, 'x, x));
    if (abs(x) < 1e-30, next);
    \\ valid set X for this U (other x' with g(x') = g(-x') and x' g(x') = U)
    my(X = []);
    for (b = 1, #xs, my(y = xs[b]); if (abs(y * subst(g, 'x, y) - U) < 1e-40 * (1 + abs(U)), X = concat(X, [y])));
    my(D = polroots('d^k - 'd^m0 - U));
    for (i = 1, k, for (j = i + 1, k,
      my(pr = D[i] * D[j]);
      for (u = 1, #X, for (v = u, #X,
        cnt++;
        my(err = abs(pr - X[u] * X[v]) / (1 + abs(pr)));
        if (err < best[1], best = [err, x, pr]);
        if (err < 1e-40, hits++; print("  HIT k=", k, " x=", x, " U=", U, " pair ", D[i], " ", D[j])))))));
  [cnt, hits, best[1]];
}
kmax = 41;
forstep (k = 5, kmax, 2, my(v = r1(k)); printf("k = %d (m0 = %d, r = 1): %d (x, pair) choices, %d exact hits, smallest relative mismatch %.3e\n", k, k - 2, v[1], v[2], v[3]));
