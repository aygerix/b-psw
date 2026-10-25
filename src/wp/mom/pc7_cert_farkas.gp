\\ exact check of the obstruction: D-vectors at the rational configurations CONF (rows) at SHAPE; find theta >= 0, not all 0,
\\ with sum theta_k D(w^(k)) = 0.  If the D's span R^8 and such theta > 0 exists, no constant alpha separates (Farkas).
read("pc7_cert_shape.gp");
read("pc7_cert_obs.gp");
{
  my(R1 = SHAPE[1], R2 = SHAPE[2], m = matsize(CONF)[1]);
  my(Dn = shsub(SYM0, R1, R2), D1 = shsub(SYM1, R1, R2), D2 = shsub(SYM2, R1, R2));
  my(A = vector(8, i, D1[i] - Dn[i]), Bp = vector(8, i, D2[i] - Dn[i] - A[i] - subst(A[i], w, v)), lam3 = [1,0,0,0,-1,0,0,0]);
  my(M = matrix(8, m));
  for(k = 1, m, my(ws = CONF[k, ]);
    for(i = 1, 8, M[i, k] = Dn[i] + sum(j = 1, 7, subst(A[i], w, ws[j])) + sum(j = 1, 7, sum(l = j + 1, 7, substvec(Bp[i], [w, v], [ws[j], ws[l]])))
       + lam3[i] * sum(j = 1, 7, sum(l = j + 1, 7, sum(q = l + 1, 7, ws[j]*ws[l]*ws[q])))));
  print("configurations ", m, ", rank of D-matrix ", matrank(M));
  \\ normalise columns (positive scaling does not change the cone) and look at the kernel
  my(Mn = matrix(8, m, i, k, M[i, k] / vecsum(abs(M[, k])~)));
  my(Kr = matker(Mn));
  print("kernel dimension ", matsize(Kr)[2]);
  if(matsize(Kr)[2] == 1, my(th = Kr[, 1]); if(vecmax(th~) <= 0, th = -th);
     print("kernel vector signs: ", apply(sign, th~)); print("all positive: ", vecmin(th~) > 0));
  \\ greedy: try subsets of 9 columns of largest weight
  my(sel = vector(min(9, m), k, k), Ms = matrix(8, #sel, i, k, Mn[i, sel[k]]), Ks = matker(Ms));
  if(matsize(Ks)[2] == 1, my(th = Ks[, 1]); if(vecmax(th~) <= 0, th = -th);
     print("first 9 configs: kernel signs ", apply(sign, th~), "  rank ", matrank(Ms));
     print("EXACT FARKAS WITNESS (no constant alpha at this shape): ", vecmin(th~) > 0 && matrank(Ms) == 8));
}
