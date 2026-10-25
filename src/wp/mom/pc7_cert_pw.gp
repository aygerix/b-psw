\\ pc7_cert_pw.gp -- pointwise boundary functionals E1..E8 in 3-body form, written for Python:
\\ E = [dV(r1), dV(r2), dV'(r1), dV'(r2), dW(1), dW(-1), dW'(1), dW'(-1)], dV = V - c0 G^6, dW = W - c0 H^6.
\\ E_i(w_1..w_7) = kap_i + sum_j a_i(w_j) + sum_{i<j} b_i(w_i, w_j) + lam3_i e3(w), lam3 = [1,1,0,0,-1,-1,0,0].
\\ Output kept POLYNOMIAL in (r1, r2, e1, e2, lam, mu); the user substitutes e1 = r1 + r2, e2 = r1 r2, lam, mu.
read("pc7_cert_dec.gp");
r1; r2;
Eof(ws) = {my(B = Bof(ws));
  [B[1] + B[2]*r1, B[1] + B[2]*r2, B[3] + B[4]*r1, B[3] + B[4]*r2, B[5] + B[6], B[5] - B[6], B[7] + B[8], B[7] - B[8]]};
{
  my(E0 = Eof([]), E1 = Eof([w]), E2 = Eof([w, v]), f = "pc7_cert_pw.txt");
  my(A = vector(8, i, E1[i] - E0[i]));
  my(Bp = vector(8, i, E2[i] - E0[i] - A[i] - subst(A[i], w, v)));
  my(E3 = Eof([w, v, u]), lam3 = [1,1,0,0,-1,-1,0,0]);
  my(chk = vector(8, i, E3[i] - (E0[i] + A[i] + subst(A[i], w, v) + subst(A[i], w, u) + Bp[i] + subst(Bp[i], v, u) + subst(subst(Bp[i], v, u), w, v) + lam3[i]*w*v*u)));
  print("3-body check (must be zeros): ", chk);
  write(f, "kap = ", E0);
  for(i = 1, 8, write(f, "a", i, " = ", Vec(A[i] / w) ));
  for(i = 1, 8, write(f, "b", i, " = ", Bp[i]));
  print("degrees in w: ", vector(8, i, poldegree(A[i], w)));
}
