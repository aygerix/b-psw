\\ pc7_cert_shape.gp -- EXACT per-shape data for the pointwise functionals E1..E8 (rational shapes r1 <= r2 in (-1, 1)):
\\ (B-basis: lam3 = [1,0,0,0,-1,0,0,0]; E-basis: lam3 = [1,1,0,0,-1,-1,0,0])
\\ E_i = kap_i + sum_j a_i(w_j) + sum_{i<j} b_i(w_i, w_j) + lam3_i e3,  b_i(w, v) = sum_{n=2}^{6} c_{i,n} sum_{a+b=n, a,b>=1} w^a v^b.
\\ shdata(R1, R2) = [kap (8), A (8 x 9: coeff of w^1..w^9), C (8 x 5: c_2..c_6)]; dump(file, list of shapes) writes them.
read("pc7_cert_dec.gp");
r1; r2;
Eof(ws) = {my(B = Bof(ws));
  [B[1] + B[2]*r1, B[1] + B[2]*r2, B[3] + B[4]*r1, B[3] + B[4]*r2, B[5] + B[6], B[5] - B[6], B[7] + B[8], B[7] - B[8]]};
shsub(P, R1, R2) = {my(Dl = (1 - R1^2)*(1 - R2^2)); substvec(P, [r1, r2, e1, e2, lam, mu], [R1, R2, R1 + R2, R1*R2, 2*(R1 + R2)/Dl, 2*(1 - (R1 + R2)^2 + R1*R2)/Dl])};
\\ symbolic forms once
\\ D-basis: coefficients t^0..t^3 of (V - c0 G^6) mod H^2 and of (W - c0 H^6) mod G^2.  Equivalent to B1..B8 = 0 when r1 != r2,
\\ and the CORRECT (rank-4) form at a double root r1 = r2, where the B-encoding (Delta mod H, Delta' mod H) loses Delta'''(r).
Dof(ws) = {my(V = Vof(ws), W = Wof(ws), R1 = (V - c0*G^6) % H^2, R2 = (W - c0*H^6) % G^2);
  concat(vector(4, i, polcoef(R1, i - 1, t)), vector(4, i, polcoef(R2, i - 1, t)))};
BASIS = 2;  \\ 2: D-basis; 1: B-basis (B1..B8 of q22_moments.gp); 0: pointwise E-basis.  lam3 = [1,0,0,0,-1,0,0,0] for B and D
mkSYM() = if(BASIS == 2, SYM0 = Dof([]); SYM1 = Dof([w]); SYM2 = Dof([w, v]), if(BASIS, SYM0 = Bof([]); SYM1 = Bof([w]); SYM2 = Bof([w, v]), SYM0 = Eof([]); SYM1 = Eof([w]); SYM2 = Eof([w, v])));
mkSYM();
shdata(R1, R2) = {my(E0 = shsub(SYM0, R1, R2), E1 = shsub(SYM1, R1, R2), E2 = shsub(SYM2, R1, R2));
  my(A = matrix(8, 9, i, m, polcoef(E1[i] - E0[i], m, w)));
  my(C = matrix(8, 5, i, n, polcoef(polcoef(E2[i] - E0[i] - (E1[i] - E0[i]) - subst(E1[i] - E0[i], w, v), 1, w), n, v)));
  \\ verify Hankel structure exactly
  for(i = 1, 8, my(Bp = E2[i] - E0[i] - (E1[i] - E0[i]) - subst(E1[i] - E0[i], w, v));
    my(T = sum(n = 2, 6, C[i, n - 1] * sum(a = 1, n - 1, w^a * v^(n - a))));
    if(Bp != T, error("Hankel structure fails at ", [R1, R2, i])));
  [E0, A, C]};
dump(f, L) = {for(k = 1, #L, my(S = shdata(L[k][1], L[k][2]));
  write(f, "S ", L[k][1], " ", L[k][2]);
  write(f, "K ", S[1]);
  for(i = 1, 8, write(f, "A ", S[2][i, ]));
  for(i = 1, 8, write(f, "C ", S[3][i, ])))};
