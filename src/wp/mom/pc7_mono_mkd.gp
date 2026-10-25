\\ d-FORMULATION of the order-M (2,2) pencil system, k = 7: unknowns d0..d8 with V = G^6 + sum d_i t^i (no power sums),
\\ W = phiVW(V) (connection map, pc7_mono_phi.gp, needs il = 1/lam), E1..EM (elementary symmetric functions of the c_j).
\\ Equations: bilinear V-side and W-side order-M conditions (pc7_mono_bil.gp).
\\ Shape: light-cone coordinates hp = H(1), hm = H(-1), bp = 1/hp, bm = 1/hm; lam = bp - bm, il lam = 1 (off the slice),
\\ y (e1^2 - 4 e2) = 1 (r1 != r2).
default(parisizemax, 2000000000);
read("pc7_mono_phi.gp"); read("pc7_mono_bil.gp");
dsys(k, M, e1, e2, lam, mu, il) = {
  my(dv = vector(2*k - 5, i, eval(Str("d", i - 1))), V = (t^2 - 1)^(k - 1) + sum(i = 1, 2*k - 5, dv[i] * t^(i - 1)));
  my(W = phiVW(k, V, e1, e2, lam, mu, il), Ev = vector(M, m, eval(Str("E", m))));
  bil(k, e1, e2, lam, mu, M, Ev, [V, W])};
mkd(k, M, fname, ch) = {
  my(t0 = getabstime(), E = dsys(k, M, e1, e2, lam, mu, il));
  my(S1 = (hm - hp)/2, S2 = (hp + hm)/2 - 1, SL = bp - bm, SM = bp + bm - (bp - bm)*(hm - hp)/2);
  E = vector(#E, i, subst(subst(subst(subst(E[i], e1, S1), e2, S2), lam, SL), mu, SM));
  my(vs = strjoin(concat(vector(2*k - 5, i, Str("d", 2*k - 5 - i)), vector(M, m, Str("E", M + 1 - m))), ","));
  write1(fname, Str("ring R = ", ch, ", (y,il,", vs, ",bp,bm,hp,hm), dp;\noption(redSB);\nideal I =\n"));
  for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
  write1(fname, Str(hp*bp - 1, ",\n", hm*bm - 1, ",\n", il*(bp - bm) - 1, ",\n", y*(S1^2 - 4*S2) - 1, ";\n"));
  write1(fname, "int t0 = timer; ideal J = std(I); \"dim:\"; dim(J); \"size:\"; size(J); \"seconds:\"; timer - t0; if (size(J) < 10) { J; }\nquit;\n");
  print(fname, ": ", #E, " eqs, total length ", sum(i = 1, #E, #Str(E[i])), " (", getabstime() - t0, " ms)");
};
\\ fixed-shape control
mkdf(k, M, r1, r2, fname, ch) = {
  my(E1 = r1 + r2, E2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2), LA = 2*E1/Dl, E = dsys(k, M, E1, E2, LA, 2*(1 - E1^2 + E2)/Dl, 1/LA));
  my(vs = strjoin(concat(vector(2*k - 5, i, Str("d", 2*k - 5 - i)), vector(M, m, Str("E", M + 1 - m))), ","));
  write1(fname, Str("ring R = ", ch, ", (", vs, "), dp;\noption(redSB);\nideal I =\n"));
  for(i = 1, #E, write1(fname, Str(E[i], if(i < #E, ",\n", ";\n"))));
  write1(fname, "ideal J = std(I); \"dim:\"; dim(J); \"size:\"; size(J); if (size(J) < 10) { J; }\nquit;\n");
};
mkdf(7, 1, 1/3, -1/5, "pc7_mono_DF1.sing", 32003);
mkdf(7, 2, 1/3, -1/5, "pc7_mono_DF2.sing", 32003);
mkdf(7, 3, 1/3, -1/5, "pc7_mono_DF3.sing", 32003);
mkd(7, 3, "pc7_mono_D3_p.sing", 32003);
mkd(7, 3, "pc7_mono_D3_0.sing", 0);
mkd(7, 2, "pc7_mono_D2_p.sing", 32003);
