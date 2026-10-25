\\ order-M system in LIGHT-CONE shape coordinates: hp = H(1) = (1 - r1)(1 - r2), hm = H(-1) = (1 + r1)(1 + r2), bp = 1/hp, bm = 1/hm:
\\   e1 = (hm - hp)/2, e2 = (hp + hm)/2 - 1, Dl = hp hm, lam = bp - bm, mu = bp + bm - (bp - bm)(hm - hp)/2   (all polynomial).
\\ Relations hp bp = 1, hm bm = 1 (Dl != 0); Rabinowitsch y (hm - hp)(e1^2 - 4 e2) = 1 (off the slice e1 = 0, r1 != r2).
read("pc7_mono_hs.gp");
mkh(k, M, fname, ch, ord = "dp") = {
  my(Pv = vector(M, m, eval(Str("P", m))), E1 = (hm - hp)/2, E2 = (hp + hm)/2 - 1, LA = bp - bm, MU = bp + bm - (bp - bm)*(hm - hp)/2);
  my(E = hsym(k, E1, E2, LA, MU, M, Pv));
  E = vector(#E, i, my(q = E[i]); q / content(q));
  my(vs = strjoin(concat(vector(2*k - 1, a, Str("p", 2*k - a)), vector(M, m, Str("P", M + 1 - m))), ","));
  write1(fname, Str("ring R = ", ch, ", (y,", vs, ",bp,bm,hp,hm), ", ord, ";\noption(redSB);\nideal I =\n"));
  for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
  write1(fname, Str(hp*bp - 1, ",\n", hm*bm - 1, ",\n", y*(hm - hp)*(E1^2 - 4*E2) - 1, ";\n"));
  write1(fname, "int t0 = timer; ideal J = std(I); \"dim:\"; dim(J); \"size:\"; size(J); \"seconds:\"; timer - t0; if (size(J) < 10) { J; }\nquit;\n");
  print(fname, ": total length ", sum(i = 1, #E, #Str(E[i])));
};
mkh(7, 3, "pc7_mono_H3_p.sing", 32003);
mkh(7, 3, "pc7_mono_H3_0.sing", 0);
\\ control on a line: hm = hp + 4/15 (e1 = 2/15), char p -- should reproduce the unit ideal of L3 quickly
mkhl(k, M, fname, ch) = {
  my(Pv = vector(M, m, eval(Str("P", m))), HM = hp + 4/15, E1 = (HM - hp)/2, E2 = (hp + HM)/2 - 1, LA = bp - bm, MU = bp + bm - (bp - bm)*(HM - hp)/2);
  my(E = hsym(k, E1, E2, LA, MU, M, Pv)); E = vector(#E, i, my(q = E[i]); q / content(q));
  my(vs = strjoin(concat(vector(2*k - 1, a, Str("p", 2*k - a)), vector(M, m, Str("P", M + 1 - m))), ","));
  write1(fname, Str("ring R = ", ch, ", (y,", vs, ",bp,bm,hp), dp;\noption(redSB);\nideal I =\n"));
  for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
  write1(fname, Str(hp*bp - 1, ",\n", HM*bm - 1, ",\n", y*(E1^2 - 4*E2) - 1, ";\n"));
  write1(fname, "int t0 = timer; ideal J = std(I); \"dim:\"; dim(J); \"size:\"; size(J); \"seconds:\"; timer - t0; if (size(J) < 10) { J; }\nquit;\n");
};
mkhl(7, 3, "pc7_mono_HL3_p.sing", 32003);
