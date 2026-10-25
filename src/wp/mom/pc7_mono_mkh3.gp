\\ light-cone shape coordinates, substitution done in gp on the final (e1, e2, lam, mu) equations:
\\   e1 -> (hm - hp)/2, e2 -> (hp + hm)/2 - 1, lam -> bp - bm, mu -> bp + bm - (bp - bm)(hm - hp)/2;
\\   relations hp bp = 1, hm bm = 1 (Dl = hp hm != 0); Rabinowitsch y (hm - hp)(e1^2 - 4 e2) = 1.
default(parisizemax, 2000000000);
read("pc7_mono_hs.gp");
mkh(k, M, fname, ch) = {
  my(Pv = vector(M, m, eval(Str("P", m))), E = hsym(k, e1, e2, lam, mu, M, Pv), t0 = getabstime());
  my(S1 = (hm - hp)/2, S2 = (hp + hm)/2 - 1, SL = bp - bm, SM = bp + bm - (bp - bm)*(hm - hp)/2);
  E = vector(#E, i, my(q = subst(subst(subst(subst(E[i], e1, S1), e2, S2), lam, SL), mu, SM)); q);
  my(vs = strjoin(concat(vector(2*k - 1, a, Str("p", 2*k - a)), vector(M, m, Str("P", M + 1 - m))), ","));
  write1(fname, Str("ring R = ", ch, ", (y,", vs, ",bp,bm,hp,hm), dp;\noption(redSB);\nideal I =\n"));
  for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
  write1(fname, Str(hp*bp - 1, ",\n", hm*bm - 1, ",\n", y*(hm - hp)*(S1^2 - 4*S2) - 1, ";\n"));
  write1(fname, "int t0 = timer; ideal J = std(I); \"dim:\"; dim(J); \"size:\"; size(J); \"seconds:\"; timer - t0; if (size(J) < 10) { J; }\nquit;\n");
  print(fname, ": total length ", sum(i = 1, #E, #Str(E[i])), " (", getabstime() - t0, " ms)");
};
mkh(7, 3, "pc7_mono_H3_p.sing", 32003);
mkh(7, 3, "pc7_mono_H3_0.sing", 0);
