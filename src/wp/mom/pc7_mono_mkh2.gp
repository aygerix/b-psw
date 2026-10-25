\\ light-cone shape coordinates via a Singular ring map (fast): the ideal is generated in (e1, e2, lam, mu) and mapped by
\\   e1 -> (hm - hp)/2, e2 -> (hp + hm)/2 - 1, lam -> bp - bm, mu -> bp + bm - (bp - bm)(hm - hp)/2;  relations hp bp = 1, hm bm = 1;
\\   Rabinowitsch y (hm - hp)(e1^2 - 4 e2) = 1.
read("pc7_mono_hs.gp");
mkh(k, M, fname, ch) = {
  my(Pv = vector(M, m, eval(Str("P", m))), E = hsym(k, e1, e2, lam, mu, M, Pv));
  my(vs = strjoin(concat(vector(2*k - 1, a, Str("p", 2*k - a)), vector(M, m, Str("P", M + 1 - m))), ","));
  write1(fname, Str("ring R0 = ", ch, ", (", vs, ",lam,mu,e1,e2), dp;\nideal I0 =\n"));
  for(i = 1, #E, write1(fname, Str(E[i], if(i < #E, ",\n", ";\n"))));
  write1(fname, Str("ring R = ", ch, ", (y,", vs, ",bp,bm,hp,hm), dp;\noption(redSB);\n"));
  write1(fname, Str("map phi = R0, ", vs, ", bp - bm, bp + bm - (bp - bm)*(hm - hp)/2, (hm - hp)/2, (hp + hm)/2 - 1;\n"));
  write1(fname, "ideal I = phi(I0), hp*bp - 1, hm*bm - 1, y*(hm - hp)*((hm - hp)^2/4 - 2*(hp + hm) + 4) - 1;\n");
  write1(fname, "int t0 = timer; ideal J = std(I); \"dim:\"; dim(J); \"size:\"; size(J); \"seconds:\"; timer - t0; if (size(J) < 10) { J; }\nquit;\n");
};
mkh(7, 3, "pc7_mono_H3_p.sing", 32003);
mkh(7, 3, "pc7_mono_H3_0.sing", 0);
mkh(7, 4, "pc7_mono_H4_p.sing", 32003);
