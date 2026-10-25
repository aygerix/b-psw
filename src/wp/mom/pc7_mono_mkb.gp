\\ BILINEAR order-M system (pc7_mono_bil.gp: linear in E_m = e_m(c)), light-cone shape coordinates (see pc7_mono_mkh3.gp).
default(parisizemax, 2000000000);
read("pc7_mono_bil.gp");
mkb(k, M, fname, ch, lc = 1) = {
  my(Ev = vector(M, m, eval(Str("E", m))), E = bil(k, e1, e2, lam, mu, M, Ev), t0 = getabstime());
  my(S1 = (hm - hp)/2, S2 = (hp + hm)/2 - 1, SL = bp - bm, SM = bp + bm - (bp - bm)*(hm - hp)/2);
  my(vs = strjoin(concat(vector(2*k - 1, a, Str("p", 2*k - a)), vector(M, m, Str("E", M + 1 - m))), ","));
  if(lc,
    E = vector(#E, i, subst(subst(subst(subst(E[i], e1, S1), e2, S2), lam, SL), mu, SM));
    write1(fname, Str("ring R = ", ch, ", (y,", vs, ",bp,bm,hp,hm), dp;\noption(redSB);\nideal I =\n"));
    for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
    write1(fname, Str(hp*bp - 1, ",\n", hm*bm - 1, ",\n", y*(hm - hp)*(S1^2 - 4*S2) - 1, ";\n")),
    my(Dl = 1 - e1^2 + 2*e2 + e2^2);
    write1(fname, Str("ring R = ", ch, ", (y,", vs, ",lam,mu,e1,e2), dp;\noption(redSB);\nideal I =\n"));
    for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
    write1(fname, Str(lam*Dl - 2*e1, ",\n", mu*Dl - 2*(1 - e1^2 + e2), ",\n", y*e1*Dl*(e1^2 - 4*e2) - 1, ";\n")));
  write1(fname, "int t0 = timer; ideal J = std(I); \"dim:\"; dim(J); \"size:\"; size(J); \"seconds:\"; timer - t0; if (size(J) < 10) { J; }\nquit;\n");
  print(fname, ": total length ", sum(i = 1, #E, #Str(E[i])), " (", getabstime() - t0, " ms)");
};
mkb(7, 3, "pc7_mono_B3h_p.sing", 32003, 1);
mkb(7, 3, "pc7_mono_B3e_p.sing", 32003, 0);
mkb(7, 3, "pc7_mono_B3h_0.sing", 0, 1);
\\ fixed-shape control (bilinear form must also be unit at a fixed generic shape)
mkbf(k, M, r1, r2, fname, ch) = {
  my(E1 = r1 + r2, E2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2), Ev = vector(M, m, eval(Str("E", m))));
  my(E = bil(k, E1, E2, 2*E1/Dl, 2*(1 - E1^2 + E2)/Dl, M, Ev));
  my(vs = strjoin(concat(vector(2*k - 1, a, Str("p", 2*k - a)), vector(M, m, Str("E", M + 1 - m))), ","));
  write1(fname, Str("ring R = ", ch, ", (", vs, "), dp;\noption(redSB);\nideal I =\n"));
  for(i = 1, #E, write1(fname, Str(E[i], if(i < #E, ",\n", ";\n"))));
  write1(fname, "ideal J = std(I); \"dim:\"; dim(J); \"size:\"; size(J); if (size(J) < 10) { J; }\nquit;\n");
};
mkbf(7, 2, 1/3, -1/5, "pc7_mono_BF2.sing", 32003);
mkbf(7, 3, 1/3, -1/5, "pc7_mono_BF3.sing", 32003);
