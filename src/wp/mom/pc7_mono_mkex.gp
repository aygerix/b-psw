\\ Singular scripts for the exact bilinear system.
default(parisizemax, 2000000000);
read("pc7_mono_exact.gp");
hdr(k, ch, extra) = {Str("ring R = ", ch, ", (y,z,", strjoin(vector(2*k - 2, i, Str("v", 2*k - 2 - i)), ","), ",",
  strjoin(vector(2*k - 2, i, Str("w", 2*k - 2 - i)), ","), ",", strjoin(vector(k, m, Str("E", k + 1 - m)), ","), extra, "), dp;\noption(redSB);\nideal I =\n")};
tail0 = "int t0 = timer; ideal J = std(I); \"dim:\"; dim(J); \"size:\"; size(J); \"seconds:\"; timer - t0; if (size(J) < 10) { J; }\nquit;\n";
\\ fixed shape
mkexf(k, r1, r2, fname, ch) = {
  my(E1 = r1 + r2, E2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2), E = exsys(k, E1, E2, 2*E1/Dl, 2*(1 - E1^2 + E2)/Dl));
  E = select(x -> x != 0, E);
  write1(fname, hdr(k, ch, "")); for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
  write1(fname, Str("z*E", k, " - 1, y - 1;\n")); write1(fname, tail0);
  print(fname, ": ", #E, " nonzero equations");
};
\\ light-cone shape, full 2-parameter; slice NOT excluded; y (e1^2 - 4 e2) = 1, z E_k = 1, hp bp = hm bm = 1
mkexh(k, fname, ch) = {
  my(S1 = (hm - hp)/2, S2 = (hp + hm)/2 - 1, SL = bp - bm, SM = bp + bm - (bp - bm)*(hm - hp)/2);
  my(E = exsys(k, S1, S2, SL, SM)); E = select(x -> x != 0, E);
  write1(fname, hdr(k, ch, ",bp,bm,hp,hm")); for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
  write1(fname, Str(hp*bp - 1, ",\n", hm*bm - 1, ",\n", y*(S1^2 - 4*S2) - 1, ",\nz*E", k, " - 1;\n")); write1(fname, tail0);
  print(fname, ": ", #E, " nonzero equations, total length ", sum(i = 1, #E, #Str(E[i])));
};





\\ line e1 = c (hm = hp + 2c) in light-cone coordinates
mkexl(k, c, fname, ch) = {
  my(HM = hp + 2*c, S1 = c, S2 = (hp + HM)/2 - 1, SL = bp - bm, SM = bp + bm - (bp - bm)*(HM - hp)/2);
  my(E = exsys(k, S1, S2, SL, SM)); E = select(x -> x != 0, E);
  write1(fname, hdr(k, ch, ",bp,bm,hp")); for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
  write1(fname, Str(hp*bp - 1, ",\n", HM*bm - 1, ",\n", y*(S1^2 - 4*S2) - 1, ",\nz*E", k, " - 1;\n")); write1(fname, tail0);
  print(fname, ": ", #E, " nonzero equations, total length ", sum(i = 1, #E, #Str(E[i])));
};
\\ original symbols with relations (lam Dl = 2 e1, mu Dl = 2(1 - e1^2 + e2)), y Dl (e1^2 - 4 e2) = 1
mkexo(k, fname, ch) = {
  my(Dl = 1 - e1^2 + 2*e2 + e2^2, E = exsys(k, e1, e2, lam, mu)); E = select(x -> x != 0, E);
  write1(fname, hdr(k, ch, ",lam,mu,e1,e2")); for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
  write1(fname, Str(lam*Dl - 2*e1, ",\n", mu*Dl - 2*(1 - e1^2 + e2), ",\n", y*Dl*(e1^2 - 4*e2) - 1, ",\nz*E", k, " - 1;\n")); write1(fname, tail0);
  print(fname, ": ", #E, " nonzero equations, total length ", sum(i = 1, #E, #Str(E[i])));
};
mkexl(7, 1/15, "pc7_mono_XL_p.sing", 32003);
mkexl(7, 0, "pc7_mono_XLs_p.sing", 32003);
mkexo(7, "pc7_mono_XO_p.sing", 32003);
mkexo(7, "pc7_mono_XO_0.sing", 0);
mkexo(5, "pc7_mono_XO5_0.sing", 0);
