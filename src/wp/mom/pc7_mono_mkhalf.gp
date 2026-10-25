\\ HALF systems of the exact bilinear form: V-side only (F == 1 mod every f_j - 1) or W-side only.
default(parisizemax, 2000000000);
read("pc7_mono_exact.gp");
hh(k, ch, side, extra) = {my(nm = if(side == 1, "v", "w"));
  Str("ring R = ", ch, ", (y,z,", strjoin(vector(2*k - 2, i, Str(nm, 2*k - 2 - i)), ","), ",", strjoin(vector(k, m, Str("E", k + 1 - m)), ","), extra, "), dp;\noption(redSB);\nideal I =\n")};
tl = "int t0 = timer; ideal J = std(I); \"dim:\"; dim(J); \"size:\"; size(J); \"seconds:\"; timer - t0; if (size(J) < 10) { J; }\nquit;\n";
mkhf(k, side, r1, r2, fname, ch) = {
  my(E1 = r1 + r2, E2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2), E = exsys(k, E1, E2, 2*E1/Dl, 2*(1 - E1^2 + E2)/Dl));
  E = if(side == 1, E[1..4*k], E[4*k+1..8*k]); E = select(x -> x != 0, E);
  write1(fname, hh(k, ch, side, "")); for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
  write1(fname, Str("z*E", k, " - 1, y - 1;\n")); write1(fname, tl);
};
mkho(k, side, fname, ch) = {
  my(Dl = 1 - e1^2 + 2*e2 + e2^2, E = exsys(k, e1, e2, lam, mu));
  E = if(side == 1, E[1..4*k], E[4*k+1..8*k]); E = select(x -> x != 0, E);
  write1(fname, hh(k, ch, side, ",lam,mu,e1,e2")); for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
  write1(fname, Str(lam*Dl - 2*e1, ",\n", mu*Dl - 2*(1 - e1^2 + e2), ",\n", y*Dl*(e1^2 - 4*e2) - 1, ",\nz*E", k, " - 1;\n")); write1(fname, tl);
  print(fname, ": ", #E, " eqs, length ", sum(i = 1, #E, #Str(E[i])));
};
\\ mkhf(7, 1, 1/3, -1/5, "pc7_mono_HVF.sing", 32003);
\\ mkhf(7, 2, 1/3, -1/5, "pc7_mono_HWF.sing", 32003);
\\ mkhf(7, 1, -1/2, 1/2, "pc7_mono_HVFs.sing", 32003);
\\ mkhf(7, 2, -1/2, 1/2, "pc7_mono_HWFs.sing", 32003);
\\ mkho(7, 1, "pc7_mono_HVO_p.sing", 32003);
\\ mkho(7, 2, "pc7_mono_HWO_p.sing", 32003);
\\ mkho(7, 1, "pc7_mono_HVO_0.sing", 0);
\\ mkho(7, 2, "pc7_mono_HWO_0.sing", 0);
\\ slice-excluded versions: y e1 Dl (e1^2 - 4 e2) = 1
mkhoe(k, side, fname, ch) = {
  my(Dl = 1 - e1^2 + 2*e2 + e2^2, E = exsys(k, e1, e2, lam, mu));
  E = if(side == 1, E[1..4*k], E[4*k+1..8*k]); E = select(x -> x != 0, E);
  write1(fname, hh(k, ch, side, ",lam,mu,e1,e2")); for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
  write1(fname, Str(lam*Dl - 2*e1, ",\n", mu*Dl - 2*(1 - e1^2 + e2), ",\n", y*e1*Dl*(e1^2 - 4*e2) - 1, ",\nz*E", k, " - 1;\n")); write1(fname, tl);
  print(fname, ": ", #E, " eqs, length ", sum(i = 1, #E, #Str(E[i])));
};
mkhoe(7, 1, "pc7_mono_HVE_p.sing", 32003);
mkhoe(7, 2, "pc7_mono_HWE_p.sing", 32003);
mkhoe(7, 1, "pc7_mono_HVE_0.sing", 0);
mkhoe(7, 2, "pc7_mono_HWE_0.sing", 0);
mkhoe(5, 1, "pc7_mono_HVE5_0.sing", 0);
