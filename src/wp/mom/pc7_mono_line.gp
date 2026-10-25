\\ order-M systems on a LINE of shapes: e1 fixed rational, e2 free; lam, mu variables with lam Dl = 2 e1, mu Dl = 2(1 - e1^2 + e2);
\\ Rabinowitsch y Dl (e1^2 - 4 e2) = 1.  Elimination to e2.
read("pc7_mono_hs.gp");
mkl(k, M, E1, fname, ch) = {
  my(Pv = vector(M, m, eval(Str("P", m))), E = hsym(k, E1, e2, lam, mu, M, Pv), Dl = 1 - E1^2 + 2*e2 + e2^2);
  my(vs = strjoin(concat(vector(2*k - 1, a, Str("p", 2*k - a)), vector(M, m, Str("P", M + 1 - m))), ","));
  write1(fname, Str("ring R = ", ch, ", (y,", vs, ",lam,mu,e2), dp;\noption(redSB);\nideal I =\n"));
  for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
  write1(fname, Str(lam*Dl - 2*E1, ",\n", mu*Dl - 2*(1 - E1^2 + e2), ",\n", y*Dl*(E1^2 - 4*e2) - 1, ";\n"));
  write1(fname, "int t0 = timer; ideal J = std(I); \"dim:\"; dim(J); \"vdim:\"; vdim(J); \"size:\"; size(J); \"seconds:\"; timer - t0;\n");
  write1(fname, Str("ideal E = eliminate(J, y*", strjoin(concat(vector(2*k - 1, a, Str("p", a)), vector(M, m, Str("P", m))), "*"), "*lam*mu); \"elim to e2:\"; E; \"seconds:\"; timer - t0; if (size(E) > 0) { factorize(E[1]); }\nquit;\n"));
};
mkl(7, 2, 2/15, "pc7_mono_L2.sing", 32003);
mkl(7, 3, 2/15, "pc7_mono_L3.sing", 32003);
mkl(7, 2, 2/15, "pc7_mono_L2_0.sing", 0);
