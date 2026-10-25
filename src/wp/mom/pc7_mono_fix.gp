\\ fixed-shape order-M systems (exact rational shape) for Singular: unknowns p1..p9, P1..PM only.
read("pc7_mono_hs.gp");
mkf(k, M, r1, r2, fname, ch) = {
  my(E1 = r1 + r2, E2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2), lam0 = 2*E1/Dl, mu0 = 2*(1 - E1^2 + E2)/Dl);
  my(Pv = vector(M, m, eval(Str("P", m))), E = hsym(k, E1, E2, lam0, mu0, M, Pv));
  my(vs = ""); forstep(a = 2*k - 1, 1, -1, vs = Str(vs, "p", a, ",")); forstep(m = M, 1, -1, vs = Str(vs, "P", m, ","));
  vs = strjoin(concat(vector(2*k - 1, a, Str("p", 2*k - a)), vector(M, m, Str("P", M + 1 - m))), ",");
  write1(fname, Str("ring R = ", ch, ", (", vs, "), dp;\noption(redSB);\nideal I =\n"));
  for(i = 1, #E, write1(fname, Str(E[i], if(i < #E, ",\n", ";\n"))));
  write1(fname, "int t0 = timer; ideal J = std(I); \"dim:\"; dim(J); \"vdim:\"; vdim(J); \"size:\"; size(J); \"seconds:\"; timer - t0; if (size(J) < 15) { J; }\nquit;\n");
};
mkf(7, 1, 1/3, -1/5, "pc7_mono_F1.sing", 32003);
mkf(7, 2, 1/3, -1/5, "pc7_mono_F2.sing", 32003);
mkf(7, 3, 1/3, -1/5, "pc7_mono_F3.sing", 32003);
mkf(7, 2, 2/7, 3/5, "pc7_mono_F2b.sing", 32003);
mkf(7, 2, 1/3, -1/5, "pc7_mono_F2_0.sing", 0);
