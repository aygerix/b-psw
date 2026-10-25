\\ order-M system with e1, e2 both symbolic, slice e1 = 0 EXCLUDED: Rabinowitsch y e1 Dl (e1^2 - 4 e2) = 1.
\\ Unit ideal => no (2,2) pencil (k = 7, any nonzero c_j) off the symmetric slice, for Dl != 0, r1 != r2.
read("pc7_mono_hs.gp");
mke(k, M, fname, ch) = {
  my(Pv = vector(M, m, eval(Str("P", m))), E = hsym(k, e1, e2, lam, mu, M, Pv), Dl = 1 - e1^2 + 2*e2 + e2^2);
  my(vs = strjoin(concat(vector(2*k - 1, a, Str("p", 2*k - a)), vector(M, m, Str("P", M + 1 - m))), ","));
  write1(fname, Str("ring R = ", ch, ", (y,", vs, ",lam,mu,e1,e2), dp;\noption(redSB);\nideal I =\n"));
  for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
  write1(fname, Str(lam*Dl - 2*e1, ",\n", mu*Dl - 2*(1 - e1^2 + e2), ",\n", y*e1*Dl*(e1^2 - 4*e2) - 1, ";\n"));
  write1(fname, "int t0 = timer; ideal J = std(I); \"dim:\"; dim(J); \"size:\"; size(J); \"seconds:\"; timer - t0; if (size(J) < 10) { J; }\nquit;\n");
};
mke(7, 3, "pc7_mono_E3_p.sing", 32003);
mke(7, 3, "pc7_mono_E3_0.sing", 0);
mke(7, 4, "pc7_mono_E4_p.sing", 32003);
