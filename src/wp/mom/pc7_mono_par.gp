\\ order-M system over the parameter field Q(e1) (resp. F_p(e1)): Singular ring (ch, e1), unknowns y, p9..p1, PM..P1, lam, mu, e2.
\\ 1 in the ideal over Q(e1)  <=>  some nonzero g(e1) lies in the ideal (no solutions for e1 outside the roots of g).
read("pc7_mono_hs.gp");
mkp(k, M, fname, ch) = {
  my(Pv = vector(M, m, eval(Str("P", m))), E = hsym(k, e1, e2, lam, mu, M, Pv), Dl = 1 - e1^2 + 2*e2 + e2^2);
  my(vs = strjoin(concat(vector(2*k - 1, a, Str("p", 2*k - a)), vector(M, m, Str("P", M + 1 - m))), ","));
  write1(fname, Str("ring R = (", ch, ",e1), (y,", vs, ",lam,mu,e2), dp;\noption(redSB);\nideal I =\n"));
  for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
  write1(fname, Str(lam*Dl - 2*e1, ",\n", mu*Dl - 2*(1 - e1^2 + e2), ",\n", y*Dl*(e1^2 - 4*e2) - 1, ";\n"));
  write1(fname, "int t0 = timer; ideal J = std(I); \"dim:\"; dim(J); \"size:\"; size(J); \"seconds:\"; timer - t0; if (size(J) < 10) { J; }\nquit;\n");
};
mkp(7, 3, "pc7_mono_P3_p.sing", 32003);
mkp(7, 3, "pc7_mono_P3_0.sing", 0);
