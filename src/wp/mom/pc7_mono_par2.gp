\\ PARAMETRIC order-M systems (k = 7): one shape coordinate in the coefficient field.
\\  (a) field Q(e1) (or F_p(e1)): vars y, p9..p1, PM..P1, lam, mu, e2; relations lam Dl = 2 e1, mu Dl = 2(1 - e1^2 + e2); y Dl (e1^2 - 4 e2) = 1.
\\  (b) field Q(e2): vars y, p's, P's, lam, mu, e1; Rabinowitsch y e1 Dl (e1^2 - 4 e2) = 1 (slice excluded).
\\ Only p1..p9 occur (V, W have degree 12).
read("pc7_mono_hs.gp");
mkpar(k, M, par, fname, ch) = {
  my(Pv = vector(M, m, eval(Str("P", m))), E = hsym(k, e1, e2, lam, mu, M, Pv), Dl = 1 - e1^2 + 2*e2 + e2^2, np = 2*k - 5);
  my(vs = strjoin(concat(vector(np, a, Str("p", np + 1 - a)), vector(M, m, Str("P", M + 1 - m))), ","));
  my(other = if(par == "e1", "e2", "e1"), rab = if(par == "e1", y*Dl*(e1^2 - 4*e2) - 1, y*e1*Dl*(e1^2 - 4*e2) - 1));
  write1(fname, Str("ring R = (", ch, ",", par, "), (y,", vs, ",lam,mu,", other, "), dp;\noption(redSB);\nideal I =\n"));
  for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
  write1(fname, Str(lam*Dl - 2*e1, ",\n", mu*Dl - 2*(1 - e1^2 + e2), ",\n", rab, ";\n"));
  write1(fname, "int t0 = timer; ideal J = std(I); \"dim:\"; dim(J); \"size:\"; size(J); \"seconds:\"; timer - t0; if (size(J) < 10) { J; }\nquit;\n");
};
mkpar(7, 3, "e1", "pc7_mono_Pe1_p.sing", 32003);
mkpar(7, 3, "e2", "pc7_mono_Pe2_p.sing", 32003);
mkpar(7, 3, "e1", "pc7_mono_Pe1_0.sing", 0);
mkpar(7, 3, "e2", "pc7_mono_Pe2_0.sing", 0);
