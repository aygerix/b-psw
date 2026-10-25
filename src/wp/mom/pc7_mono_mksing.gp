\\ write the order-M (2,2)-pencil system (k = 7) for Singular: unknowns p1..p9 (w-moments, free), P1..PM (c-moments, free),
\\ lam, mu, r1, r2; relations lam Dl = 2 e1, mu Dl = 2 (1 - e1^2 + e2); Rabinowitsch y Dl (r1 - r2) = 1.
read("pc7_mono_hier.gp");
mks(k, M, fname, ch) = {
  my(Pv = vector(M, m, eval(Str("P", m))), E = heqs(k, r1, r2, lam, mu, M, Pv), e1 = r1 + r2, e2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2));
  my(vs = ""); forstep(a = 2*k - 1, 1, -1, vs = Str(vs, "p", a, ",")); forstep(m = M, 1, -1, vs = Str(vs, "P", m, ","));
  write1(fname, Str("ring R = ", ch, ", (y,", vs, "lam,mu,r1,r2), dp;\noption(redSB); option(prot);\nideal I =\n"));
  for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
  write1(fname, Str(lam*Dl - 2*e1, ",\n", mu*Dl - 2*(1 - e1^2 + e2), ",\n", y*Dl*(r1 - r2) - 1, ";\n"));
  write1(fname, "int t0 = timer; ideal J = std(I); \"dim:\"; dim(J); \"size:\"; size(J); \"seconds:\"; timer - t0;\n");
  write1(fname, "if (size(J) < 40) { J; } else { ideal E2 = eliminate(J, y*p1*p2*p3*p4*p5*p6*p7*p8*p9*lam*mu"); for(m = 1, M, write1(fname, Str("*P", m))); write1(fname, "); \"elim to r1,r2:\"; E2; }\nquit;\n");
  my(L = vector(#E, i, #Str(E[i]))); print(fname, ": ", #E, " equations, string lengths ", L);
};
mks(7, 2, "pc7_mono_M2_p.sing", 32003);
mks(7, 3, "pc7_mono_M3_p.sing", 32003);
mks(7, 3, "pc7_mono_M3_0.sing", 0);
mks(7, 4, "pc7_mono_M4_p.sing", 32003);
