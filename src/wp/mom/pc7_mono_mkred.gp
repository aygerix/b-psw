\\ Singular script for the REDUCED order-M system (pc7_mono_red.gp), shape symbolic, slice excluded (y lam Dl (e1^2 - 4 e2) = 1).
read("pc7_mono_red.gp");
mkr(k, M, fname, ch) = {
  my(Pv = vector(M, m, eval(Str("P", m))), R = redsys(k, M, e1, e2, lam, mu, Pv), E = R[1], av = R[3], Dl = 1 - e1^2 + 2*e2 + e2^2);
  E = vector(#E, i, numerator(E[i]));
  my(vs = strjoin(concat(vector(#av, i, Str(av[i])), vector(M, m, Str("P", M + 1 - m))), ","));
  write1(fname, Str("ring R = ", ch, ", (y,", vs, ",lam,mu,e1,e2), dp;\noption(redSB);\nideal I =\n"));
  for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
  write1(fname, Str(lam*Dl - 2*e1, ",\n", mu*Dl - 2*(1 - e1^2 + e2), ",\n", y*lam*Dl*(e1^2 - 4*e2) - 1, ";\n"));
  write1(fname, "int t0 = timer; ideal J = std(I); \"dim:\"; dim(J); \"size:\"; size(J); \"seconds:\"; timer - t0; if (size(J) < 10) { J; }\nquit;\n");
  print(fname, ": ", #E, " eqs, lengths ", vector(#E, i, #Str(E[i])), " denominators ", vector(#R[1], i, factor(denominator(R[1][i]))[,1]~));
};
mkr(7, 3, "pc7_mono_R3_p.sing", 32003);
mkr(7, 3, "pc7_mono_R3_0.sing", 0);
