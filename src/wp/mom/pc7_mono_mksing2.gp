\\ Singular scripts for the symmetric order-M system (k = 7): unknowns p1..p9, P1..PM free; lam, mu, e1, e2 with
\\ lam Dl = 2 e1, mu Dl = 2(1 - e1^2 + e2), Dl = 1 - e1^2 + 2 e2 + e2^2; Rabinowitsch y Dl (e1^2 - 4 e2) = 1 (r1 != r2, H != G).
read("pc7_mono_hs.gp");
mks(k, M, fname, ch, stdonly = 0) = {
  my(Pv = vector(M, m, eval(Str("P", m))), E = hsym(k, e1, e2, lam, mu, M, Pv), Dl = 1 - e1^2 + 2*e2 + e2^2);
  my(vs = ""); forstep(a = 2*k - 1, 1, -1, vs = Str(vs, "p", a, ",")); forstep(m = M, 1, -1, vs = Str(vs, "P", m, ","));
  write1(fname, Str("ring R = ", ch, ", (y,", vs, "lam,mu,e1,e2), dp;\noption(redSB);\nideal I =\n"));
  for(i = 1, #E, write1(fname, Str(E[i], ",\n")));
  write1(fname, Str(lam*Dl - 2*e1, ",\n", mu*Dl - 2*(1 - e1^2 + e2), ",\n", y*Dl*(e1^2 - 4*e2) - 1, ";\n"));
  write1(fname, "int t0 = timer; ideal J = std(I); \"dim:\"; dim(J); \"size:\"; size(J); \"seconds:\"; timer - t0;\n");
  write1(fname, "if (size(J) < 60) { J; } else { ideal E2 = eliminate(J, y*p1*p2*p3*p4*p5*p6*p7*p8*p9*lam*mu"); for(m = 1, M, write1(fname, Str("*P", m))); write1(fname, "); \"elim to e1,e2:\"; E2; }\nquit;\n");
  print(fname, ": ", #E, " equations, total length ", sum(i = 1, #E, #Str(E[i])));
};
my(t0 = getabstime());
mks(7, 2, "pc7_mono_S2_p.sing", 32003); print("t ", getabstime() - t0);
mks(7, 3, "pc7_mono_S3_p.sing", 32003); print("t ", getabstime() - t0);
mks(7, 3, "pc7_mono_S3_0.sing", 0); print("t ", getabstime() - t0);
mks(7, 4, "pc7_mono_S4_p.sing", 32003); print("t ", getabstime() - t0);
