\\ k = 7 (2,2) pencil: boundary system with lam, mu substituted (denominator Dl cleared), p's free; Singular elimination to (e1, e2)
read("q22_moments.gp");
mkel(k, fname, ch) = {
  my(B = bnd(k, e1, e2, lam, mu), np = 2*k - 1, Dl = 1 - e1^2 + 2*e2 + e2^2);
  my(vs = ""); forstep(a = np, 1, -1, vs = Str(vs, "p", a, ","));
  write1(fname, Str("ring R = ", ch, ", (", vs, "lam,mu,e1,e2), dp;\noption(redSB);\nideal I =\n"));
  for(i = 1, #B, write1(fname, Str(B[i], ",\n")));
  write1(fname, "lam*(1 - e1^2 + 2*e2 + e2^2) - 2*e1,\nmu*(1 - e1^2 + 2*e2 + e2^2) - 2*(1 - e1^2 + e2);\n");
  write1(fname, "int t0 = timer; ideal E = eliminate(I, "); forstep(a = np, 1, -1, write1(fname, Str("p", a, "*"))); write1(fname, "lam*mu);\n");
  write1(fname, "\"elimination ideal in (e1,e2):\"; E; \"seconds:\"; timer - t0;\n");
  \\ saturate by Dl (degenerate locus H = G excluded): Dl invertible
  write1(fname, "ring S = (", ch, "), (e1,e2,y), dp; ideal E2 = imap(R, E); E2 = E2, y*(1 - e1^2 + 2*e2 + e2^2) - 1; ideal F = eliminate(std(E2), y); \"after removing Dl = 0:\"; F; quit;\n");
};
mkel(7, "q22e_k7_p.sing", 32003);
mkel(7, "q22e_k7_0.sing", 0);
mkel(5, "q22e_k5_0.sing", 0);
