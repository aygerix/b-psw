\\ write the (2,2)-pencil boundary system for k members (power-sum form) as a Singular script
read("q22_moments.gp");
mksing(k, fname, quad = 0) = {
  my(B = if(quad, bnd(k, 0, e2, 0, mu, 1), bnd(k, e1, e2, lam, mu)), np = 2*k - 1);
  my(vs = if(quad, "aa,e1,e2,lam,mu", "e1,e2,lam,mu")); for(a = 1, np, vs = Str(vs, ",p", a));
  write1(fname, Str("ring R = 0, (", vs, "), dp;\nideal I =\n"));
  for(i = 1, #B, write1(fname, Str(B[i], ",\n")));
  if(quad, write1(fname, "e1, lam, aa*(1 + e2) - 2,\n"));
  write1(fname, "lam*(1 - e1^2 + 2*e2 + e2^2) - 2*e1,\nmu*(1 - e1^2 + 2*e2 + e2^2) - 2*(1 - e1^2 + e2);\n");
  write1(fname, "option(redSB); ideal J = std(I); \"dim:\"; dim(J); \"size:\"; size(J); if (size(J) < 30) { J; }\n");
  write1(fname, "\"elim to e1,e2:\"; ideal E = eliminate(J, lam*mu"); for(a = 1, np, write1(fname, Str("*p", a))); write1(fname, "); E;\nquit;\n");
};
mksing(5, "q22m_k5.sing"); mksing(5, "q22m_k5q.sing", 1); mksing(7, "q22m_k7.sing"); mksing(7, "q22m_k7q.sing", 1);
