\\ pc7_mono: decompose the k = 7 boundary expressions B1..B8 (q22_moments.gp) into monomial symmetric functions m_lambda(w),
\\ lambda with <= 3 parts. p_a -> x^a + y^a + z^a; coefficient of x^a y^b z^c (a >= b >= c) = coefficient of m_(a,b,c).
\\ Output: pc7_mono_tab.txt, one line per (i, a, b, c): coefficient as a polynomial in e1, e2, lam, mu.
read("q22_moments.gp");
{
  my(k = 7, B = bnd(k, e1, e2, lam, mu), out = "pc7_mono_tab.txt");
  my(P = vector(2*k, a, x^a + y^a + z^a));
  write1(out, "");
  for(i = 1, 8,
    my(Q = B[i]);
    for(a = 1, 2*k, Q = subst(Q, eval(Str("p", a)), P[a]));
    for(a = 0, 9, for(b = 0, a, for(c = 0, b,
      my(C = polcoef(polcoef(polcoef(Q, a, x), b, y), c, z));
      if(C != 0, write(out, i, " ", a, " ", b, " ", c, " ", C))))));
  print("done");
}
