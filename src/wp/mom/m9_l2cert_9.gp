\\ driver: LEVEL-2 Krawczyk certificates for k = 9 (both cases), large-S numerical starts (m9_lev2.py / m9_climb.py)
\p 160
read("m9_kraw2.gp");
{
  my(cases = [[9, 0, [-1.5022143532841334, -4.123908001403858, 0.5139699491914088, -6.9488109387036845, 0.4557359591171056, 1.9898717712453684, 8.384115532862122, 0.0300074710797625, 0.2116870195358351],
                      [2.1953624588934777, 4.160359143007728, 0.1638235274005706, 12.93704174320486, 0.1504054295230128, 11.332431645079225, 13.503906596594524, 0.8344275287407846, 16.583303646561404]],
               [8, 1, [7.383602524018461, -0.281220887363323, 3.7362686543302996, -6.113476957120817, -1.9894291170655451, 0.9158133556559549, -1.4477281105089814, -2.7740334229151604],
                      [10.187126685349437, 7.09708766558695, 3.420275946834393, 9.28262974093801, 1.3912245328824524, 13.036011292752544, 1.249721086956731, 4.659440794003997]]]);
  foreach(cases, c, my(m = c[1], lin = c[2], U0 = c[3], V0 = c[4], E = lev2eqs(m, lin),
      ALL = concat(vector(m, i, eval(Str("u", i))), vector(m, i, eval(Str("v", i)))), X0 = concat(U0, V0), r, ok2);
    print("==== LEVEL-2 Krawczyk certificate: k = 9, cubic members ", m, ", linear ", lin);
    print("  residuals of the 6 equations at the numerical start: ", precision(vector(6, j, substvec(E[j], ALL, X0)), 5));
    r = krawczykN(E, ALL, X0, 10^8, 6);
    ok2 = domcheck(m, r[2], r[6], U0, V0, 10^8);
    print("  CERTIFICATE ", if(r[1] && ok2, "VALID: a real admissible point satisfying orders 0, 1, 2 at t = -1 for BOTH halves EXISTS", "NOT VALID")));
}
