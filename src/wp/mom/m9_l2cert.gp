\\ driver: level-2 Krawczyk certificate, k = 11, 10 cubics + 2t+1 (numerical start from m9_lev2.py, trial 0)
\p 160
read("m9_kraw2.gp");
{
  my(m = 10, lin = 1, U0 = [-1.1244994674169049, 0.5084007433164041, 0.5853893824692474, 0.4272437189095763, 0.4658918552026194, 2.038708594947349, 0.4585185161079019, 0.4019156901346113, 1.3113093911144942, -0.0813873425213251],
     V0 = [0.6072983191125251, 0.0333256491510538, 0.0433496341649076, 0.0433484130288553, 0.0394985389031685, 0.5215586183929436, 1.5089365600228837, 0.0467015076228214, 0.5069716337365633, 0.3615705499023993]);
  my(E = lev2eqs(m, lin), ALL = concat(vector(m, i, eval(Str("u", i))), vector(m, i, eval(Str("v", i)))), X0 = concat(U0, V0), r, ok2);
  print("==== LEVEL-2 Krawczyk certificate: k = 11, 10 cubics + 2t+1");
  print("  residuals of the 6 equations at the numerical start: ", precision(vector(6, j, substvec(E[j], ALL, X0)), 5));
  r = krawczykN(E, ALL, X0, 10^8, 6);
  ok2 = domcheck(m, r[2], r[6], U0, V0, 10^8);
  print("  CERTIFICATE ", if(r[1] && ok2, "VALID: a real admissible point satisfying orders 0, 1, 2 at t = -1 for BOTH halves EXISTS", "NOT VALID"));
}
