\\ driver: level-2 Krawczyk certificate, k = 11, 11 cubics (numerical start from m9_lev2.py 11 0 2 16 67, trial 1)
\p 160
read("m9_kraw2.gp");
{
  my(m = 11, lin = 0, U0 = [2.0298944246460873, 0.5480388113614812, 0.6231902878442425, 0.5887296959880828, 0.927924080059095, -1.0242027689426119, 0.4665236748732248, 0.3582115404525946, 0.3961937258451305, 0.0207343600735967, 0.5011414135503955],
     V0 = [0.6888134530440482, 0.0523635305344567, 0.0644535936710554, 0.0559001369441409, 0.1556751115219284, 0.6881720530293298, 0.0535931296099456, 0.0676680940799653, 0.0539578840893951, 0.1474180567735788, 1.6515759852140512]);
  my(E = lev2eqs(m, lin), ALL = concat(vector(m, i, eval(Str("u", i))), vector(m, i, eval(Str("v", i)))), X0 = concat(U0, V0), r, ok2);
  print("==== LEVEL-2 Krawczyk certificate: k = 11, 11 cubics");
  print("  residuals of the 6 equations at the numerical start: ", precision(vector(6, j, substvec(E[j], ALL, X0)), 5));
  r = krawczykN(E, ALL, X0, 10^8, 6);
  ok2 = domcheck(m, r[2], r[6], U0, V0, 10^8);
  print("  CERTIFICATE ", if(r[1] && ok2, "VALID: a real admissible point satisfying orders 0, 1, 2 at t = -1 for BOTH halves EXISTS", "NOT VALID"));
}
