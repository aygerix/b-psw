read("m9_lib.gp");
{
  my(r = mform(9, 0), Z = r[1], Zr = [reflect(Z[1]), reflect(Z[2])], Z0p = (Z[1] + Zr[1])/2, Z1p = (Z[2] + Zr[2])/2, Lc = 70*Z0p + 9*Z1p);
  print("poldegree(Lc, M1_4) = ", poldegree(Lc, M1_4), ";  Lc - Lc|_{M1_4=0} = ", Lc - subst(Lc, M1_4, 0));
  print("poldegree(Z0+, M1_4) = ", poldegree(Z0p, M1_4), "; Z0+ - Z0+|_{M1_4=0} = ", Z0p - subst(Z0p, M1_4, 0));
  \\ symmetric slice: odd symbols -> 0
  my(od = select(z -> eval(strsplit(Str(z), "_")[2]) % 2 == 1, syms(Lc)), Ls = substvec(Lc, od, vector(#od)), Z0s = substvec(Z0p, od, vector(#od)));
  print("70Z0+ + 9Z1+ on the symmetric slice: ", Ls);
  print("Z0+ on the symmetric slice: ", Z0s);
}
quit;
