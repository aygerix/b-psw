\\ vrf line: symbolic check of the mechanism statement (uses the claimant's m9_lib.gp moment form)
read("m9_lib.gp");
{
  my(r = mform(9, 0), Z = r[1], Zr = [reflect(Z[1]), reflect(Z[2])], Z0p = (Z[1] + Zr[1])/2, Z1p = (Z[2] + Zr[2])/2, Lc = 70*Z0p + 9*Z1p);
  print("symbols in Z0+: ", syms(Z0p));
  print("symbols in 70Z0+ + 9Z1+: ", syms(Lc));
  print("coefficient of M1_4 in 70Z0+ + 9Z1+: ", polcoef(Lc, 1, M1_4), "   in Z0+: ", polcoef(Z0p, 1, M1_4));
  \\ even (b even) symbols are linear in y except b = 4 (quadratic), b = 6 (cubic), ...
  my(ev = select(z -> eval(strsplit(Str(z), "_")[2]) % 2 == 0 && eval(strsplit(Str(z), "_")[2]) >= 4, syms(Lc)));
  print("even symbols with b >= 4 in 70Z0+ + 9Z1+: ", ev);
}
quit;
