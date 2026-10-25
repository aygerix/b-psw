\\ m9vmf_neg.gp -- negative controls (sensitivity of m9vmf checks) and k = 5, 7 check of the author's mform (Lam route)
\\ against the independent direct-division E-side.
read("m9_lib.gp");
read("m9vmf_fun.gp");
{
  setrand(7);
  my(CS = parseout("m9_mform.out"));
  my(pert = [[4, M1_0^2/1000, "Z0+ += S^2/1000"], [6, M3_1/7, "Z0- += M3_1/7"], [5, 1/1024, "Z1+ += 1/1024"],
             [7, M1_1*M1_0/3, "Z1- += S*M11/3"], [4, M1_0^3/99, "Z0+ += S^3/99"]]);
  foreach([1, 3], ci, my(cs0 = CS[ci]);
    print("== negative controls, case ", cs0[1..3]);
    foreach(pert, p, my(cs = cs0); cs[p[1]] += p[2];
      print("  ", p[3], ":  (A) mismatching points ", partA(cs, 5) ,"/5;  (B, r=3) residual zero? ", partB(cs, 3) == [0,0,0,0])));
  \\ k = 5, 7 via author's mform (Lam route), compared with independent E-side
  foreach([[5,0],[4,1],[7,0],[6,1],[9,0],[8,1]], c, my(m = c[1], lin = c[2], k = m + lin, Z = mform(m, lin)[1]);
    my(cs = [k, m, lin, Z[1], Z[2], 0, 0], R = max((k-1)\3, vecmax([mdeg(Z[1]), mdeg(Z[2])])));
    print("== author mform(", m, ",", lin, "), k = ", k, ": R = ", R, ";  (A) mismatches ", partA(cs, 20), "/20;  (B) r = ", R, "..", R+1, ": ",
      vector(2, j, if(k - 1 - (R + j - 1) - lin >= 0, partB(cs, R + j - 1), "skip"))));
}
