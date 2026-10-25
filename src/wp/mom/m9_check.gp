read("m9_moments.gp");
{
  foreach([[5,0],[4,1],[7,0],[6,1],[9,0],[8,1],[11,0],[10,1]], c,
    my(t0 = getabstime(), r = checkid(c[1], c[2]));
    print("m=", c[1], " lin=", c[2], " identity residuals: ", r, "   (", getabstime() - t0, " ms)"));
  for(k = 9, 11, if(k % 2 == 0, next); my(Ls = Lsym(k), TF = Tform(k, Ls));
    print("k=", k, " T(-1) = ", TF[1]); print("k=", k, " T'(-1) = ", TF[2]));
  my(k = 11, uv = ellv(u, v, 11)); for(n = 3, 10, print("ell_", n, " = ", factor(uv[n])));
}
