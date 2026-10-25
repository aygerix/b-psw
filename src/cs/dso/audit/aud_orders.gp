\\ Independent recomputation of (ord_p 2, z_D(p), (D/p)) for every census row, with a different z algorithm:
\\ z_D(p) = multiplicative order of gamma = alpha/beta (alpha, beta the roots of x^2 - x + Q), computed with
\\ fforder in F_p (split) or F_{p^2} (inert); ramified p: z = p (checked by U_p == 0 mod p).
\\ Also writes the rows as "p o z s" for aud_search.c.
\\ usage: echo 'IN="dso_X10000.txt"; DD=5; X=10000; OUTP="rows.txt"; read("aud_orders.gp")' | gp -q
QQ = (1 - DD)/4;
R = read(IN);
U(k, m) = lift((Mod([1, -QQ; 1, 0], m)^k)[2,1]);
zg(p) = {
  my(kr = kronecker(DD, p));
  if(kr == 0, if(U(p, p) != 0, error("ramified check ", p)); return(p));
  if(kr == 1,
    my(r = polrootsmod(x^2 - x + QQ, p)); return(znorder(Mod(r[1], p) / Mod(r[2], p))),
    my(a = ffgen(Mod(1, p)*(x^2 - x + QQ), 'a)); return(fforder(a / (1 - a))));
}
bad = 0; ni = 0; ns = 0; nr = 0; maxp = 0; fo = 0;
{
for(i = 1, #R,
  my(p = R[i][1], o = R[i][2], z = R[i][3], k = R[i][4]);
  if(!isprime(p), bad++; print("not prime ", p));
  my(o2 = if(p == 2, 0, znorder(Mod(2, p))), z2 = zg(p), k2 = kronecker(DD, p));
  if(o2 != o || z2 != z || k2 != k, bad++; print("MISMATCH ", R[i], " recomputed ", [o2, z2, k2]));
  if(o2 > X || z2 > X, bad++; print("ORDER > X ", R[i]));
  \\ definition check of z: U_z == 0, U_{z/r} != 0
  if(U(z2, p) != 0, bad++; print("U_z != 0 at ", p));
  my(fz = factor(z2)[,1]); for(j = 1, #fz, if(U(z2/fz[j], p) == 0, bad++; print("z not minimal at ", p)));
  if(k2 == -1, ni++, k2 == 1, ns++, nr++);
  maxp = max(maxp, p);
  write(OUTP, p, " ", o2, " ", z2, " ", k2));
}
print("rows ", #R, " inert ", ni, " split ", ns, " ramified ", nr, " max p ", maxp, " mismatches/bad ", bad);
