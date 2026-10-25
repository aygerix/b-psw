\\ Independent Lucas-Wieferich check: for p in {1093, 3511} and every Selfridge D (5, -7, 9, -11, ...) with |D| <= DMAX,
\\ p not dividing D Q: z(p) = least divisor m of p - (D/p) with U_m(1,Q) == 0 mod p; hit iff U_z == 0 mod p^2.
\\ usage: echo 'DMAX=1000; read("aud_lw.gp")' | gp -q
U(k, Q, m) = lift((Mod([1, -Q; 1, 0], m)^k)[2,1]);
{
foreach([1093, 3511], p,
  my(hits = List(), nt = 0, sk = List(), D = 5);
  while(abs(D) <= DMAX,
    my(Q = (1 - D)/4);
    if(D % p == 0 || Q % p == 0, listput(sk, D),
      nt++;
      my(kr = kronecker(D, p), z = 0);
      fordiv(p - kr, m, if(U(m, Q, p) == 0, z = m; break));
      if(z == 0, error("no z ", D));
      if(U(z, Q, p^2) == 0, listput(hits, [D, z, kr, issquare(D)])));
    D = if(D > 0, -D - 2, -D + 2));
  print("p = ", p, ": D tested ", nt, ", skipped ", Vec(sk), ", hits [D, z, (D/p), issquare] ", Vec(hits)));
}
