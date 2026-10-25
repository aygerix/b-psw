\\ j2 = v2(z_D(M_q)) (inert, via V-doubling chain) or split data, for the Selfridge sequence D (non-square) |D| <= 103.
read("oc_lib.gp");
mexp = [3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423];
j2chain(p, D, q) = { my(Q = Mod(QofD(D), p), V = Mod(1, p), Qk = Q);
  for(r = 0, q - 1, if(V == 0, return(r + 1)); V = V^2 - 2*Qk; Qk = Qk^2); -1; }
\\ split: z | p - 1 = 2 m (m odd): j2 = 0 if U_m == 0 else 1
j2split(p, D, q) = { my(m = (p - 1)/2, a = Mod(Mod(1,p)*x, x^2 - x + QofD(D))^m); if(polcoeff(lift(lift(a)),1) == 0, 0, 1); }
{
D = 5; while(abs(D) <= 103,
  if(!issquare(D), for(i = 1, #mexp, q = mexp[i]; p = 2^q - 1; e = kronecker(D, p);
      j = if(e == 0, -9, if(e == -1, j2chain(p, D, q), j2split(p, D, q)));
      print(D, " ", q, " ", e, " ", j)));
  D = if(D > 0, -D - 2, -D + 2));
}
