\\ AUDIT: rational M = p/q (|p|, q <= 300, gcd 1) where G5(y;M) has a rational root / splits completely over Q.
G(a,M) = (7*M^3+57*M^2+57*M+7)*a^3 - (98*M^3+308*M^2+98*M)*a^2 + (392*M^3+392*M^2)*a - 448*M^3;
{
  my(n1 = 0, n3 = 0, L3 = List(), L1 = List());
  for(q = 1, 300, for(p = -300, 300, if(p == 0 || gcd(p,q) != 1, next);
    my(M = p/q, P = G('y, M));
    if (poldegree(P) < 3, next);
    my(r = nfroots(, P));
    if (#r >= 1, n1++; listput(L1, [M, r]));
    if (#r == 3, n3++; listput(L3, [M, r]))));
  print("M with >= 1 rational root: ", n1);
  print(L1);
  print("M with full splitting: ", n3, "  ", L3);
}
