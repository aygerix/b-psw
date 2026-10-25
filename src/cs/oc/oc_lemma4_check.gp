\\ Lemma OC4 check: inert p == 3 mod 4, (Q/p) = 1, s^2 = Q mod p.
\\ Claim: v2(z(p)) = v2(p+1) - 1  iff  ((2Q + s)/p) = -1 ; and the symbol does not depend on the sign of s.
\\ Also Lemma OC3 sanity: for inert p, min(v2(ord_p 2), v2(z(p))) <= 1.
read("oc_fam.gp");
Dl = [5,-7,-11,13,-15,17,-19,-23,29,-31,37,41,-43,-47];
{
for(k = 1, #Dl, D = Dl[k]; Q = QofD(D); ntest = 0; nbad = 0; nsign = 0; n3 = 0; n3bad = 0;
  forprime(p = 3, 3*10^6, if((2*D*Q) % p == 0, next); if(kronecker(D, p) != -1, next);
    j1 = j1of(p); j2 = j2of(p, D, p + 1); n3++; if(min(j1, j2) > 1, n3bad++);
    if(p % 4 != 3 || kronecker(Q, p) != 1, next);
    s = sqrt(Mod(Q, p)); a = kronecker(lift(2*Q + s), p); b = kronecker(lift(2*Q - s), p);
    if(a != b, nsign++);
    full1 = (j2 == valuation(p + 1, 2) - 1);
    ntest++; if(full1 != (a == -1), nbad++));
  printf("D=%4d Q=%3d  OC4 tested %6d  mismatches %d  sign-dependence %d | OC3 inert %6d violations %d\n", D, Q, ntest, nbad, nsign, n3, n3bad));
}
