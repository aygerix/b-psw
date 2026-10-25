\\ Lemma OC4 on Mersenne primes: for every inert (D, M_q) with (Q/M_q) = 1 in the j2 table,
\\ t := q - 1 - j2 ; predicted t = 0 iff ((2Q + s)/M_q) = -1, s = sqrt(Q) mod M_q.
\\ Closed forms (Theorem M7 and its D = -31 analogue):
\\   D = -7  (Q = 2): s = 2^((q+1)/2), symbol = (2^m + 1 / 7) with m = (q-3)/2  -> always -1 for q == 5 mod 6 (q >= 7)
\\   D = -31 (Q = 8): s = 2^((q+3)/2), symbol = (2^m + 1 mod 31 / 31) with m = (q-5)/2 (q >= 9)
read("oc_lib.gp");
tab = readvec("out/j2table.gp.txt");
{
ntest = 0; nbad = 0; c7 = 0; c7bad = 0; c31 = 0; c31bad = 0;
for(i = 1, #tab, [D, q, e, j] = tab[i]; if(e != -1 || q < 5, next); Q = QofD(D); p = 2^q - 1;
  if(kronecker(Q, p) != 1, if(j != q, print("Lemma18 violation ", D, " ", q)); next);
  s = sqrt(Mod(Q, p)); sym = kronecker(lift(2*Q + s), p); t = q - 1 - j;
  ntest++; if((t == 0) != (sym == -1), nbad++; print("OC4 MISMATCH D=", D, " q=", q, " t=", t, " sym=", sym));
  if(D == -7, c7++; m = (q-3)/2; pred = if(q >= 7, kronecker(lift(Mod(2,7)^m + 1), 7), sym); if(pred != sym || t != 0, c7bad++));
  if(D == -31 && q >= 9, c31++; m = (q-5)/2; pred = kronecker(lift(Mod(2,31)^m + 1), 31); if(pred != sym, c31bad++)));
printf("OC4 on Mersenne table: %d inert (D,q) with (Q/p)=1 tested, %d mismatches\n", ntest, nbad);
printf("D=-7 closed form (t = 0 always): %d cases, %d failures;  D=-31 closed form: %d cases, %d failures\n", c7, c7bad, c31, c31bad);
}
