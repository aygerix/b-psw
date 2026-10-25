\\ Line OC, family 1: Mersenne primes M_q. Fermat side forced: ord = q, j1 = 0.
\\ For each Selfridge D: inert/split; if inert, p + 1 = 2^q so z = 2^j2 (odd part 1: free Lucas side);
\\ j2 found by the V-doubling chain V_{2^r} (first r with V_{2^r} == 0 gives z = 2^(r+1)).
\\ Prediction (Lemma 18): (Q/p) = -1 => j2 = q exactly; (Q/p) = +1 => j2 <= q - 1 (Kummer depth t = q-1-j2).
read("oc_lib.gp");
mexp = [3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701];
if(type(QMAX) != "t_INT", QMAX = 21701);
j2chain(p, D, q) = { my(Q = Mod(QofD(D), p), V = Mod(1, p), Qk = Q);
  for(r = 0, q - 1, if(V == 0, return(r + 1)); V = V^2 - 2*Qk; Qk = Qk^2);
  if(V == 0, return(q + 1)); -1; }   \\ -1: no zero found (would mean z odd => z = 1, impossible)
{
for(i = 1, #mexp, q = mexp[i]; if(q > QMAX, break); p = 2^q - 1; t0 = getabstime();
  s = Str("q=", q, " q%4=", q%4, " q%3=", q%3, " |");
  for(k = 1, #SD, D = SD[k]; e = kronecker(D, p);
    if(e == 0, s = Str(s, " D", D, ":ram"); next);
    if(e == 1, s = Str(s, " D", D, ":split"); next);
    j = j2chain(p, D, q); qr = kronecker(QofD(D), p);
    s = Str(s, " D", D, ":j2=", j, "(Q/p=", qr, ")"));
  print(s, " | ms=", getabstime() - t0));
}
