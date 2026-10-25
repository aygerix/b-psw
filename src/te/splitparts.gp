\\ Split parts S for Theorem E: every product S of admissible split primes q <= XMAX (Lemma 17a: odd(ord) coprime to
\\ odd(z), min(v2 ord, v2 z) <= 1) that can divide a BPSW psp: all q in one 2-adic class (v2 ord, v2 z) and the odd
\\ parts cross-coprime (lcm of odd ords coprime to lcm of odd z's). Output "D S A B" with A = lcm ord_q 2,
\\ B = lcm z_q (FULL orders; (W2) of SPEC.md: n == 1 mod A, n == -1 mod B). S = 1 is listed as "D 1 1 1".
\\ Input: ROWS = [[D, q, ord, z], ...] from the Theorem B scans (src/us/rows_1e11.gp, rows_2p40_d5.gp).
\\ Usage: gp -q splitparts.gp  (reads ../us/rows_1e11.gp and ../us/rows_2p40_d5.gp; writes splitparts_<XMAX>.txt)
oddp(x) = x >> valuation(x, 2);
adm(r) = min(valuation(r[3], 2), valuation(r[4], 2)) <= 1 && gcd(oddp(r[3]), oddp(r[4])) == 1;
parts(D, X) = {my(R = select(r -> r[1] == D && r[2] <= X && adm(r), ROWS), out = List([[1, 1, 1]]),
    cls = Set(apply(r -> [valuation(r[3], 2), valuation(r[4], 2)], R)));
  foreach(cls, c, my(L = select(r -> [valuation(r[3], 2), valuation(r[4], 2)] == c, R), k = #L, stack = List([[0, 1, 1, 1, 1, 1]]));
    while(#stack, my(t = stack[#stack]); listpop(stack);
      for(i = t[1] + 1, k, my(q = L[i][2], S = t[2]*q, Ao = lcm(t[3], oddp(L[i][3])), Bo = lcm(t[4], oddp(L[i][4])),
          A = lcm(t[5], L[i][3]), B = lcm(t[6], L[i][4]));
        if(gcd(Ao, Bo) > 1, next);
        \\ a split prime may not divide another member's order on the opposite side (n == 0 mod q vs n == +-1)
        if(gcd(S, A*B) > 1, next);
        listput(stack, [i, S, Ao, Bo, A, B]); listput(out, [S, A, B]))));
  Vec(out)};
read("../us/rows_1e11.gp"); R9 = ROWS; read("../us/rows_2p40_d5.gp"); ROWS = concat(R9, ROWS);
XMAX = 10^10;
{
  my(f = Str("splitparts_1e10.txt"), tot = 0);
  foreach([5, -7, -11, 13, -15, 17, -19, -23, 29, -31], D,
    my(P = parts(D, XMAX)); tot += #P;
    foreach(P, s, write(f, D, " ", s[1], " ", s[2], " ", s[3]));
    print("D = ", D, ": ", #P, " split parts (incl. S = 1), max S = ", vecmax(apply(s -> s[1], P))));
  print("total ", tot);
}
