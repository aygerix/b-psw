\\ Calibration of the Kummer null on random inert primes (no forcing, no residue-class conditioning).
\\ For odd q in Q: among inert p with q | p+1 (Lucas torus, element alpha/beta) count q | (p+1)/z(p);
\\ among inert p with q | p-1 (Fermat side, element 2) count q | (p-1)/ord_p(2). Kummer predicts 1/q in both.
fibm(n, p) = lift((Mod([1,1;1,0], p)^n)[1,2]);
setrand(777);
Q = [3, 7, 11, 13, 17, 19, 23];
nL = vector(#Q); hL = vector(#Q); nF = vector(#Q); hF = vector(#Q);
{
  my(cnt = 0);
  while (cnt < NS,
    my(p = randomprime([2^39, 2^40]));
    if (kronecker(5, p) != -1, next);
    cnt++;
    for (i = 1, #Q, my(q = Q[i]);
      if ((p + 1) % q == 0, nL[i]++; if (fibm((p + 1) / q, p) == 0, hL[i]++));
      if ((p - 1) % q == 0, nF[i]++; if (Mod(2, p)^((p - 1) / q) == 1, hF[i]++))));
  for (i = 1, #Q, my(q = Q[i]);
    printf("q=%2d  Lucas: %6d with q|p+1, q|index %6d, rate %.4f (1/q = %.4f, z %+.2f)   Fermat: %6d, %6d, rate %.4f (z %+.2f)\n",
      q, nL[i], hL[i], hL[i] / nL[i], 1. / q, (hL[i] - nL[i] / q) / sqrt(nL[i] / q * (1 - 1. / q)),
      nF[i], hF[i], hF[i] / nF[i], (hF[i] - nF[i] / q) / sqrt(nF[i] / q * (1 - 1. / q))));
}
quit;
