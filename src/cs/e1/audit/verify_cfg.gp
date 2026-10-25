\\ Independent PARI verification of a primal configuration (lower bound on N_max(b)).  Input: D, j1, j2, E1, E2 as
\\ [[q, e], ...] and a list of primes.  For each p: isprime, (D/p) != 0, ord_p 2 = 2^j1 * o1 with o1 | E1 and z_D(p) =
\\ 2^j2 * o2 with o2 | E2 (z from 2x2 matrix powers, independent of all line code).  bits = log2 phi(E1 E2).
Uk(k, Q, p) = lift((Mod([1, -Q; 1, 0], p)^k)[2, 1]);
zrank(p, Q, m) = { my(z = m, f = factor(m)[, 1]); if (Uk(m, Q, p) != 0, error("U_m")); for (i = 1, #f, my(l = f[i]); while (z % l == 0 && Uk(z / l, Q, p) == 0, z /= l)); z; }
vcfg(D, j1, j2, E1, E2, L) =
{
  my(Q = (1 - D) / 4, e1 = prod(i = 1, #E1, E1[i][1]^E1[i][2]), e2 = prod(i = 1, #E2, E2[i][1]^E2[i][2]), ok = 0, bad = 0, ninert = 0);
  if (gcd(e1, e2) != 1, error("rows on both sides"));
  for (i = 1, #L, my(p = L[i], o, z, js);
    if (!isprime(p), bad++; next);
    js = kronecker(D, p); if (js == 0, bad++; next);
    o = znorder(Mod(2, p)); z = zrank(p, Q, p - js);
    if (valuation(o, 2) == j1 && valuation(z, 2) == j2 && e1 % (o >> j1) == 0 && e2 % (z >> j2) == 0, ok++; ninert += (js == -1), bad++));
  print("D = ", D, " class (", j1, ",", j2, "): primes listed ", #L, ", verified usable ", ok, " (inert ", ninert, "), failed ", bad,
        "; bits = log2 phi(E1 E2) = ", log(eulerphi(e1 * e2)) / log(2));
}
