\\ R2 reproduction: which Grantham-Alford primes (OEIS A018188) are usable for a given instance (L1, L2)?
\\ Usable (instance proof, src/pom/instance.py): p == 3 mod 8, (5/p) = -1, 2^L1 == -1 mod p (v2(ord) = 1, odd(ord) | L1),
\\ and gamma^(2 L2) == -1, i.e. alpha^(4 L2) == -1 in Z[phi]/p (v2 z = 2, odd(z) | L2). Input via read() of a gp file
\\ that sets L1, L2, inst (instance primes) and ga (A018188 primes).
fit(p) = { my(a);
  if (p % 8 != 3 || kronecker(5, p) != -1, return(0));
  if (Mod(2, p)^L1 != Mod(-1, p), return(0));
  a = Mod(Mod(1, p) * x, Mod(1, p) * (x^2 - x - 1))^(4 * L2);
  a == Mod(Mod(-1, p), Mod(1, p) * (x^2 - x - 1)) };
