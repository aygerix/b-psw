\\ exact checks of the four usable forced primes named in the E4 result, and a numerical check of the reciprocity remark
fibm(n, p) = lift((Mod([1,1;1,0], p)^n)[1,2]);
zent(p) = my(N = p + 1, fa = factor(N), n = N); for (i = 1, #fa~, my(q = fa[i,1]); for (j = 1, fa[i,2], if (fibm(n/q, p) == 0, n /= q, break))); n;
ph(f) = my(c = 1); fordiv(f, d, c *= fibonacci(d)^moebius(f/d)); c;
show(p) = {
  my(o = znorder(Mod(2, p)), z = zent(p));
  print("p=", p, " bits ", #binary(p), " isprime(proven) ", isprime(p), " (5/p) ", kronecker(5, p), " p mod 8 ", p % 8);
  print("   p-1 = ", factor(p - 1), "   p+1 = ", factor(p + 1));
  print("   ord2 = ", o, " = ", factor(o), "   z = ", z, " = ", factor(z), "   class (v2 ord2, v2 z) = (", valuation(o,2), ",", valuation(z,2), ")");
  print("   odd(ord2) primes mod 4: ", apply(q -> q % 4, factor(o >> valuation(o,2))[,1]~), "  odd(z) primes mod 4: ", apply(q -> q % 4, factor(z >> valuation(z,2))[,1]~));
}
show(2^107 - 1); show(2^127 - 1);
show(18104700793); print("   divides Phi*_167: ", ph(167) % 18104700793 == 0);
show(486077598803); print("   divides Phi*_1596: ", ph(1596) % 486077598803 == 0);
\\ reciprocity remark: p inert, p | Phi*_f, odd q | p-1. In K = Q(sqrt5, zeta_m, zeta_q) the prime P0 of L = Q(sqrt5, zeta_m)
\\ (residue field F_{p^2}) with phi == zeta_m mod P0 splits into q-1 primes P_a (zeta_q -> w^a, w a fixed q-th root of 1 in F_p).
\\ (2/P_a)_q = zeta_q^{c_a} where w^{a c_a} = 2^{(p^2-1)/q}. Print s (2^{(p^2-1)/q} = w^s), and sum_a c_a mod q.
rec(p, q) = {
  my(g = znprimroot(p), w = g^((p - 1)/q), t = Mod(2, p)^((p^2 - 1)/q), s = -1, S = 0);
  for (k = 0, q - 1, if (w^k == t, s = k; break));
  for (a = 1, q - 1, S += s * lift(Mod(a, q)^-1));
  print("p=", p, " q=", q, ": 2 is a q-th power mod p? ", ispower(Mod(2,p), q), "; symbol exponent s on P_1 = ", s, "; exponent of the product over the ", q - 1, " conjugates = ", S % q);
}
rec(13, 3);
rec(18104700793, 23); rec(18104700793, 163); rec(18104700793, 311);
