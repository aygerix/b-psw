\\ E2 doubly-primitive census, D = 5 (P = 1, Q = -1, Fibonacci).
\\ All primes p with ord_p(2) <= X and z(p) <= X, where z = Fibonacci entry point.
\\ Usage:  echo 'X=10000; OUT="dso_X1e4.txt"; read("dso_census.gp")' | gp -q -D parisizemax=3000000000
\\ Method (complete): A = prod_{d<=X} Phi_d(2), B = prod_{e<=X} Phi_e(alpha,beta) (homogenised
\\ cyclotomic of the Fibonacci sequence, = prod_{k|e} F_k^mu(e/k)).  A prime p divides A iff
\\ ord_p(2) <= X (p | Phi_d(2) means ord = d/p^j <= d), and p | B iff z(p) <= X.  So the prime
\\ support of G = gcd(A, B) is exactly the census.  G is split by gcd with each Phi_d(2) and factored.
\\ run with: gp -q -D parisizemax=3000000000 (setting it inside read() aborts the read)
tt() = getabstime();
T0 = tt();
A = vecprod(vector(X, d, polcyclo(d, 2)));
TA = tt(); print("A: ", (exponent(A)+1), " bits, ", round((TA-T0)/100.)/10., " s");
FF = vector(X, k, fibonacci(k));
Phif(e) = {my(num = 1, den = 1); fordiv(e, k, my(m = moebius(e/k)); if(m == 1, num *= FF[k], m == -1, den *= FF[k])); num/den};
B = vecprod(vector(X, e, Phif(e)));
TB = tt(); print("B: ", (exponent(B)+1), " bits, ", round((TB-TA)/100.)/10., " s");
G = gcd(A, B);
TG = tt(); print("G = gcd(A,B): ", (exponent(G)+1), " bits, ", round((TG-TB)/100.)/10., " s");
A = 0; B = 0;
\\ split G by ord: g_d = gcd(R, Phi_d(2)); factor g_d.
LX = lcm(vector(X, i, i));
fibmod(k, p) = lift((Mod([1,1;1,0], p)^k)[1,2]);
\\ entry point z(p) given a multiple m of it whose prime factors are <= X
zdesc(p, m) = {my(f = factor(m)); for(i = 1, #f~, my(r = f[i,1]); for(j = 1, f[i,2], if(fibmod(m/r, p) == 0, m /= r, break))); m};
PR = Map(); R = G; nbig = 0;
{for(d = 1, X,
  my(g = gcd(R, polcyclo(d, 2)));
  if(g > 1,
    my(f = factor(g)[,1]);
    for(i = 1, #f, my(p = f[i]);
      if(!isprime(p), error("non-prime factor ", p));
      if(!mapisdefined(PR, p), mapput(PR, p, d));
      while(R % p == 0, R /= p))));
}
if(R != 1, error("G not exhausted, cofactor ", R));
TF = tt(); print("split+factor: ", #PR, " primes, ", round((TF-TG)/100.)/10., " s");
rows = List();
{foreach(Vec(PR), p,
  if(p == 2, next);
  my(o, z, kr = kronecker(5, p));
  o = if(p == 3, 2, znorder(Mod(2, p), mapget(PR, p)));   \\ PR(p) = d with p | Phi_d(2), so ord | d
  if(p == 5, z = 5, z = zdesc(p, gcd(p - kr, LX)));
  if(fibmod(z, p) != 0 || o > X || z > X, error("order check failed at ", p));
  listput(rows, [p, o, z, kr]));
}
rows = vecsort(Vec(rows), 1);
TZ = tt(); print("orders: ", round((TZ-TF)/100.)/10., " s; rows ", #rows, "; max p ", rows[#rows][1]);
write(OUT, rows);
print("done total ", round((TZ-T0)/100.)/10., " s");
